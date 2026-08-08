/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    NatureDSP_Baseband library. FIR filters and Related Functions
    Interpolating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

#include "firinterp_common.h"

/*-------------------------------------------------------------------------
Interpolating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with interpolation using real
IR stored in vector h. The complex data input is stored in vector x. The
filter output result is stored in vector y. The filter calculates N*D complex
output samples using M*D coefficients and requires last N+M-1 samples in the
delay line.

Representation:
firinterp   16-bit signed fixed-point format
            Filter coefficients are Q15
            Number of fractional bits for input/output samples is user-difined
firinterpf  IEEE-754 Std. single precision floating-point format for filter 
            coefficients and input/output samples

Parameters:
Input:
D           Interpolation ratio 
N           Length of input sample block
M           Length of subfilter. Total length of filter is M*D
h[M*D]      Filter coefficients; h[0] is to be multiplied by the newest 
            sample,Q15
x[N]        Input complex samples
Output:
y[N*D]      Output complex samples

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
N           Multiple of 8 (firinterp) or 4 (firinterpf)
M           2,4,8 or a positive multiple of 16 for D=2,3,4,6,12; or 
            a positive multiple of 8 for other D
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
subfilter lengths M=2,4,8,16 and 32 and interpolation factors D=2,3 and 4,
in any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firinterp[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "NatureDSP_Baseband_fir.h"
#include "common.h"

#define ALIGNED_CIRCULAR_WRAP(ptr, begin, end)  \
  assert(ptr<=(end));                           \
  if (ptr==(end)) ptr = (begin);

/* Filter processing generic function for all D*/
void filter_proc_dx_mx_8n (  
                               firinterp_handle_t handle,
                                 int16_t *  restrict  y,
                           const int16_t *  restrict  x,
                           const int16_t *  restrict  coef,
                                 int16_t *  restrict  delayLine,
                                      int   M,
                                      int   N,
                                      int   D
                          )
{
  const int exp = XT_NSA(D-1) - 16; // scale factor: ceil(log2(D))
  int d,n,m;
  int32_t rnd;
  const xb_vecNx16 *  restrict pX  = (const xb_vecNx16 *) x;
  const xb_vecNx16 *  restrict pH  = (const xb_vecNx16 *) coef;
        xb_vecNx16 *  restrict pD  = (      xb_vecNx16 *) delayLine;
        xb_vecNx16 *  restrict pD_ = (      xb_vecNx16 *) (delayLine + 2*(M+(-M&(BBE_SIMD_WIDTH/2-1))+BBE_SIMD_WIDTH));
        int32_t    *  restrict pYW = (      int32_t *) y;

  xb_vecNx16 y0;
  xb_vecNx16 d0, d1;
  uint32_t   c0, c1, c2, c3, c4, c5, c6, c7;
  xb_vecNx40 A0;
  vsaN       shft;
  valign     h_align;
  xb_vecNx16 CoefVec;
  NASSERT(N>0 && N%8==0);
  NASSERT(M%8==0);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(coef);
  NASSERT_ALIGN32(delayLine);
  
  shft = BBE_MOVVSA32(exp);  
  rnd = 1L<<(exp-1);
  
  __Pragma("ymemory(pH)")
  __Pragma("ymemory(pD)")
  __Pragma("ymemory(pD_)")
  for (d=0;d<D;d++)
  {
    pX  = (const xb_vecNx16 *) x;
    
    pD  = (      xb_vecNx16 *) delayLine;
    pD_ = (      xb_vecNx16 *) (delayLine + 2*(M+(-M&(BBE_SIMD_WIDTH/2-1))+BBE_SIMD_WIDTH));
    __Pragma("loop_count min=1");
    for (m=0;m<(2*M/(BBE_SIMD_WIDTH));m++)
    {
      BBE_LVNX16_IP(d0,pD_,2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(d0,pD,2*BBE_SIMD_WIDTH);
    }
    
    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
      xb_vecNx16 p00, p01, p02, p03, p04, p05, p06, p07;
      pD = (xb_vecNx16 *)delayLine;
      BBE_LVNX16_IP(d1, pX, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_X(d1, pD, 2 * 2 * M);
      pH = (const xb_vecNx16 *)(coef + M*d);
      h_align = BBE_LAVNX16_PP(pH);
     // A0 = BBE_MOVWA32(rnd);
      d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
      d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);

      BBE_LAVNX16_XP(CoefVec, h_align, pH, BBE_SIMD_WIDTH);
      c0 = BBE_EXTRNX16C(CoefVec, 0);
      c1 = BBE_EXTRNX16C(CoefVec, 1);
      c2 = BBE_EXTRNX16C(CoefVec, 2);
      c3 = BBE_EXTRNX16C(CoefVec, 3);

      BBE_SELPCNX16I(p01, p00, d1, d0, 1);
      BBE_SELPCNX16I(p03, p02, d1, d0, 3);
      BBE_SELPCNX16I(p05, p04, d1, d0, 5);
      BBE_SELPCNX16I(p07, p06, d1, d0, 7);

      BBE_SVNX16_IP(d1, pD, 1 * 2 * BBE_SIMD_WIDTH);

      A0 = BBE_MULRNX16PR(p07, p06, c3,rnd);
      BBE_MULANX16PR(A0, p05, p04, c2);
      BBE_MULANX16PR(A0, p03, p02, c1);
      BBE_MULANX16PR(A0, p01, p00, c0);

      for (m = 0; m<(M / (BBE_SIMD_WIDTH / 2) - 1); m++)
      {
        d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
        d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);

        BBE_LAVNX16_XP(CoefVec, h_align, pH, BBE_SIMD_WIDTH);
        c0 = BBE_EXTRNX16C(CoefVec, 0);
        c1 = BBE_EXTRNX16C(CoefVec, 1);
        c2 = BBE_EXTRNX16C(CoefVec, 2);
        c3 = BBE_EXTRNX16C(CoefVec, 3);

        BBE_SELPCNX16I(p01, p00, d1, d0, 1);
        BBE_SELPCNX16I(p03, p02, d1, d0, 3);
        BBE_SELPCNX16I(p05, p04, d1, d0, 5);
        BBE_SELPCNX16I(p07, p06, d1, d0, 7);

        BBE_SVNX16_IP(d1, pD, 1 * 2 * BBE_SIMD_WIDTH);

        BBE_MULANX16PR(A0, p07, p06, c3);
        BBE_MULANX16PR(A0, p05, p04, c2);
        BBE_MULANX16PR(A0, p03, p02, c1);
        BBE_MULANX16PR(A0, p01, p00, c0);
      }
      y0 = BBE_PACKVNX40(A0, shft);
      c0 = BBE_EXTRNX16C(y0, 0);
      c1 = BBE_EXTRNX16C(y0, 1);
      c2 = BBE_EXTRNX16C(y0, 2);
      c3 = BBE_EXTRNX16C(y0, 3);
      c4 = BBE_EXTRNX16C(y0, 4);
      c5 = BBE_EXTRNX16C(y0, 5);
      c6 = BBE_EXTRNX16C(y0, 6);
      c7 = BBE_EXTRNX16C(y0, 7);
      pYW[8 * n*D + 0 * D + d] = (int32_t)c0;
      pYW[8 * n*D + 1 * D + d] = (int32_t)c1;
      pYW[8 * n*D + 2 * D + d] = (int32_t)c2;
      pYW[8 * n*D + 3 * D + d] = (int32_t)c3;
      pYW[8 * n*D + 4 * D + d] = (int32_t)c4;
      pYW[8 * n*D + 5 * D + d] = (int32_t)c5;
      pYW[8 * n*D + 6 * D + d] = (int32_t)c6;
      pYW[8 * n*D + 7 * D + d] = (int32_t)c7;
    }
  }
  
  pD  = (xb_vecNx16 *) delayLine;
  pD_ = (xb_vecNx16 *)(delayLine + 2 * (M + (-M&(BBE_SIMD_WIDTH / 2 - 1)) + BBE_SIMD_WIDTH));
  __Pragma("loop_count min=1")
  for (m=0;m<(2*M/(BBE_SIMD_WIDTH));m++)
  {
    BBE_LVNX16_IP(d0,pD, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(d0,pD_,2*BBE_SIMD_WIDTH);
  }
}
const tFirFxdxns interp_dx_mx_8n  = {&firinterp_gen,filter_proc_dx_mx_8n};

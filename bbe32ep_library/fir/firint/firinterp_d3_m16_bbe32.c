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

/* Filter processing function for D=3 M=16 and N=8*n */
static void filter_proc_3d_16_8n (         
                              void* handle,
                                 int16_t *  restrict  y,
                           const int16_t *  restrict  x,
                           const int16_t *  restrict  coef,
                                 int16_t *  restrict  delayLine,
                                      int   M,
                                      int   N,
                                      int   D
                          )
{
  int n;
  const xb_vecNx16 *  restrict pX = (const xb_vecNx16 *)x;
  const        int *  restrict pH = (const int        *)coef;
  xb_vecNx16 *  restrict pD = (xb_vecNx16 *)delayLine;
  xb_vecNx16 *  restrict pY = (xb_vecNx16 *)y;

  xb_vecNx16 x0, x1, y0, y1, y2;
  xb_vecNx16 d0, d1, d2;
  xb_vecNx16 z0, z1, z2;
  uint32_t   c00, c01, c10, c11, c20, c21;
  uint32_t   c02, c03, c12, c13, c22, c23;
  uint32_t   c04, c05, c14, c15, c24, c25;
  uint32_t   c06, c07, c16, c17, c26, c27;
  xb_vecNx40 A0, A1, A2;
  vsaN       shft = BBE_MOVVSA32(13);
  NASSERT(N % 8 == 0);
  NASSERT(M == 16);
  NASSERT(D == 3);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(coef);
  NASSERT_ALIGN32(delayLine);

  if (N<=0) return;
  __Pragma("loop_count min=1");
  for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
  {
    xb_vecNx16 p00, p01, p02, p03, p04, p05, p06, p07;
    xb_vecNx16 p10, p11, p12, p13, p14, p15, p16, p17;
    unsigned tmp;

        
    BBE_L32IP(tmp, pH, 4);  c00 = tmp;
    BBE_L32IP(tmp, pH, 4);  c10 = tmp;
    BBE_L32IP(tmp, pH, 4);  c20 = tmp;
    BBE_L32IP(tmp, pH, 4);  c01 = tmp;
    BBE_L32IP(tmp, pH, 4);  c11 = tmp;
    BBE_L32IP(tmp, pH, 4);  c21 = tmp;
    BBE_L32IP(tmp, pH, 4);  c02 = tmp;
    BBE_L32IP(tmp, pH, 4);  c12 = tmp;
    BBE_L32IP(tmp, pH, 4);  c22 = tmp;
    BBE_L32IP(tmp, pH, 4);  c03 = tmp;
    BBE_L32IP(tmp, pH, 4);  c13 = tmp;
    BBE_L32IP(tmp, pH, 4);  c23 = tmp;
    BBE_L32IP(tmp, pH, 4);  c04 = tmp;
    BBE_L32IP(tmp, pH, 4);  c14 = tmp;
    BBE_L32IP(tmp, pH, 4);  c24 = tmp;
    BBE_L32IP(tmp, pH, 4);  c05 = tmp;
    BBE_L32IP(tmp, pH, 4);  c15 = tmp;
    BBE_L32IP(tmp, pH, 4);  c25 = tmp;


    c06 = XT_L32I(pH, 0);
    c16 = XT_L32I(pH, 4);
    c26 = XT_L32I(pH, 8);
    c07 = XT_L32I(pH, 12);
    c17 = XT_L32I(pH, 16);
    c27 = XT_L32I(pH, 20);

    d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
    d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(d2, pX, 2 * BBE_SIMD_WIDTH);

    pH-=18;

    BBE_SELPCNX16I(p01, p00, d1, d0, 1);
    BBE_SELPCNX16I(p03, p02, d1, d0, 3);
    BBE_SELPCNX16I(p05, p04, d1, d0, 5);
    p07 = d1;
    p06 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_14);

    BBE_SELPCNX16I(p11, p10, d2, d1, 1);
    BBE_SELPCNX16I(p13, p12, d2, d1, 3);
    BBE_SELPCNX16I(p15, p14, d2, d1, 5);
    p17 = d2;
    p16 = BBE_SELNX16I(d2, d1, BBE_SELI_ROTATE_RIGHT_14);
    BBE_SVNX16_I(d1, pD, 0 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(d2, pD, 1 * 2 * BBE_SIMD_WIDTH);

    A0 = BBE_MULNX16PR(p16, p17, c00);
    A1 = BBE_MULNX16PR(p16, p17, c10);
    A2 = BBE_MULNX16PR(p16, p17, c20);
        
    BBE_MULANX16PR(A0, p14, p15, c01);
    BBE_MULANX16PR(A1, p14, p15, c11);
    BBE_MULANX16PR(A2, p14, p15, c21);
    BBE_MULANX16PR(A0, p12, p13, c02);
    BBE_MULANX16PR(A1, p12, p13, c12);
    BBE_MULANX16PR(A2, p12, p13, c22);
    BBE_MULANX16PR(A0, p10, p11, c03);
    BBE_MULANX16PR(A1, p10, p11, c13);
    BBE_MULANX16PR(A2, p10, p11, c23);
        

    BBE_MULANX16PR(A0, p06, p07, c04);
    BBE_MULANX16PR(A1, p06, p07, c14);
    BBE_MULANX16PR(A2, p06, p07, c24);
    BBE_MULANX16PR(A0, p04, p05, c05);
    BBE_MULANX16PR(A1, p04, p05, c15);
    BBE_MULANX16PR(A2, p04, p05, c25);
    BBE_MULANX16PR(A0, p02, p03, c06);
    BBE_MULANX16PR(A1, p02, p03, c16);
    BBE_MULANX16PR(A2, p02, p03, c26);
    BBE_MULANX16PR(A0, p00, p01, c07);
    BBE_MULANX16PR(A1, p00, p01, c17);
    BBE_MULANX16PR(A2, p00, p01, c27);

    y0 = BBE_PACKVNX40(A0, shft);
    y1 = BBE_PACKVNX40(A1, shft);
    y2 = BBE_PACKVNX40(A2, shft);
        

    BBE_DSELNX16I(x1, x0, y1, y0, BBE_DSELI_INTERLEAVE_2);
        
    BBE_DSELNX16I(z1, z0, y2, x0, BBE_DSELI_INTERLEAVE_C3_STEP_0);

    BBE_DSELNX16I_H(z1, z2, y2, x1, BBE_DSELI_INTERLEAVE_C3_STEP_1);

    BBE_SVNX16_IP(z0, pY, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z1, pY, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z2, pY, 2 * BBE_SIMD_WIDTH);
  }
}
const tFirFxdxns interp_3d_16_8n   ={&firinterp_dx, filter_proc_3d_16_8n };

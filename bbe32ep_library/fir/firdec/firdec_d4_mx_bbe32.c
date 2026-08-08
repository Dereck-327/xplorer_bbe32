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
    Decimating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "firdec_common.h"

/*-------------------------------------------------------------------------
Decimating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with decimation using real IR 
stored in vector h. The complex data input is stored in vector x. The filter
output result is stored in vector y. The filter calculates N output samples
using M coefficients and requires last D*N+M-1 samples in the delay line.

NOTE:
To avoid aliasing, the IR should be synthesized in such a way that filter pass
band is limited by input sample frequency divided by 2*D.

Representation:
firdec   16-bit signed fixed-point format
         Filter coefficients are Q15
         Number of fractional bits for input/output samples is user-difined
firdecf  IEEE-754 Std. single precision floating-point format for filter 
         coefficients and input/output samples

Parameters:
Input:
D        Decimation factor
N        Length of output sample block
M        Length of filter
h[M]     Filter coefficients; h[0] is to be multiplied by the newest 
         sample
x[N*D]   Input complex samples
Output:
y[N]     Output complex samples

Restrictions:
x,y      Must not overlap
x,y      Aligned on 32-byte boundary
N        Multiple of 8 (firdec) or 4 (firdecf)
M        2,4,8 or a positive multiple of 16 for D=2,3,4; or 
         a positive multiple of 16 for D>4
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
filter lengths M=2,4,8,16 and 32 and decimation factors D=2,3 and 4, in
any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firdec[f]_init returns NULL handle.
-------------------------------------------------------------------------*/


static void firdec_proc_D4_MXodd(int16_t * restrict y,
  const int16_t * restrict x,
  const int16_t * restrict coef,
  int16_t * restrict delayLine,
  int M, int N, int D)
{
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 * restrict pD;
  xb_vecNx16 * restrict pD1;
  xb_vecNx16 * restrict pH;
  const xb_vecNx16 *          X;

  xb_vecNx40 w0, w1;

  xb_vecNx16 cf;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 y0, y1, y2, y3;
  xb_vecNx16 y4, y5, y6, y7;
  xb_vecNx16 y8, y9, y10, y11;

  xb_vecNx16 p08, p09, p0a, p0b;
  xb_vecNx16 p18, p19, p1a, p1b;
  xb_vecNx16 p28, p29, p2a, p2b;
  xb_vecNx16 p38, p39, p3a, p3b;

  uint32_t   q00, q01;
  uint32_t   q10, q11;
  uint32_t   q20, q21;
  uint32_t   q30, q31;

  int n, m;

  NASSERT(N>0 && !(N & 7));

  NASSERT(D == 4 && !(M & 15));

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(coef);
  NASSERT_ALIGN32(delayLine);

  Y = (xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;
  for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
  {
    w0 = 0;
    w1 = 0;
    pH = (xb_vecNx16*)coef;
    cf = BBE_LVNX16_I(pH, 0);
    pD = (const xb_vecNx16 *)(delayLine);
    pD1 = (xb_vecNx16 *)(delayLine);

    // Load 16x4 input samples, CQ15
    BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, X, 2 * BBE_SIMD_WIDTH);
    // deinterleave inputs
    BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);

    BBE_SVNX16_X(x0, pD1, 2 * 2 * (M + 16));
    BBE_SVNX16_X(x1, pD1, 2 * 2 * (M + 16) + 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_X(x2, pD1, 2 * 2 * (M + 16) + 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_X(x3, pD1, 2 * 2 * (M + 16) + 3 * 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(y0, pD, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, pD, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, pD, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y3, pD, 2 * BBE_SIMD_WIDTH);

    __Pragma("ymemory(pD)");
    __Pragma("ymemory(pH)");
    for (m = 0; m<(M / (BBE_SIMD_WIDTH * 2)); m++)
    {
      // Load 8x2 input samples, CQ15
      BBE_LVNX16_IP(y4, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y5, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y6, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y7, pD, -1 * 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);

      // coefficients bank 0
      q00 = BBE_EXTRNX16C(cf, 6);
      q01 = BBE_EXTRNX16C(cf, 7);
      // coefficients bank 3
      q10 = BBE_EXTRNX16C(cf, 4);
      q11 = BBE_EXTRNX16C(cf, 5);
      // coefficients bank 2
      q20 = BBE_EXTRNX16C(cf, 2);
      q21 = BBE_EXTRNX16C(cf, 3);
      // coefficients bank 1
      q30 = BBE_EXTRNX16C(cf, 0);
      q31 = BBE_EXTRNX16C(cf, 1);

      BBE_SELPCNX16I(p09, p08, y4, y0, 5);
      BBE_SELPCNX16I(p0b, p0a, y4, y0, 7);

      BBE_SELPCNX16I(p19, p18, y5, y1, 4);
      BBE_SELPCNX16I(p1b, p1a, y5, y1, 6);

      BBE_SELPCNX16I(p29, p28, y6, y2, 4);
      BBE_SELPCNX16I(p2b, p2a, y6, y2, 6);

      BBE_SELPCNX16I(p39, p38, y7, y3, 4);
      BBE_SELPCNX16I(p3b, p3a, y7, y3, 6);

      BBE_MULANX16PR(w0, p09, p08, q00);
      BBE_MULANX16PR(w1, p0b, p0a, q01);

      BBE_MULANX16PR(w0, p19, p18, q30);
      BBE_MULANX16PR(w1, p1b, p1a, q31);

      BBE_MULANX16PR(w0, p29, p28, q20);
      BBE_MULANX16PR(w1, p2b, p2a, q21);

      BBE_MULANX16PR(w0, p39, p38, q10);
      BBE_MULANX16PR(w1, p3b, p3a, q11);

      BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);

      // coefficients bank 0
      q00 = BBE_EXTRNX16C(cf, 6);
      q01 = BBE_EXTRNX16C(cf, 7);
      // coefficients bank 3
      q10 = BBE_EXTRNX16C(cf, 4);
      q11 = BBE_EXTRNX16C(cf, 5);
      // coefficients bank 2
      q20 = BBE_EXTRNX16C(cf, 2);
      q21 = BBE_EXTRNX16C(cf, 3);
      // coefficients bank 1
      q30 = BBE_EXTRNX16C(cf, 0);
      q31 = BBE_EXTRNX16C(cf, 1);


      y0 = y4;
      y1 = y5;
      BBE_LVNX16_IP(y2, pD, 2 * 2 * BBE_SIMD_WIDTH);
      y3 = BBE_LVNX16_I(pD, -2 * BBE_SIMD_WIDTH);

      y8 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
      y9 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
      y10 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH);
      y11 = BBE_LVNX16_I(pD, 3 * 2 * BBE_SIMD_WIDTH);

      BBE_SELPCNX16I(p09, p08, y8, y4, 1);
      BBE_SELPCNX16I(p0b, p0a, y8, y4, 3);

      BBE_SELPCNX16I(p19, p18, y9, y5, 0);
      BBE_SELPCNX16I(p1b, p1a, y9, y5, 2);

      BBE_SELPCNX16I(p29, p28, y10, y6, 0);
      BBE_SELPCNX16I(p2b, p2a, y10, y6, 2);

      BBE_SELPCNX16I(p39, p38, y11, y7, 0);
      BBE_SELPCNX16I(p3b, p3a, y11, y7, 2);

      BBE_MULANX16PR(w0, p09, p08, q00);
      BBE_MULANX16PR(w1, p0b, p0a, q01);

      BBE_MULANX16PR(w0, p19, p18, q30);
      BBE_MULANX16PR(w1, p1b, p1a, q31);

      BBE_MULANX16PR(w0, p29, p28, q20);
      BBE_MULANX16PR(w1, p2b, p2a, q21);

      BBE_MULANX16PR(w0, p39, p38, q10);
      BBE_MULANX16PR(w1, p3b, p3a, q11);

      // update delay line 
      BBE_SVNX16_IP(y0, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y1, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y2, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y3, pD1, 2 * BBE_SIMD_WIDTH);
    }

    {
      // Load 8x2 input samples, CQ15
      BBE_LVNX16_IP(y4, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y5, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y6, pD, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(y7, pD, -1 * 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);

      // coefficients bank 0
      q00 = BBE_EXTRNX16C(cf, 6);
      q01 = BBE_EXTRNX16C(cf, 7);
      // coefficients bank 3
      q10 = BBE_EXTRNX16C(cf, 4);
      q11 = BBE_EXTRNX16C(cf, 5);
      // coefficients bank 2
      q20 = BBE_EXTRNX16C(cf, 2);
      q21 = BBE_EXTRNX16C(cf, 3);
      // coefficients bank 1
      q30 = BBE_EXTRNX16C(cf, 0);
      q31 = BBE_EXTRNX16C(cf, 1);


      BBE_SELPCNX16I(p09, p08, y4, y0, 5);
      BBE_SELPCNX16I(p0b, p0a, y4, y0, 7);

      BBE_SELPCNX16I(p19, p18, y5, y1, 4);
      BBE_SELPCNX16I(p1b, p1a, y5, y1, 6);

      BBE_SELPCNX16I(p29, p28, y6, y2, 4);
      BBE_SELPCNX16I(p2b, p2a, y6, y2, 6);

      BBE_SELPCNX16I(p39, p38, y7, y3, 4);
      BBE_SELPCNX16I(p3b, p3a, y7, y3, 6);

      BBE_MULANX16PR(w0, p09, p08, q00);
      BBE_MULANX16PR(w1, p0b, p0a, q01);

      BBE_MULANX16PR(w0, p19, p18, q30);
      BBE_MULANX16PR(w1, p1b, p1a, q31);

      BBE_MULANX16PR(w0, p29, p28, q20);
      BBE_MULANX16PR(w1, p2b, p2a, q21);

      BBE_MULANX16PR(w0, p39, p38, q10);
      BBE_MULANX16PR(w1, p3b, p3a, q11);

      // update delay line
      BBE_SVNX16_IP(y4, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y5, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y6, pD1, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y7, pD1, 2 * BBE_SIMD_WIDTH);
    }

    //
    // Save 8 output samples.
    //
    w0 = BBE_ADDNX40(w0, w1);
    // CQ15 <- CQ30 - 15 w/ rounding and saturation.
    y0 = BBE_PACKQNX40(w0);

    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
  }

}

static void firdec_proc_D4_MXeven(int16_t * restrict y,
  const int16_t * restrict x,
  const int16_t * restrict coef,
  int16_t * restrict delayLine,
  int M, int N, int D)
{
  xb_vecNx16 * restrict Y;
  xb_vecNx16 * restrict pD;
  xb_vecNx16 * restrict pH;
  const xb_vecNx16 *          X;

  xb_vecNx40 w0, w1;

  xb_vecNx16 cf;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 y0, y1, y2, y3;
  xb_vecNx16 y4, y5, y6, y7;

  xb_vecNx16 p08, p09, p0a, p0b;
  xb_vecNx16 p18, p19, p1a, p1b;
  xb_vecNx16 p28, p29, p2a, p2b;
  xb_vecNx16 p38, p39, p3a, p3b;

  uint32_t   q00, q01;
  uint32_t   q10, q11;
  uint32_t   q20, q21;
  uint32_t   q30, q31;

  int n, m;

  NASSERT(N>0 && !(N & 7));

  NASSERT(D == 4 && !(M & 31));

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(coef);
  NASSERT_ALIGN32(delayLine);

  Y = (xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;

  for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
  {
    w0 = 0;
    w1 = 0;
    pH = (xb_vecNx16*)coef;
    cf = BBE_LVNX16_I(pH, 0);
    pD = (xb_vecNx16 *)(delayLine);

    // Load 8x4 input samples, CQ15
    BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, X, 2 * BBE_SIMD_WIDTH);
    // deinterleave inputs
    BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);

    BBE_SVNX16_X(x0, pD, 2 * 2 * M);
    BBE_SVNX16_X(x1, pD, 2 * 2 * M + 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_X(x2, pD, 2 * 2 * M + 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_X(x3, pD, 2 * 2 * M + 3 * 2 * BBE_SIMD_WIDTH);

    y0 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 0);
    y1 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 2 * BBE_SIMD_WIDTH);
    y2 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 4 * BBE_SIMD_WIDTH);
    y3 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 6 * BBE_SIMD_WIDTH);

    __Pragma("ymemory(pD)");
    __Pragma("ymemory(pH)");
    __Pragma("loop_count min=1");
    for (m = 0; m<(M / (BBE_SIMD_WIDTH * 2)); m++)
    {
      // Load 8x2 input samples, CQ15
      y4 = BBE_LVNX16_I(pD, 4 * 2 * BBE_SIMD_WIDTH);
      y5 = BBE_LVNX16_I(pD, 5 * 2 * BBE_SIMD_WIDTH);
      y6 = BBE_LVNX16_I(pD, 6 * 2 * BBE_SIMD_WIDTH);
      y7 = BBE_LVNX16_I(pD, 7 * 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);

      // coefficients bank 0
      q00 = BBE_EXTRNX16C(cf, 6);
      q01 = BBE_EXTRNX16C(cf, 7);
      // coefficients bank 3
      q10 = BBE_EXTRNX16C(cf, 4);
      q11 = BBE_EXTRNX16C(cf, 5);
      // coefficients bank 2
      q20 = BBE_EXTRNX16C(cf, 2);
      q21 = BBE_EXTRNX16C(cf, 3);
      // coefficients bank 1
      q30 = BBE_EXTRNX16C(cf, 0);
      q31 = BBE_EXTRNX16C(cf, 1);

      BBE_SELPCNX16I(p09, p08, y4, y0, 1);
      BBE_SELPCNX16I(p0b, p0a, y4, y0, 3);

      BBE_SELPCNX16I(p19, p18, y5, y1, 0);
      BBE_SELPCNX16I(p1b, p1a, y5, y1, 2);

      BBE_SELPCNX16I(p29, p28, y6, y2, 0);
      BBE_SELPCNX16I(p2b, p2a, y6, y2, 2);

      BBE_SELPCNX16I(p39, p38, y7, y3, 0);
      BBE_SELPCNX16I(p3b, p3a, y7, y3, 2);

      BBE_MULANX16PR(w0, p09, p08, q00);
      BBE_MULANX16PR(w1, p0b, p0a, q01);

      BBE_MULANX16PR(w0, p19, p18, q30);
      BBE_MULANX16PR(w1, p1b, p1a, q31);

      BBE_MULANX16PR(w0, p29, p28, q20);
      BBE_MULANX16PR(w1, p2b, p2a, q21);

      BBE_MULANX16PR(w0, p39, p38, q10);
      BBE_MULANX16PR(w1, p3b, p3a, q11);

      BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);

      // coefficients bank 0
      q00 = BBE_EXTRNX16C(cf, 6);
      q01 = BBE_EXTRNX16C(cf, 7);
      // coefficients bank 3
      q10 = BBE_EXTRNX16C(cf, 4);
      q11 = BBE_EXTRNX16C(cf, 5);
      // coefficients bank 2
      q20 = BBE_EXTRNX16C(cf, 2);
      q21 = BBE_EXTRNX16C(cf, 3);
      // coefficients bank 1
      q30 = BBE_EXTRNX16C(cf, 0);
      q31 = BBE_EXTRNX16C(cf, 1);

      BBE_SELPCNX16I(p09, p08, y4, y0, 5);
      BBE_SELPCNX16I(p0b, p0a, y4, y0, 7);

      BBE_SELPCNX16I(p19, p18, y5, y1, 4);
      BBE_SELPCNX16I(p1b, p1a, y5, y1, 6);

      BBE_SELPCNX16I(p29, p28, y6, y2, 4);
      BBE_SELPCNX16I(p2b, p2a, y6, y2, 6);

      BBE_SELPCNX16I(p39, p38, y7, y3, 4);
      BBE_SELPCNX16I(p3b, p3a, y7, y3, 6);

      BBE_MULANX16PR(w0, p09, p08, q00);
      BBE_MULANX16PR(w1, p0b, p0a, q01);

      BBE_MULANX16PR(w0, p19, p18, q30);
      BBE_MULANX16PR(w1, p1b, p1a, q31);

      BBE_MULANX16PR(w0, p29, p28, q20);
      BBE_MULANX16PR(w1, p2b, p2a, q21);

      BBE_MULANX16PR(w0, p39, p38, q10);
      BBE_MULANX16PR(w1, p3b, p3a, q11);

      // update delay line
      BBE_SVNX16_IP(y4, pD, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y5, pD, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y6, pD, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(y7, pD, 2 * BBE_SIMD_WIDTH);

      y0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
      y1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
      y2 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH);
      y3 = BBE_LVNX16_I(pD, 3 * 2 * BBE_SIMD_WIDTH);
    }
    //
    // Save 8 output samples.
    //
    w0 = BBE_ADDNX40(w0, w1);
    // CQ15 <- CQ30 - 15 w/ rounding and saturation.
    y0 = BBE_PACKQNX40(w0);

    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
  }

}

const tFirFxdxns firdec_4d_xeven_8n = { &firdec_alloc_d2_mx, firdec_proc_D4_MXeven };
const tFirFxdxns firdec_4d_xodd_8n = { &firdec_alloc_d2_mx, firdec_proc_D4_MXodd };


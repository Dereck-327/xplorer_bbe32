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
  NatureDSP_Baseband library. FFT
    Mixed radix inverse FFT on complex data, auto scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#include "fft_tw.h"
/*-------------------------------------------------------------------------
Mixed radix inverse FFT on complex data, auto scaling
  
Description: These functions make inverse FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120,144,180,192,216,240,288,300,324,360,
384,432,480,540,576,600,648,720,768,864,900,960,972,1080,1152,1200,1536.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of CINFFT_SCRATCH_SIZE(N) bytes
  Input:            
    S           Required input/output buffer size may exceed actual data size:
                S >= N. Use CINFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[S]        N complex samples of input spectrum
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[S]        N complex samples of output signal
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

int cinfft12 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp )
{
  NASSERT_ALIGN32(pScr);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  xb_vecNx16 w, x0, x1, tmp, y0, y1;
  xb_vecNx16 *pw = (xb_vecNx16 *)fft12_tw;
  xb_vecNx40 acc0;
  int scaling = 4 - bexp;
  vsaN s = BBE_MOVVSA32((int)15 + scaling);  // W12 in Q15 format

  x0 = BBE_LVNX16_I((xb_vecNx16 *)x, 0);
  x1 = BBE_LVNX16_I((xb_vecNx16 *)x, 2 * BBE_SIMD_WIDTH);

  tmp = BBE_REPNX16C(x0, 0);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); acc0 = BBE_MULRNX16J(tmp, w, s);
  tmp = BBE_REPNX16C(x0, 1);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 2);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 3);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  tmp = BBE_REPNX16C(x0, 4);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 5);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 6);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 7);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  tmp = BBE_REPNX16C(x1, 0);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 1);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 2);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 3);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  y0 = BBE_PACKVNX40(acc0, s);
  BBE_SVNX16_I(y0, (xb_vecNx16*)y, 0);

  tmp = BBE_REPNX16C(x0, 0);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); acc0 = BBE_MULRNX16J(tmp, w, s);
  tmp = BBE_REPNX16C(x0, 1);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 2);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 3);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  tmp = BBE_REPNX16C(x0, 4);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 5);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 6);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x0, 7);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  tmp = BBE_REPNX16C(x1, 0);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 1);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 2);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);
  tmp = BBE_REPNX16C(x1, 3);  BBE_LVNX16_IP(w, pw, BBE_SIMD_WIDTH * 2); BBE_MULANX16J(acc0, tmp, w);

  y1 = BBE_PACKVNX40(acc0, s);
  BBE_SVNX16_I(y1, (xb_vecNx16*)y, 2 * BBE_SIMD_WIDTH);

  return scaling; 

} /* cinfft12() */

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
#define IS_INV_FFT
#include "fft_common.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(int, cinfft96, (void * restrict pScr,complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp))
#else
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

int cinfft96 (void * restrict pScr,complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp)
{
  int scaling;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  scaling = R1_DFT4_L64_16(fft96_tw1, (int16_t*)x, (int16_t*)y, 96, bexp);
  scaling += R2_DFT4xI4(fft96_tw2, (int16_t*)y, (int16_t*)x, 96, BBE_RRANGE());
  scaling += R2_DFT6xIN_6((int16_t*)x, (int16_t*)y, 96, BBE_RRANGE());
  return scaling;
} /* cinfft96() */
#endif

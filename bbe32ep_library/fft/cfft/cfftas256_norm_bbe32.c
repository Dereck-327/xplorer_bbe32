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
    Radix-2 forward FFT on complex data, auto scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#include "fft_common.h"
#if !(HAVE_FFT && 1)
DISCARD_FUN(int, cfftas256_norm, (complex_fract16 * restrict y, complex_fract16 * restrict x))
#else
/*-------------------------------------------------------------------------
Radix-2 forward FFT on complex data, auto scaling

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=4..15. Functions with _norm suffix expect input data to be
normalized, i.e. the minimum number of redundant sign bits over x[]
(a.k.a the common block exponent) should be zero. Neglecting to normalize
data leads to significant loss in transform quality. On the contrary, regular
variants with no _norm suffix allow for non-zero common block exponent, but
they appear slightly slower due to internal data normalization.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Complex input signal
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Output spectrum samples
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data, with positive numbers corresponding to the right
                shift. _norm-suffixed functions return strictly positive
                values, while for regular variants the total shift amount is
                bi-directional.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

int cfftas256_norm ( complex_fract16 * restrict y, complex_fract16 * restrict x           )
{
  int scaling;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  scaling = R1_DFT4_L64_16_norm(fft256_tw1, (int16_t*)x, (int16_t*)y, 256);
  scaling += R2_DFT4xI4(fft256_tw2, (int16_t*)y, (int16_t*)x, 256, BBE_RRANGE());
  scaling += R3_DFT4xIv(fft256_tw3, (int16_t*)x, (int16_t*)y, 256, 16, BBE_RRANGE());
  scaling += R2_DFT4xIN_4((int16_t*)y, (int16_t*)y, 256, BBE_RRANGE());

  return scaling;
} /* cfftas256_norm() */
#endif

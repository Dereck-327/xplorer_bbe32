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
    Radix-2 forward FFT on complex data with reduced twiddle table, auto scaling
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
DISCARD_FUN(int, tfft8192, (complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp))
#else

/*-------------------------------------------------------------------------
Radix-2 forward FFT on complex data with reduced twiddle table, auto scaling

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Complex input signal
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]      Output spectrum samples
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

int tfft8192 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp )
{
  int scaling;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  scaling = R1_tDFT4_L64_16(tfft8192_tw1,(int16_t*)x, (int16_t*)y, 8192, bexp);
  scaling += R2_tDFT4xI4_U2(tfft8192_tw2,(int16_t*)y, (int16_t*)x, 8192, BBE_RRANGE());
  scaling += R3_DFT4xIv(tfft8192_tw3, (int16_t*)x, (int16_t*)y, 8192, 16, BBE_RRANGE());
  scaling += R3_DFT4xIv(tfft8192_tw4, (int16_t*)y, (int16_t*)x, 8192, 64, BBE_RRANGE());
  scaling += R3_DFT4xIv_2U(tfft8192_tw5, (int16_t*)x, (int16_t*)y, 8192, 256, BBE_RRANGE());
  scaling += R2_DFT8xIN_8((int16_t*)y, (int16_t*)y, 8192, BBE_RRANGE());

  return scaling;
} /* tfft8192() */
#endif

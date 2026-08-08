/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
NatureDSP_Baseband library. FFT
Radix-2 forward FFT on complex_float data
C code optimized for BBE32
IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Tables of a twiddles */
#include "fft_fp_tw.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"

/*-------------------------------------------------------------------------
Radix-2 inverse floating point FFT forming real data

Description: These functions make inverse FFT forming real data of length
N=2^n, n=4..15. Algorithm exploits the symmetry properties of the FFT:
the input spectrum is modified in such a way that a complex-valued inverse FFT
of half the original size that is applied to the transformed spectrum actually
results in real data.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:
    x[(N/2+1)]    Input spectrum samples
  Output:
    y[N]          Real output signal
  Returned value:
                  zero
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(int, rifft32f, (float32_t * restrict y, complex_float * restrict x) )
#else
int rifft32f(float32_t * restrict y, complex_float * restrict x)
{
    const complex_float *twiddle_table = (const complex_float*)rfftf32_twd1;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT(x != (complex_float*)y);

    rifft_spec_conv_fp(x, twiddle_table, 32);
    cifft16f((complex_float*)y, x);

     return 0;
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU

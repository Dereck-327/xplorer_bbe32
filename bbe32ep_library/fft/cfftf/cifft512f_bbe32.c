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


#define IS_INV_FFTF 1
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"

/*-------------------------------------------------------------------------
Radix-2 floating point inverse FFT on complex data

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=4..15. 

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Input spectrum samples
  Output:            
    y[N]        Complex output signal
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(int, cifft512f, (complex_float * restrict y, complex_float * restrict x) )
#else
int cifft512f(complex_float * restrict y, complex_float * restrict x)
{
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT(x != y);

    first_stage_iDFT4_FP_v2((const complex_float*)cfftf512_twd1, x, y, 512);
    stage_inner_R3_iDFT4xIv_FP((const complex_float*)cfftf512_twd2, y, x, 512, 4);
    stage_inner_merged_iDFT4_FP((const complex_float*)cfftf512_twd3, x, y, 512, 16);
    stage_last_iDFT8xIN_8_FP(y, y, 512);
    return 0;
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU


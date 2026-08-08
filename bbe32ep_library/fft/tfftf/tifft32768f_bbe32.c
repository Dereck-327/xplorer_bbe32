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
#include "fft_tw.h"
#define IS_INV_FFTF 1
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"

/* Tables of a twiddles */
#include "fft_fp_tw.h"


/*-------------------------------------------------------------------------
Radix-2 floating point inverse FFT on complex data with reduced twiddle table

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Input spectrum samples
  Output:            
    y[N]        Complex output signal
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/



#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(int, tifft32768f, (complex_float * restrict y, complex_float * restrict x) )
#else
int tifft32768f(complex_float * restrict y, complex_float * restrict x)
{


    stage_first_tifft_IDFT4_FP((const complex_float*)tfftf32768_twd1, x, y, 32768);
    stage_inner_R3_iDFT4xIv_FP((const complex_float*)cfftf32768_twd2, y, x, 32768, 4);
    stage_inner_iDFT4_v16_FP((const complex_float*)cfftf32768_twd3, x, y, 32768, 16);
    stage_inner_R3_iDFT4xIv_FP((const complex_float*)cfftf32768_twd4, y, x, 32768, 64);
    stage_inner_R3_iDFT4xIv_FP((const complex_float*)cfftf32768_twd5, x, y, 32768, 256);
    stage_inner_R3_iDFT4xIv_FP((const complex_float*)cfftf32768_twd6, y, x, 32768, 1024);
    stage_last_iDFT8xIN_8_FP(x, y, 32768);
    return 0;
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU

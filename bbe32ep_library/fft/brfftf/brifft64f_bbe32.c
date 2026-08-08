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

/* Portable data types. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Common utility declarations. */
#include "common.h"
/* Twiddles tables for float point FFT */
#include "fft_fp_tw.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"


/*-------------------------------------------------------------------------
Blockwise inverse floating point FFT forming real data

These functions make inverse FFT forming real data on L blocks.
FFT implementation for real signal exploits the symmetry properties of the 
Fourier Transform: first, a blockwise complex FFT of half the original size
is applied to input data, and then the resulting spectrum undergoes a 
conversion procedure which results in complex spectrum of real input data.
NOTES:
1.  Bit-reversing permutation is done here.
2.  IFFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after
    the call.
3.  Spectrum of real signal is conjugate-symmetric about the midpoint, thus 
    real inverse FFT functions require only the first half of spectrum and 
    Nyquist frequency bin (altogether N/2+1 complex points) to reconstruct 
    the time domain signal. Moreover, the DC component and Nyquist 
    frequency component have zero imaginary part, thus they can be packed 
    in a single 'complex' value in input data. So, the m-th row of the input 
    matrix shall contain the following values:
    - DC component of the signal in x[m][0].re
    - Nyquist frequency component in x[m][0].im
    - first half of the complex spectrum in x[m][1]...x[m][(N/2)-1]
These data are used to reconstruct N values that are stored into the m-th 
row of real output matrix y[m][N]. The reconstruction is accomplished 
through a special conversion of input spectrum, such that subsequent 
invocation of blockwise complex inverse FFT of size N/2 actually produces 
the desired real signal.

Representation: floating point

Parameters:
  Input:            
    x[L][N/2] half of input spectrum, see note above
  Output:          
    y[L][N]   Output real samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(void, brifft64f, (float32_t * restrict y, complex_float * restrict x, int L) )
#else
void brifft64f(float32_t * restrict y, complex_float * restrict x, int L)
{
  
    const  int N = 64;
    complex_float *twiddle_table = (complex_float*)rfftf64_twd1;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((uintptr_t)x != (uintptr_t)y);
    NASSERT(L>0);

    blkrifft_spec_conv_packed_fp(x, twiddle_table, N, L);
    bcifft32f((complex_float*)y, x, L);
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU





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
Blockwise radix-2 floating point forward FFT on real data

Description: These functions make forward real FFT on L blocks, each of N=2^n
complex samples, where n=4..7. 
FFT implementation for real signal exploits the symmetry properties of the 
Fourier Transform: first, a blockwise complex FFT of half the original size 
is applied to input data, and then the resulting spectrum undergoes a 
conversion procedure which results in complex spectrum of real input data.
NOTES:
1. Bit-reversing permutation is done here. 
2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
   the call. 
3. As the complex spectrum of a real FFT is conjugate symmetric about the 
   midpoint, the brfftf functions generate only the first (N/2)+1 points 
   of the FFT, with the first point being the DC component, and the last 
   point ÿ the Nyquist frequency component. For real signal these two 
   spectral components have zero imaginary part, thus they can be packed 
   in a single 'complex' value. Finally, the the m-th row of the output 
   matrix will contain the following values:
    - DC component of the signal in y[m][0].re
    - Nyquist frequency component in y[m][0].im
    - first half of the complex spectrum in y[m][1]...y[m][(N/2)-1]
   These data are used to reconstruct N values that are stored into the 
   m-th row of real output matrix y[m][N]. The reconstruction is 
   accomplished through a special conversion of input spectrum, such 
   that subsequent invocation of blockwise complex inverse FFT of size 
   N/2 actually produces the desired real signal.

Representation: floating point


Parameters:
  Input:            
    x[L][N]   Real input signal
  Output:          
    y[L][N/2] Half of output spectrum, see note above
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(void, brfft64f, (complex_float * restrict y, float32_t * restrict x, int L) )
#else
void brfft64f(complex_float * restrict y, float32_t * restrict x, int L)
{
 
    complex_float *twiddle_table = (complex_float*)rfftf64_twd1;


    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((uintptr_t)x != (uintptr_t)y);
    NASSERT(L>0);

    bcfft32f(y, (complex_float*)x, L);
    blkrfft_spec_conv_packed_fp(y, twiddle_table, 64, L);

}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU



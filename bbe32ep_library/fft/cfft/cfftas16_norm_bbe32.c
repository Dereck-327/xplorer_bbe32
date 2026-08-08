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
#include "fft_tw.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#if !(HAVE_FFT && 1)
DISCARD_FUN(int, cfftas16_norm, (complex_fract16* restrict y, complex_fract16 * restrict x ))
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
inline_ ATTRIBUTE_ALWAYS_INLINE int cfft16as_norm_loc( int16_t * y,
                                                  int16_t * x                                               )
{
        xb_vecNx16 * Y;
  const xb_vecNx16 * X;
  const xb_vecNx16 * C;

  xb_vecNx16 p0, p1;

  xb_vecNx16 tw0, tw1;
  xb_vecNx16 t0, t1, t2;
  vselN      sel_2of8_0, sel_2of8_1;


  int shiftSum;
  int shift0, shift1;
  int nsa;

  NASSERT_ALIGN32( y );
  NASSERT_ALIGN32( x );

  NASSERT_ALIGN32( fft16_tw );

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file

  X = (const xb_vecNx16*)x;
  C = (const xb_vecNx16*)fft16_tw;
  Y = (      xb_vecNx16*)y;

  t0 = 0;

  BBE_MOVSAV( t0 );
  BBE_MOVSBV( t0 );

  BBE_LVNX16_IP(  t0, C, +1*(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw0, C, +1*(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw1, C, +1*(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP(  t1, C, +1*(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP(  t2, C, +1*(2*BBE_SIMD_WIDTH) );

  sel_2of8_0 = BBE_MOVVSELNX16( t1, 0 );
  sel_2of8_1 = BBE_MOVVSELNX16( t2, 0 );

  //----------------------------------------------------------------------------
  // Stage 1

  {
    xb_vecNx16 a0, a1, b0, b1;

    BBE_LVNX16_IP( a0, X, +(2*BBE_SIMD_WIDTH) );
    BBE_LVNX16_IP( a1, X, +(2*BBE_SIMD_WIDTH) );

    //
    // 1st scaling stage (input data scaling).
    //

   

 
    shiftSum = 3;

    //
    // Radix-4 w/ twiddles.
    //

    BBE_FFTWMODE( 1 );

    b0 = BBE_FFTADDSSRNX16( a0, a1 );
    b1 = BBE_FFTSUBSSRNX16( a0, a1 );

    // Is safe after NSA detection, because t0 holds either 1 or -1j.
    b1 = BBE_MULNX16CPACKQ( b1, t0 );

    a0 = BBE_SELNX16I( b1, b0, BBE_SELI_EXTRACT_LO_HALVES );
    a1 = BBE_SELNX16I( b1, b0, BBE_SELI_EXTRACT_HI_HALVES );

    BBE_FFTWMODE( 2 );

    // Radix-2! States A and B hold zeros. We use ADD4 because it allows shift
    // amounts exceeding 1.
    b0 = BBE_FFTADD4SABNX16( a0, a1, 0, 0 );
    b1 = BBE_FFTADD4SABNX16( a0, a1, 2, 0 );

    p0 = BBE_MULNX16CPACKQ( b0, tw0 );
    p1 = BBE_MULNX16CPACKQ( b1, tw1 );
  }

  //----------------------------------------------------------------------------
  // Stage 2 and digit reverse permutation

  {
    xb_vecNx16 a0, a1, b0, b1;

    //
    // 2nd scaling stage.
    //

    BBE_WRANGE( 4 );

    BBE_RANGENX16( p0 );
    BBE_RANGENX16( p1 );

    nsa = BBE_RRANGE();

    shift0 = shift1 = (nsa == 0);
      
    XT_MOVEQZ( shift0, 1, nsa     );
    XT_MOVEQZ( shift1, 1, nsa - 1 );

    shiftSum += shift0 + shift1;

    //
    // Radix-4 w/o twiddles.
    //

    a0 = BBE_SELNX16( p1, p0, sel_2of8_0 );
    a1 = BBE_SELNX16( p1, p0, sel_2of8_1 );

    BBE_FFTWMODE( shift0 );

    b0 = BBE_FFTADDSSRNX16( a0, a1 );
    b1 = BBE_FFTSUBSSRNX16( a0, a1 );

    b1 = BBE_MULNX16CPACKQ( b1, t0 );

    a0 = BBE_SELNX16I( b1, b0, BBE_SELI_EXTRACT_LO_HALVES );
    a1 = BBE_SELNX16I( b1, b0, BBE_SELI_EXTRACT_HI_HALVES );

    BBE_FFTWMODE( shift1 );

    b0 = BBE_FFTADDSSRNX16( a0, a1 );
    b1 = BBE_FFTSUBSSRNX16( a0, a1 );

    BBE_WRANGE( 4 );

    BBE_SVRNX16_IP( b0, Y, +(2*BBE_SIMD_WIDTH) );
    BBE_SVRNX16_IP( b1, Y, +(2*BBE_SIMD_WIDTH) );


  }

  return (shiftSum);

} // cfft16as_loc()

int cfftas16_norm ( complex_fract16* restrict y, complex_fract16 * restrict x           )
{
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  return cfft16as_norm_loc((int16_t*)y, (int16_t*)x);
} /* cfftas16_norm() */
#endif

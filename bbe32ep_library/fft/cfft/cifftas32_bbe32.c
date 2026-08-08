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
    Radix-2 inverse FFT on complex data, auto scaling
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
DISCARD_FUN(int, cifftas32, (complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp))
#else
/*-------------------------------------------------------------------------
Radix-2 inverse FFT on complex data, auto scaling

Description: These functions make inverse FFT on complex data of power of 2
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
    x[N]        Input spectrum samples
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Complex output signal
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data, with positive numbers corresponding to the right
                shift. _norm-suffixed functions return strictly positive
                values, while for regular variants the total shift amount is
                bi-directional.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
inline_ ATTRIBUTE_ALWAYS_INLINE int cifftas32_loc( int16_t * y,
                                                  const int16_t * x,
                                                  int * bexp )
{
  //
  // DFT32 decomposition:
  //
  // y = DFT32*x;
  //
  // DFT32 -> ( DFT2 x I16 )*( L8_2 x I4 )*
  //          ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )*
  //          T32_8*( DFT4 x I8 )
  //

        xb_vecNx16 * Y;
  const xb_vecNx16 * X;
  const xb_vecNx16 * TWD;

  xb_vecNx16 p0, p1, p2, p3;
  xb_vecNx16 q0, q1, q2, q3;

  xb_vecNx16 tw01, tw02, tw03;
  xb_vecNx16 tw11, tw12, tw13;

  vsaN  vsa0;

  int shift, shiftSum, nsa;

  NASSERT_ALIGN32( y );
  NASSERT_ALIGN32( x );

  NASSERT_ALIGN32( fft32_tw );

  NASSERT( BBE_SIMD_WIDTH == 16 );

  X   = (const xb_vecNx16*)x;
  Y   = (      xb_vecNx16*)y;
  TWD = (const xb_vecNx16*)fft32_tw;

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file

  // T32_8, 8..31
  BBE_LVNX16_IP( tw01, TWD, +(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw02, TWD, +(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw03, TWD, +(2*BBE_SIMD_WIDTH) );
  // T8_2 x I4, 8..31
  BBE_LVNX16_IP( tw11, TWD, +(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw12, TWD, +(2*BBE_SIMD_WIDTH) );
  BBE_LVNX16_IP( tw13, TWD, +(2*BBE_SIMD_WIDTH) );

  //----------------------------------------------------------------------------
  // Stage 1: T32_8*( DFT4 x I8 )

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    a0 = BBE_LVNX16_I( X, +0*(2*BBE_SIMD_WIDTH) );
    a1 = BBE_LVNX16_I( X, +1*(2*BBE_SIMD_WIDTH) );
    a2 = BBE_LVNX16_I( X, +2*(2*BBE_SIMD_WIDTH) );
    a3 = BBE_LVNX16_I( X, +3*(2*BBE_SIMD_WIDTH) );

    //
    // 1st scaling stage.
    //

    vsa0 = BBE_MOVVSA32( *bexp );

    a0 = BBE_SLLNX16( a0, vsa0 );
    a1 = BBE_SLLNX16( a1, vsa0 );
    a2 = BBE_SLLNX16( a2, vsa0 );
    a3 = BBE_SLLNX16( a3, vsa0 );

    shiftSum = 3 - (*bexp);

    BBE_FFTWMODE( 0x10 | 3 );

    //
    // DFT4 x I8
    //

    BBE_MOVSAV( a2 );
    BBE_MOVSBV( a3 );

    b0 = BBE_FFTADD4SABNX16( a0, a1, 0, 0 );
    b3 = BBE_FFTADD4SABNX16( a0, a1, 1, 0 );
    b2 = BBE_FFTADD4SABNX16( a0, a1, 2, 0 );
    b1 = BBE_FFTADD4SABNX16( a0, a1, 3, 0 );

    //
    // T32_8
    //

    p0 = b0;
    p1 = BBE_MULNX16JPACKQ( b1, tw01 );
    p2 = BBE_MULNX16JPACKQ( b2, tw02 );
    p3 = BBE_MULNX16JPACKQ( b3, tw03 );
  }

  //----------------------------------------------------------------------------
  // Stage 2: ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //
    // 2nd scaling stage
    //

    BBE_WRANGE( 4 );

    BBE_RANGENX16( p0 );
    BBE_RANGENX16( p1 );
    BBE_RANGENX16( p2 );
    BBE_RANGENX16( p3 );

    nsa = BBE_RRANGE();

    shift = 3 - nsa;

    XT_MOVLTZ( shift, 0, shift );

    shiftSum += shift;

    BBE_FFTWMODE( 0x10 | shift );

    // I2 x L16_8
    BBE_DSELNX16I( a1, a0, p1, p0, BBE_DSELI_INTERLEAVE_2 );
    BBE_DSELNX16I( a3, a2, p3, p2, BBE_DSELI_INTERLEAVE_2 );

    // L16_8 x I2
    BBE_DSELNX16I( b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4 );
    BBE_DSELNX16I( b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4 );

    //
    // DFT4 x I8
    //

    BBE_MOVSAV( b2 );
    BBE_MOVSBV( b3 );

    a0 = BBE_FFTADD4SABNX16( b0, b1, 0, 0 );
    a3 = BBE_FFTADD4SABNX16( b0, b1, 1, 0 );
    a2 = BBE_FFTADD4SABNX16( b0, b1, 2, 0 );
    a1 = BBE_FFTADD4SABNX16( b0, b1, 3, 0 );

    //
    // T8_2 x I4
    //

    q0 = a0;
    q1 = BBE_MULNX16JPACKQ( a1, tw11 );
    q2 = BBE_MULNX16JPACKQ( a2, tw12 );
    q3 = BBE_MULNX16JPACKQ( a3, tw13 );
  }

  //----------------------------------------------------------------------------
  // Stage 3: ( DFT2 x I16 )*( L8_2 x I4 )

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //
    // 3rd scaling stage
    //

    BBE_WRANGE( 4 );

    BBE_RANGENX16( q0 );
    BBE_RANGENX16( q1 );
    BBE_RANGENX16( q2 );
    BBE_RANGENX16( q3 );

    nsa = BBE_RRANGE();

    shift = ( nsa == 0 );

    shiftSum += shift;

    BBE_FFTWMODE( shift );

    // L8_2 x I4
    a0 = BBE_SELNX16I( q1, q0, BBE_SELI_EXTRACT_LO_HALVES );
    a1 = BBE_SELNX16I( q3, q2, BBE_SELI_EXTRACT_LO_HALVES );
    a2 = BBE_SELNX16I( q1, q0, BBE_SELI_EXTRACT_HI_HALVES );
    a3 = BBE_SELNX16I( q3, q2, BBE_SELI_EXTRACT_HI_HALVES );

    // DFT2 x I16
    b0 = BBE_FFTADDSSRNX16( a0, a2 );
    b1 = BBE_FFTADDSSRNX16( a1, a3 );
    b2 = BBE_FFTSUBSSRNX16( a0, a2 );
    b3 = BBE_FFTSUBSSRNX16( a1, a3 );

    BBE_WRANGE( 4 );

    BBE_SVRNX16_IP( b0, Y, +(2*BBE_SIMD_WIDTH) );
    BBE_SVRNX16_IP( b1, Y, +(2*BBE_SIMD_WIDTH) );
    BBE_SVRNX16_IP( b2, Y, +(2*BBE_SIMD_WIDTH) );
    BBE_SVRNX16_IP( b3, Y, +(2*BBE_SIMD_WIDTH) );

    (*bexp) = BBE_RRANGE();
  }

  return (shiftSum);

} // cifftas32_loc()

int cifftas32 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp )
{
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(y);

  return cifftas32_loc((int16_t*)y, (int16_t*)x, &bexp);
} /* cifftas32() */
#endif

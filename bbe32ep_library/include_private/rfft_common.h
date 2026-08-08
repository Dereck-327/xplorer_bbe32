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
  NatureDSP_Baseband library. FFT part.
    FFT/IFFT on real data
    Spectrum conversion routines for real-valued forward and inverse FFT
    kernels; optimized for BBE32EP.
  IntegrIT, 2006-2017
*/

#ifndef __RFFT_COMMON_H
#define __RFFT_COMMON_H

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#if HAVE_FFT
// In-place spectrum converter for a forward RFFT. Real-valued FFT of size N is
// performed in two steps: first, we apply a complex-valued FFT of size N/2,
// and then convert the resulting spectrum. Conversion algorithm exploits the
// symmetry properties of a real-valued Fourier transform and consists of a
// single radix-2 decimation-in-time (DIT) stage.
//
//  x[2*(N/2+1)] [in/out]  Input/output spectrum; N/2 complex samples on input,
//                         N/2+1 complex samples on output
//  twd[2*(N/4)] [in    ]  Twiddle factor table
//  bexp         [in    ]  Common block exponent over input data x[2*N/2], 0..15
//
// Returns the right shift amount that has been applied to converted
// spectrum, 0..1
//

inline_ ATTRIBUTE_ALWAYS_INLINE int rfft_spec_conv( 
                                             int16_t * x,
                                       const int16_t * twd,
                                       int N, int bexp );

// In-place spectrum converter for an inverse RFFT. Real-valued IFFT of size N
// is performed in two steps: first, we convert the spectrum, and then apply a
// complex-valued IFFT of size N/2. Resulting N/2 complex samples actually
// contain N/2 even-numbered samples (real parts) interleaved with N/2 odd-
// numbered samples (imaginary parts) of the desired real signal. Spectrum
// conversion algorithm exploits the symmetry properties of a real-valued
// Fourier transform and consists of a single radix-2 decimation-in-time (DIT)
// inversion stage.
//
//  x[2*(N/2+1)]   [in/out]  Spectrum, N/2+1 complex samples on input, N/2
//                           samples on output
//  twd[2*N/4]     [in    ]  Twiddle factor table, CQ15
//  bexp           [in/out]  Common block exponent over input/output data, 0..15
//

inline_ ATTRIBUTE_ALWAYS_INLINE int rifft_spec_conv( 
                                             int16_t * x,
                                       const int16_t * twd,
                                       int N, int * bexp );

// In-place real-to-complex spectrum converter for the forward RFFT.
int rfft_spec_conv( int16_t * x,
              const int16_t * twd,
              int N, int bexp )
{
  //
  // MATLAB code:
  //  a0 = X(1:N/4);
  //  a1 = [X(1);X(N/2-(0:N/4-2))];
  //  twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
  //  y = [a0.*(1-twd)/2+conj(a1).*(1+twd)/2; ...
  //       conj(X(N/4+1)); ...
  //       wrev(conj(a0).*(1+conj(twd))/2+a1.*(1-conj(twd))/2)];
  //
  // The optimized implementation starts from the central sample (that is
  // 1/2 of Nyquist frequency) and runs to the right and left edges of spectrum.
  //

  const xb_vecNx16 * TW;
  const xb_vecNx16 * X0;
  const xb_vecNx16 * X1;
        xb_vecNx16 * Y0;
        xb_vecNx16 * Y1;

  valign X1_va, Y1_va;

  xb_vecNx16 a0, a1;
  xb_vecNx40 b0, b1;
  xb_vecNx16 c0, c1;

  xb_vecNx16 tw0, tw1, tw2;

  vsaN vsa_cnv;

  int shift;

  int n;

  NASSERT_ALIGN32( x   );
  NASSERT_ALIGN32( twd );

  ASSERT( !( N % (2*BBE_SIMD_WIDTH) ) );

  // Twiddle table is read starting from the end.
  TW = (const xb_vecNx16*)( (uintptr_t)twd + 4*(N/2) - 2*BBE_SIMD_WIDTH );

  // Setup left/right in/out data pointers to the middle of spectrum.
  X0 = (const xb_vecNx16*)( (uintptr_t)x + 4*(N/4) - 2*BBE_SIMD_WIDTH );
  X1 = (const xb_vecNx16*)( (uintptr_t)x + 4*(N/4)                    );
  Y0 = (      xb_vecNx16*)( (uintptr_t)x + 4*(N/4) - 2*BBE_SIMD_WIDTH );

  WUR_CBEGIN( (uintptr_t)x           );
  WUR_CEND  ( (uintptr_t)x + 4*(N/2) );

  //
  // Setup the scaling shift amount.
  //

  // Data scaling shift motivation:
  //  1. Data should be normalized using the common block exponent.
  //  2. Forward spectrum conversion doesn't change data magnitude, but includes
  //     twiddle multiplication, which requires a reservation of 1 bit position.
  shift = 1 - bexp;

  // Setup for Q30->Q15 conversion coupled with data scaling
  vsa_cnv = BBE_MOVVSA32( 15 + shift );

  //
  // Transform the central sample:
  //  y(N/4+1) = conj(X(N/4+1));
  //

  // X(N/4+1)
  BBE_LPNX16_IP( a1, X1, +4 );

  // Initialize Y1 by updated X1 to prevent the compiler from scheduling
  // BBE_LANX16POS_PC and BBE_SPNX16_IP on the same cycle.
  Y1 = (xb_vecNx16*)X1;

  BBE_LANX16POS_PC( X1_va, X1 );

  // Q30 <- Q15 + 15
  b0 = BBE_UNPKQNX16( a1 );
  b0 = BBE_CONJNX40C( b0 );
  b0 = BBE_RNDADJNX40( b0, vsa_cnv );
  // Q15 <- Q30/2^shift - 15
  c1 = BBE_PACKVNX40 ( b0, vsa_cnv );

  // Y(N/4+1)
  BBE_SPNX16_I( c1, Y1, -4 );

  Y1_va = BBE_ZALIGN();

  //
  // Run to the edges of spectrum, converting 8 left-hand and 8 right-hand
  // samples at each trip.
  //

  __Pragma( "ymemory( X0 )" );
  __Pragma( "ymemory( TW )" );
  for ( n=0; n<(N/4)/(BBE_SIMD_WIDTH/2); n++ )
  {
    // a0 = X(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1));
    BBE_LVNX16_IP( a0, X0, -2*BBE_SIMD_WIDTH );
    // a1 = X(N/4+1+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2)); 
    BBE_LANX16_IC( a1, X1_va, X1 );

    // a1 = X(N/4+1+n*SIMD_WIDTH_2+(SIMD_WIDTH_2:-1:1)); 
    a1 = BBE_SHFLNX16I( a1, BBE_SHFLI_REVERSE_2 );

    // tw = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // tw1 = (1+tw)/2; Q15
    BBE_LVNX16_IP( tw1, TW, -2*BBE_SIMD_WIDTH );
    // tw0 = (1-tw)/2; Q15
    BBE_LVNX16_IP( tw0, TW, -2*BBE_SIMD_WIDTH );
    // tw2 = (1+conj(tw))/2;
    tw2 = BBE_CONJSNX16C( tw1 );

    // b0 = a0*(1-tw)/2; Q30 <- Q15*Q15
    b0 = BBE_MULRNX16C( tw0, a0, vsa_cnv );
    // b1 = conj(a0)*(1+conj(tw))/2; Q30 <- Q15*Q15
    b1 = BBE_MULRNX16J( tw2, a0, vsa_cnv );

    // b0 += conj(a1)*(1+tw)/2; Q30 <- Q15*Q15
    BBE_MULANX16J( b0, tw1, a1 );
    // b1 += a1*(1-conj(tw))/2; Q30 <- Q15*Q15
    BBE_MULANX16J( b1, a1, tw0 );

    // Q15 <- Q30/2^shift - 15 w/ rounding and saturation
    c0 = BBE_PACKVNX40( b0, vsa_cnv );
    c1 = BBE_PACKVNX40( b1, vsa_cnv );

    // Right-hand vector indices: N/4+1+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2)
    c1 = BBE_SHFLNX16I( c1, BBE_SHFLI_REVERSE_2 );

    // Store y(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1))
    BBE_SVNX16_IP( c0, Y0, -2*BBE_SIMD_WIDTH );
    // Store y(N/4+1+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2))
    BBE_SANX16_IP( c1, Y1_va, Y1 );
  }

  BBE_SANX16POS_FP( Y1_va, Y1 );

  return (shift);

} // rfft_spec_conv()

int rifft_spec_conv( int16_t * x,
               const int16_t * twd,
               int N, int * bexp )
{
  //
  // MATLAB code:
  //  a0 = x(1:N/4);
  //  a1 = x(N/2+1-(0:N/4-1));
  //  twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
  //  y = [a0.*(1-conj(twd))+conj(a1).*(1+conj(twd)); ...
  //       2*conj(x(N/4+1)); ...
  //       wrev(conj(a0).*(1+twd)+a1.*(1-twd))];
  //
  // The optimized implementation starts from the central sample (that is
  // 1/2 of Nyquist frequency) and runs to the right and left edges of spectrum.
  //

  const xb_vecNx16 * TW;
  const xb_vecNx16 * X0;
  const xb_vecNx16 * X1;
        xb_vecNx16 * Y0;
        xb_vecNx16 * Y1;

  xb_vecNx16 a0, a1, a2;
  xb_vecNx40 b0, b1;
  xb_vecNx16 c0, c1, c2;

  xb_vecNx16 tw0, tw1, tw2;

  vselN sel_sft;
  vsaN  vsa_cnv;

  int shift;

  int n;

  NASSERT_ALIGN32( x   );
  NASSERT_ALIGN32( twd );

  ASSERT( !( N % (2*BBE_SIMD_WIDTH) ) );

  ASSERT( bexp && 0<=(*bexp) && (*bexp)<=15 );

  // Twiddles table is read starting from the end.
  TW = (const xb_vecNx16*)( (uintptr_t)twd + 2*N );

  // Setup left/right in/out data pointers to the middle of spectrum.
  X0 = (const xb_vecNx16*)( (uintptr_t)x + N - 2*BBE_SIMD_WIDTH );
  X1 = (const xb_vecNx16*)( (uintptr_t)x + N                    );
  Y0 = (      xb_vecNx16*)( (uintptr_t)x + N - 2*BBE_SIMD_WIDTH );
  Y1 = (      xb_vecNx16*)( (uintptr_t)x + N                    );

  // Data scaling shift motivation:
  //  1. Data should be normalized using input common block exponent.
  //  2. Inverse spectrum conversion stage gains input data by a factor of 2,
  //     and it also includes twiddle multiplication. Altogether 2 bit 
  //     positions are to be reserved.
  shift = 2 - (*bexp);

  // Setup for Q29->Q15 conversion coupled with data scaling
  vsa_cnv = BBE_MOVVSA32( 14 + shift );

  // 
  // Setup the shift-and-reverse pattern for select
  //

  BBE_LVNX16_IP( a0, TW, -2*BBE_SIMD_WIDTH );

  sel_sft = BBE_MOVVSELNX16( a0, 0 );

  //
  // Transform the central sample:
  //  y(N/4+1) = 2*conj(x(N/4+1));
  //

  // a1 = X(N/4+(1:SIMD_WIDTH_2));
  BBE_LVNX16_IP( a1, X1, +2*BBE_SIMD_WIDTH );

  // Q29 <- 2*Q15 + 14
  b0 = BBE_UNPKQNX16 ( a1          );
  b0 = BBE_CONJNX40C ( b0          );
  b0 = BBE_RNDADJNX40( b0, vsa_cnv );
  // Q15 <- Q29/2^shift - 14
  c1 = BBE_PACKVNX40( b0, vsa_cnv );

  //
  // Run to the edges of spectrum, converting 8 left-hand and 8 right-hand
  // samples at each trip.
  //

  BBE_WRANGE( 4 );

  __Pragma( "ymemory( X0 )" );
  __Pragma( "ymemory( TW )" );
  for ( n=0; n<(N/4)/(BBE_SIMD_WIDTH/2); n++ )
  {
    // a0 = X(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1));
    BBE_LVNX16_IP( a0, X0, -2*BBE_SIMD_WIDTH );
    // a1 = X(N/4+(n+1)*SIMD_WIDTH_2+(1:SIMD_WIDTH_2)); 
    BBE_LVNX16_IP( a2, X1, +2*BBE_SIMD_WIDTH );
    // a1 = X(N/4+1+n*SIMD_WIDTH_2+(SIMD_WIDTH_2:-1:1)); 
    a1 = BBE_SELNX16( a2, a1, sel_sft );

    // tw = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // tw1 = (1+tw)/2; Q15
    BBE_LVNX16_IP( tw1, TW, -2*BBE_SIMD_WIDTH );
    // tw0 = (1-tw)/2; Q15
    BBE_LVNX16_IP( tw0, TW, -2*BBE_SIMD_WIDTH );
    // tw2 = (1+conj(tw))/2;
    tw2 = BBE_CONJSNX16C( tw1 );

    // b0 = a0*conj(1-twd); Q29 <- Q15*( 2*Q15 ) - 1
    b0 = BBE_MULRNX16J( a0, tw0, vsa_cnv );
    // b1 = conj(a0)*(1+twd); Q29 <- Q15*( 2*Q15 ) - 1
    b1 = BBE_MULRNX16J( tw1, a0, vsa_cnv );

    // b0 = b0 + conj(a1)*(1+conj(twd));
    BBE_MULANX16J( b0, tw2, a1 );
    // b1 = b1 + a1.*(1-twd);
    BBE_MULANX16C( b1, a1, tw0 );

    // Q15 <- Q29/2^shift - 14 w/ saturation and rounding
    c0 = BBE_PACKVNX40( b0, vsa_cnv );
    c2 = BBE_PACKVNX40( b1, vsa_cnv );

    // Right-hand vector indices: N/4+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2)
    c1 = BBE_SELNX16( c1, c2, sel_sft );

    // Store y(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1))
    BBE_SVRNX16_IP( c0, Y0, -2*BBE_SIMD_WIDTH );
    // Store y(N/4+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2))
    BBE_SVRNX16_IP( c1, Y1, +2*BBE_SIMD_WIDTH );

    a1 = a2;
    c1 = c2;
  }

  *bexp = BBE_RRANGE();

  return (shift);

} // rifft_spec_conv()
#endif // HAVE_FFT
#endif // __RFFT_COMMON_H

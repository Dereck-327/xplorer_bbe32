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
    Radix-2 inverse FFT forming real data, auto scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Twiddle tables declarations. */
#include "fft_tw.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(int, rifft64, (int16_t * restrict y, complex_fract16 * restrict x, int bexp))
#else
/*-------------------------------------------------------------------------
Radix-2 inverse FFT forming real data, auto scaling

Description: These functions make inverse FFT forming real data of length
N=2^n, n=4..15. Algorithm exploits the symmetry properties of the FFT:
the input spectrum is modified in such a way that a complex-valued inverse FFT
of half the original size that is applied to the transformed spectrum actually
results in real data.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:
    x[(N/2+1)]    Input spectrum samples
    bexp          Common block exponent, that is the minimum number of redundant
                  sign bits over input data x[]
  Output:
    y[N]          Real output signal
  Returned value:
                  Total shift amount applied throughout the transform to scale
                  the data. Total shift is bi-directional, with positive numbers
                  corresponding to the right shift.
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

int rifft64 ( int16_t * restrict y, complex_fract16 * restrict x, int bexp )
{
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          C;

  xb_vecNx16 tw01, tw02, tw03;
  xb_vecNx16 tw11, tw12, tw13;

  xb_vecNx16 twc00, twc01, twc02;
  xb_vecNx16 twc10, twc11, twc12;

  xb_vecNx16 p0, p1, p2, p3;
  xb_vecNx16 q0, q1, q2, q3;
  xb_vecNx16 r0, r1, r2, r3;

  vselN sel_rsft;
  vsaN  vsa0;

  int shiftSum, shift;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file.

  C = (const xb_vecNx16*)fft_tw_tab_rfft_64;

  // T32_8, 8..31
  tw01 = BBE_LVNX16_I(C, 0 * 2 * BBE_SIMD_WIDTH);
  tw02 = BBE_LVNX16_I(C, 1 * 2 * BBE_SIMD_WIDTH);
  tw03 = BBE_LVNX16_I(C, 2 * 2 * BBE_SIMD_WIDTH);
  // T8_2 x I4, 8..31
  tw11 = BBE_LVNX16_I(C, 3 * 2 * BBE_SIMD_WIDTH);
  tw12 = BBE_LVNX16_I(C, 4 * 2 * BBE_SIMD_WIDTH);
  tw13 = BBE_LVNX16_I(C, 5 * 2 * BBE_SIMD_WIDTH);
  // N=64; reshape(1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2)),8,N/32).'
  twc10 = BBE_LVNX16_I(C, 6 * 2 * BBE_SIMD_WIDTH);
  twc11 = BBE_LVNX16_I(C, 7 * 2 * BBE_SIMD_WIDTH);
  twc00 = BBE_LVNX16_I(C, 8 * 2 * BBE_SIMD_WIDTH);
  twc01 = BBE_LVNX16_I(C, 9 * 2 * BBE_SIMD_WIDTH);
  // Shift-and-reverse pattern for select.
  p0 = BBE_LVNX16_I(C, 10 * 2 * BBE_SIMD_WIDTH);

  sel_rsft = BBE_MOVVSV(p0, 0);

  //----------------------------------------------------------------------------
  // Real-to-complex spectrum converter.

  X = (const xb_vecNx16*)x;

  {
    xb_vecNx16 x0, x1, x2, x3, x4;

    xb_vecNx16 a00, a01, a10, a20, a21;
    xb_vecNx40 b00, b01, b10, b20, b21;
    xb_vecNx16 c00, c01, c10, c20, c21;

    //
    // 1st scaling stage.
    //

    // Data scaling shift motivation:
    //  1. Data should be normalized using input common block exponent.
    //  2. Inverse spectrum conversion stage gains input data by a factor of 2,
    //     and it also includes twiddle multiplication. Altogether 2 bit 
    //     positions are to be reserved.
    shift = 2 - bexp;

    // Setup for Q29->Q15 conversion coupled with data scaling
    vsa0 = BBE_MOVVSA32(14 + shift);

    shiftSum = shift;

    //
    // Spectrum conversion. MATLAB code:
    //  N = 64;
    //  a0 = x(1:N/4);
    //  a1 = x(N/4+1);
    //  a2 = x(N/2+1-(0:N/4-1));
    //  twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
    //  y = [a0.*(1-conj(twd))+conj(a2).*(1+conj(twd)); ...
    //       2*conj(a1); ...
    //       wrev(conj(a0).*(1+twd)+a2.*(1-twd))];
    //

    x0 = BBE_LVNX16_I(X, 0 * 4 * BBE_SIMD_WIDTH / 2);
    x1 = BBE_LVNX16_I(X, 1 * 4 * BBE_SIMD_WIDTH / 2);
    x2 = BBE_LVNX16_I(X, 2 * 4 * BBE_SIMD_WIDTH / 2);
    x3 = BBE_LVNX16_I(X, 3 * 4 * BBE_SIMD_WIDTH / 2);
    x4 = BBE_LPNX16_I(X, 4 * 4 * BBE_SIMD_WIDTH / 2);

    // a0 = x(1:N/4);
    a01 = x0;
    a00 = x1;
    // a1 = x(N/4+1);
    a10 = x2;
    // a2 = x(N/2+1-(0:N/4-1));
    a20 = BBE_SELNX16(x3, x2, sel_rsft);
    a21 = BBE_SELNX16(x4, x3, sel_rsft);

    // twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // twc0 = (1-twd)/2; CQ15
    // twc1 = (1+twd)/2; CQ15
    // twc2 = (1+conj(twd))/2; CQ15
    twc02 = BBE_CONJSNX16C(twc01);
    twc12 = BBE_CONJSNX16C(twc11);

    // b0 = 2*a0.*(1-conj(twd))/2; CQ29 <- 2*CQ15*CQ15 - 1
    b01 = BBE_MULRNX16J(a01, twc10, vsa0);
    b00 = BBE_MULRNX16J(a00, twc00, vsa0);
    // CQ29 <- 2*CQ15 + 14
    b10 = BBE_UNPKQNX16(a10);
    // b2 = 2*conj(a0).*(1+twd)/2; CQ29 <- 2*CQ15*CQ15 - 1
    b20 = BBE_MULRNX16J(twc01, a00, vsa0);
    b21 = BBE_MULRNX16J(twc11, a01, vsa0);

    // b0 += 2*conj(a2).*(1+conj(twd))/2;
    BBE_MULANX16J(b01, twc12, a21);
    BBE_MULANX16J(b00, twc02, a20);
    // b1 = 2*conj(a1); 
    b10 = BBE_CONJNX40C(b10);
    // b2 += 2*a2.*(1-twd)/2;
    BBE_MULANX16C(b20, twc00, a20);
    BBE_MULANX16C(b21, twc10, a21);

    b10 = BBE_RNDADJNX40(b10, vsa0);

    // CQ15 <- CQ29/2^shift0 - 14
    c01 = BBE_PACKVNX40(b01, vsa0);
    c00 = BBE_PACKVNX40(b00, vsa0);
    c10 = BBE_PACKVNX40(b10, vsa0);
    c20 = BBE_PACKVNX40(b20, vsa0);
    c21 = BBE_PACKVNX40(b21, vsa0);

    p0 = c01;
    p1 = c00;
    p2 = BBE_SELNX16(c10, c20, sel_rsft);
    p3 = BBE_SELNX16(c20, c21, sel_rsft);
  }

  //----------------------------------------------------------------------------
  // Complex-valued 16-point inverse FFT:
  //  IDFT32 = ( DFT2 x I16 )*( L8_2 x I4 )*
  //           ( T8_2' x I4 )*( IDFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )*
  //           T32_8'*( IDFT4 x I8 )

  //
  // Stage 1: T32_8'*( IDFT4 x I8 )
  //

  {
    xb_vecNx16 a0, a1, a2, a3;

    //
    // 2nd scaling stage.
    //

    // Fixed, radix-4 with twiddle multiplication.
    shiftSum += (shift = 3);

    BBE_FFTWMODE(0x00 | shift);

    //
    // IDFT4 x I8
    //

    BBE_MOVSAV(p2);
    BBE_MOVSBV(p3);

    a0 = BBE_FFTADD4SABNX16(p0, p1, 0, 0);
    a1 = BBE_FFTADD4SABNX16(p0, p1, 1, 0);
    a2 = BBE_FFTADD4SABNX16(p0, p1, 2, 0);
    a3 = BBE_FFTADD4SABNX16(p0, p1, 3, 0);

    //
    // T32_8'
    //

    q0 = a0;
    q1 = BBE_MULNX16JPACKQ(a1, tw01);
    q2 = BBE_MULNX16JPACKQ(a2, tw02);
    q3 = BBE_MULNX16JPACKQ(a3, tw03);
  }

  //
  // Stage 2: ( T8_2' x I4 )*( IDFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )
  //

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //
    // 3nd scaling stage
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(q0);
    BBE_RANGENX16(q1);
    BBE_RANGENX16(q2);
    BBE_RANGENX16(q3);

    bexp = BBE_RRANGE();

    shift = 3 - bexp;

    XT_MOVLTZ(shift, 0, shift);

    shiftSum += shift;

    BBE_FFTWMODE(0x00 | shift);

    // I2 x L16_8
    BBE_DSELNX16I(a1, a0, q1, q0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(a3, a2, q3, q2, BBE_DSELI_INTERLEAVE_2);

    // L16_8 x I2
    BBE_DSELNX16I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4);

    //
    // IDFT4 x I8
    //

    BBE_MOVSAV(b2);
    BBE_MOVSBV(b3);

    a0 = BBE_FFTADD4SABNX16(b0, b1, 0, 0);
    a1 = BBE_FFTADD4SABNX16(b0, b1, 1, 0);
    a2 = BBE_FFTADD4SABNX16(b0, b1, 2, 0);
    a3 = BBE_FFTADD4SABNX16(b0, b1, 3, 0);

    //
    // T8_2' x I4
    //

    r0 = a0;
    r1 = BBE_MULNX16JPACKQ(a1, tw11);
    r2 = BBE_MULNX16JPACKQ(a2, tw12);
    r3 = BBE_MULNX16JPACKQ(a3, tw13);
  }

  //----------------------------------------------------------------------------
  // Stage 3: ( DFT2 x I16 )*( L8_2 x I4 )

  Y = (xb_vecNx16*)y;

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //
    // 4th scaling stage
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(r0);
    BBE_RANGENX16(r1);
    BBE_RANGENX16(r2);
    BBE_RANGENX16(r3);

    bexp = BBE_RRANGE();

    shift = (bexp == 0);

    shiftSum += shift;

    BBE_FFTWMODE(shift);

    // L8_2 x I4
    a0 = BBE_SELNX16I(r1, r0, BBE_SELI_EXTRACT_LO_HALVES);
    a1 = BBE_SELNX16I(r3, r2, BBE_SELI_EXTRACT_LO_HALVES);
    a2 = BBE_SELNX16I(r1, r0, BBE_SELI_EXTRACT_HI_HALVES);
    a3 = BBE_SELNX16I(r3, r2, BBE_SELI_EXTRACT_HI_HALVES);

    // DFT2 x I16
    b0 = BBE_FFTADDSSRNX16(a0, a2);
    b1 = BBE_FFTADDSSRNX16(a1, a3);
    b2 = BBE_FFTSUBSSRNX16(a0, a2);
    b3 = BBE_FFTSUBSSRNX16(a1, a3);

    BBE_SVNX16_I(b0, Y, 0 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(b1, Y, 1 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(b2, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(b3, Y, 3 * 2 * BBE_SIMD_WIDTH);
  }
  return (shiftSum);
} /* rifft64() */
#endif

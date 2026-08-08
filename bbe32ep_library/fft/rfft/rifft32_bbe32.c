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
/* Spectrum conversion routines for real-valued FFTs. */
#include "rfft_common.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(int, rifft32, (int16_t * restrict y, complex_fract16 * restrict x, int bexp))
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

int rifft32 ( int16_t * restrict y, complex_fract16 * restrict x, int bexp )
{
  xb_vecNx16 * Y;
  const xb_vecNx16 * X;
  const xb_vecNx16 * C;

  xb_vecNx16 tw0, tw1, tw2, tw3, tw4, tw5, tw6, tw7;
  xb_vecNx16 tw_T16_8;
  xb_vecNx16 twc0, twc1, twc2;

  vselN sel_rsft;
  vsaN  vsa0;

  xb_vecNx16 p0, p1;

  int shift0, shift1;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file.

  C = (const xb_vecNx16*)fft_tw_tab_rfft_32;

  tw0 = BBE_LVNX16_I(C, 0 * 2 * BBE_SIMD_WIDTH);
  tw1 = BBE_LVNX16_I(C, 1 * 2 * BBE_SIMD_WIDTH);
  tw2 = BBE_LVNX16_I(C, 2 * 2 * BBE_SIMD_WIDTH);
  tw3 = BBE_LVNX16_I(C, 3 * 2 * BBE_SIMD_WIDTH);
  tw4 = BBE_LVNX16_I(C, 4 * 2 * BBE_SIMD_WIDTH);
  tw5 = BBE_LVNX16_I(C, 5 * 2 * BBE_SIMD_WIDTH);
  tw6 = BBE_LVNX16_I(C, 6 * 2 * BBE_SIMD_WIDTH);
  tw7 = BBE_LVNX16_I(C, 7 * 2 * BBE_SIMD_WIDTH);
  tw_T16_8 = BBE_LVNX16_I(C, 8 * 2 * BBE_SIMD_WIDTH);
  twc0 = BBE_LVNX16_I(C, 9 * 2 * BBE_SIMD_WIDTH);
  twc1 = BBE_LVNX16_I(C, 10 * 2 * BBE_SIMD_WIDTH);
  p0 = BBE_LVNX16_I(C, 11 * 2 * BBE_SIMD_WIDTH);

  sel_rsft = BBE_MOVVSV(p0, 0);

  //----------------------------------------------------------------------------
  // Real-to-complex spectrum converter.
  //

  X = (const xb_vecNx16*)x;

  {
    xb_vecNx16 x0, x1, x2;
    xb_vecNx16 a0, a1, a2;
    xb_vecNx40 b0, b1, b2;
    xb_vecNx16 c0, c1, c2;

    //
    // 1st scaling stage.
    //

    // Data scaling shift motivation:
    //  1. Data should be normalized using input common block exponent.
    //  2. Inverse spectrum conversion stage gains input data by a factor of 2,
    //     and it also includes twiddle multiplication. Altogether 2 bit 
    //     positions are to be reserved.
    shift0 = 2 - bexp;

    // Setup for Q29->Q15 conversion coupled with data scaling
    vsa0 = BBE_MOVVSA32(14 + shift0);

    //
    // Spectrum conversion. MATLAB code:
    //  N = 32;
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
    x2 = BBE_LPNX16_I(X, 2 * 4 * BBE_SIMD_WIDTH / 2);

    // a0 = x(1:N/4);
    a0 = x0;
    // a1 = x(N/4+1);
    a1 = x1;
    // a2 = x(N/2+1-(0:N/4-1));
    a2 = BBE_SELNX16(x2, x1, sel_rsft);

    // twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // twc0 = (1-twd)/2; CQ15
    // twc1 = (1+twd)/2; CQ15
    // twc2 = (1+conj(twd))/2; CQ15
    twc2 = BBE_CONJSNX16C(twc1);

    // b0 = 2*a0.*(1-conj(twd))/2; CQ29 <- 2*CQ15*CQ15 - 1
    b0 = BBE_MULRNX16J(a0, twc0, vsa0);
    // CQ29 <- 2*CQ15 + 14
    b1 = BBE_UNPKQNX16(a1);
    // b2 = 2*conj(a0).*(1+twd)/2; CQ29 <- 2*CQ15*CQ15 - 1
    b2 = BBE_MULRNX16J(twc1, a0, vsa0);

    // b0 += 2*conj(a2).*(1+conj(twd))/2;
    BBE_MULANX16J(b0, twc2, a2);
    // b1 = 2*conj(a1); 
    b1 = BBE_CONJNX40C(b1);
    // b2 += 2*a2.*(1-twd)/2;
    BBE_MULANX16C(b2, twc0, a2);

    b1 = BBE_RNDADJNX40(b1, vsa0);

    // CQ15 <- CQ29/2^shift0 - 14
    c0 = BBE_PACKVNX40(b0, vsa0);
    c1 = BBE_PACKVNX40(b1, vsa0);
    c2 = BBE_PACKVNX40(b2, vsa0);

    p0 = c0;
    p1 = BBE_SELNX16(c1, c2, sel_rsft);
  }

  //----------------------------------------------------------------------------
  // Complex-valued 16-point inverse FFT:
  //  IDFT16 = ( DFT2 x I8 )*T16_8*( I2 x IDFT8 )*L16_2
  //

  Y = (xb_vecNx16*)y;

  {
    xb_vecNx16 a0, a1;
    xb_vecNx16 b0, b1;
    xb_vecNx40 w0, w1;

    xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
    xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;

    //
    // 2nd scaling stage.
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(p0);
    BBE_RANGENX16(p1);

    bexp = BBE_RRANGE();

    // DFT16 w/o twiddles + data normalization
    shift1 = 4 - bexp;

    vsa0 = BBE_MOVVSA32(15 + shift1);

    //
    // T16_8*( I2 x DFT8 )*L16_2
    //

    a00 = BBE_REPNX16C(p0, 0); a10 = BBE_REPNX16C(p0, 1);
    a01 = BBE_REPNX16C(p0, 2); a11 = BBE_REPNX16C(p0, 3);
    a02 = BBE_REPNX16C(p0, 4); a12 = BBE_REPNX16C(p0, 5);
    a03 = BBE_REPNX16C(p0, 6); a13 = BBE_REPNX16C(p0, 7);
    a04 = BBE_REPNX16C(p1, 0); a14 = BBE_REPNX16C(p1, 1);
    a05 = BBE_REPNX16C(p1, 2); a15 = BBE_REPNX16C(p1, 3);
    a06 = BBE_REPNX16C(p1, 4); a16 = BBE_REPNX16C(p1, 5);
    a07 = BBE_REPNX16C(p1, 6); a17 = BBE_REPNX16C(p1, 7);

    w0 = BBE_MULNX16J(a00, tw0);
    w1 = BBE_MULNX16J(a10, tw0);

    BBE_MULANX16J(w0, a01, tw1); BBE_MULANX16J(w1, a11, tw1);
    BBE_MULANX16J(w0, a02, tw2); BBE_MULANX16J(w1, a12, tw2);
    BBE_MULANX16J(w0, a03, tw3); BBE_MULANX16J(w1, a13, tw3);
    BBE_MULANX16J(w0, a04, tw4); BBE_MULANX16J(w1, a14, tw4);
    BBE_MULANX16J(w0, a05, tw5); BBE_MULANX16J(w1, a15, tw5);
    BBE_MULANX16J(w0, a06, tw6); BBE_MULANX16J(w1, a16, tw6);
    BBE_MULANX16J(w0, a07, tw7); BBE_MULANX16J(w1, a17, tw7);

    w0 = BBE_RNDADJNX40(w0, vsa0);
    w1 = BBE_RNDADJNX40(w1, vsa0);

    a0 = BBE_PACKVNX40(w0, vsa0);
    a1 = BBE_PACKVNX40(w1, vsa0);

    a1 = BBE_MULNX16JPACKQ(a1, tw_T16_8);

    b0 = BBE_ADDNX16(a0, a1);
    b1 = BBE_SUBNX16(a0, a1);

    BBE_SVNX16_IP(b0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(b1, Y, 2 * BBE_SIMD_WIDTH);
  }

  return (shift0 + shift1);
} /* rifft32() */
#endif

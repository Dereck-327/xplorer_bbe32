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
    Radix-2 forward FFT on real data, auto scaling
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
DISCARD_FUN(int, rfft32, (complex_fract16 * restrict y, int16_t * restrict x, int bexp))
#else
/*-------------------------------------------------------------------------
Radix-2 forward FFT on real data, auto scaling

Description: These functions make FFT on real data of length N=2^n, n=4..15.
The algorithm exploits the symmetry properties of the FFT: first, a complex
FFT of half the original size is applied to input data, then the resulting
spectrum undergoes a postprocessing procedure which results in complex spectrum
of real input data.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:      
    x[N]          Real input signal
    bexp          Common block exponent, that is the minimum number of redundant
                  sign bits over input data x[N]
  Output:      
    y[(N/2+1)]    Output spectrum samples. 
  Returned value:
                  Total shift amount applied throughout the transform to scale
                  the data. Total shift is bi-directional, with positive numbers
                  corresponding to the right shift.
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

int rfft32 ( complex_fract16 * restrict y, int16_t * restrict x, int bexp )
{
  xb_vecNx16 * Y;
  const xb_vecNx16 * X;
  const xb_vecNx16 * C;

  valign Y_va;

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
  // Complex-valued 16-point FFT:
  //  DFT16 = ( DFT2 x I8 )*T16_8*( I2 x DFT8 )*L16_2
  //
  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)y;

 #if 1

  {
    xb_vecNx16 x0, x1;
    xb_vecNx16 a0, a1;
    xb_vecNx40 w0, w1;

    xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
    xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;

    x0 = BBE_LVNX16_I(X, 0 * 4 * BBE_SIMD_WIDTH / 2);
    x1 = BBE_LVNX16_I(X, 1 * 4 * BBE_SIMD_WIDTH / 2);

    //
    // 1st scaling stage.
    //

    // DFT16 w/o twiddles + input data normalization
    shift0 = 4 - bexp + (bexp==0);
   // shift0 = 4 - bexp ;
    vsa0 = BBE_MOVVSA32(15 + shift0);

    //
    // T16_8*( I2 x DFT8 )*L16_2
    //

    a00 = BBE_REPNX16C(x0, 0); a10 = BBE_REPNX16C(x0, 1);
    a01 = BBE_REPNX16C(x0, 2); a11 = BBE_REPNX16C(x0, 3);
    a02 = BBE_REPNX16C(x0, 4); a12 = BBE_REPNX16C(x0, 5);
    a03 = BBE_REPNX16C(x0, 6); a13 = BBE_REPNX16C(x0, 7);
    a04 = BBE_REPNX16C(x1, 0); a14 = BBE_REPNX16C(x1, 1);
    a05 = BBE_REPNX16C(x1, 2); a15 = BBE_REPNX16C(x1, 3);
    a06 = BBE_REPNX16C(x1, 4); a16 = BBE_REPNX16C(x1, 5);
    a07 = BBE_REPNX16C(x1, 6); a17 = BBE_REPNX16C(x1, 7);

    w0 = BBE_MULRNX16C(tw0, a00, vsa0);
    w1 = BBE_MULRNX16C(tw0, a10, vsa0);

    BBE_MULANX16C(w0, a01, tw1); BBE_MULANX16C(w1, a11, tw1);
    BBE_MULANX16C(w0, a02, tw2); BBE_MULANX16C(w1, a12, tw2);
    BBE_MULANX16C(w0, a03, tw3); BBE_MULANX16C(w1, a13, tw3);
    BBE_MULANX16C(w0, a04, tw4); BBE_MULANX16C(w1, a14, tw4);
    BBE_MULANX16C(w0, a05, tw5); BBE_MULANX16C(w1, a15, tw5);
    BBE_MULANX16C(w0, a06, tw6); BBE_MULANX16C(w1, a16, tw6);
    BBE_MULANX16C(w0, a07, tw7); BBE_MULANX16C(w1, a17, tw7);

    a0 = BBE_PACKVNX40(w0, vsa0);
    a1 = BBE_PACKVNX40(w1, vsa0);

    a1 = BBE_MULNX16CPACKQ(a1, tw_T16_8);

    p0 = BBE_ADDNX16(a0, a1);
    p1 = BBE_SUBNX16(a0, a1);
  }
 
  //----------------------------------------------------------------------------
  // Real-to-complex spectrum converter.

  #else
  bexp = cfftas16(y, x, bexp);
  p0 = BBE_LVNX16_I(Y, 0 * (2 * BBE_SIMD_WIDTH));
  p1 = BBE_LVNX16_I(Y, 1 * (2 * BBE_SIMD_WIDTH));
  shift0 = bexp;
  #endif
  {
    xb_vecNx16 a0, a1, a2;
    xb_vecNx40 b0, b1, b2;
    xb_vecNx16 c0, c1, c2;

    //
    // 2nd scaling stage.
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(p0);
    BBE_RANGENX16(p1);

    bexp = BBE_RRANGE();

    // Data scaling shift motivation:
    //  1. Data should be normalized using the common block exponent.
    //  2. Forward spectrum conversion doesn't change data magnitude, but includes
    //     twiddle multiplication, which requires a reservation of 1 bit position.
    shift1 = 1 - bexp;

    // Setup for Q30->Q15 conversion coupled with data scaling
    vsa0 = BBE_MOVVSA32(15 + shift1);

    //
    // Spectrum conversion. MATLAB code:
    //  N = 32;
    //  a0 = X(1:N/4);
    //  a1 = X(N/4+1);
    //  a2 = [X(1);X(N/2-(0:N/4-2))];
    //  twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
    //  y = [a0.*(1-twd)/2+conj(a2).*(1+twd)/2; ...
    //       conj(a1); ...
    //       wrev(conj(a0).*(1+conj(twd))/2+a2.*(1-conj(twd))/2)];
    //

    // a0 = X(N/4-(SIMD_WIDTH_2:-1:1));
    a0 = p0;
    // a1 = X(N/4+1)
    a1 = p1;
    // a2 = X(N/4+1+(SIMD_WIDTH_2:-1:1));
    a2 = BBE_SELNX16(p0, p1, sel_rsft);

    // tw = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // twc0 = (1-tw)/2; CQ15
    // twc1 = (1+tw)/2; CQ15
    // twc2 = (1+conj(tw))/2; CQ15
    twc2 = BBE_CONJSNX16C(twc1);

    // b0 = a0*(1-tw)/2; Q30 <- Q15*Q15
    b0 = BBE_MULNX16C(twc0, a0);
    // b1 = conj(a1);
    b1 = BBE_UNPKQNX16(a1);
    b1 = BBE_CONJNX40C(b1);
    // b2 = conj(a0)*(1+conj(tw))/2; Q30 <- Q15*Q15
    b2 = BBE_MULNX16J(twc2, a0);

    // b0 += conj(a2)*(1+tw)/2; Q30 <- Q15*Q15
    BBE_MULANX16J(b0, twc1, a2);
    // b2 += a2*(1-conj(tw))/2; Q30 <- Q15*Q15
    BBE_MULANX16J(b2, a2, twc0);

    b0 = BBE_RNDADJNX40(b0, vsa0);
    b1 = BBE_RNDADJNX40(b1, vsa0);
    b2 = BBE_RNDADJNX40(b2, vsa0);

    // Q15 <- Q30/2^shift - 15 w/ rounding and saturation
    c0 = BBE_PACKVNX40(b0, vsa0);
    c1 = BBE_PACKVNX40(b1, vsa0);
    c2 = BBE_PACKVNX40(b2, vsa0);

    // Right-hand vector indices: N/4+1+n*SIMD_WIDTH_2+(1:SIMD_WIDTH_2)
    c2 = BBE_SHFLNX16I(c2, BBE_SHFLI_REVERSE_2);

    Y_va = BBE_ZALIGN();

    BBE_SAVNX16_XP(c0, Y_va, Y, 8 * 4);
    BBE_SAVNX16_XP(c1, Y_va, Y, 1 * 4);
    BBE_SAVNX16_XP(c2, Y_va, Y, 8 * 4);

    BBE_SAVNX16POS_FP(Y_va, Y);
  }

  return (shift0 + shift1);

} /* rfft32() */
#endif

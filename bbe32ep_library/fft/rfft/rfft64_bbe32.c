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

#if !(HAVE_FFT && 1)
DISCARD_FUN(int, rfft64, (complex_fract16 * restrict y, int16_t * restrict x, int bexp))
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

int rfft64 ( complex_fract16 * restrict y, int16_t * restrict x, int bexp )
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
  xb_vecNx16 s0, s1, s2, s3, s4;

  vselN sel_rsft;
  vsaN  vsa0;

  int shift, shiftSum;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file.

  C = (const xb_vecNx16*)fft_tw_tab_rfft_64;

  // T32_8, 8..31
  BBE_LVNX16_IP(tw01, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw02, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw03, C, +4 * BBE_SIMD_WIDTH / 2);
  // T8_2 x I4, 8..31
  BBE_LVNX16_IP(tw11, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw12, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw13, C, +4 * BBE_SIMD_WIDTH / 2);
  // N=64; reshape(1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2)),8,N/32).'
  BBE_LVNX16_IP(twc10, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(twc11, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(twc00, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(twc01, C, +4 * BBE_SIMD_WIDTH / 2);
  // Shift-and-reverse pattern for select.
  BBE_LVNX16_IP(p0, C, +2 * BBE_SIMD_WIDTH);

  sel_rsft = BBE_MOVVSV(p0, 0);

  //----------------------------------------------------------------------------
  // Complex-valued 32-point FFT:
  //  DFT32 = ( DFT2 x I16 )*( L8_2 x I4 )*
  //          ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )*
  //          T32_8*( DFT4 x I8 )

  //
  // Stage 1: T32_8*( DFT4 x I8 )
  //

  X = (const xb_vecNx16*)x;

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    a0 = BBE_LVNX16_I(X, +0 * (2 * BBE_SIMD_WIDTH));
    a1 = BBE_LVNX16_I(X, +1 * (2 * BBE_SIMD_WIDTH));
    a2 = BBE_LVNX16_I(X, +2 * (2 * BBE_SIMD_WIDTH));
    a3 = BBE_LVNX16_I(X, +3 * (2 * BBE_SIMD_WIDTH));

    //
    // 1st scaling stage.
    //

    vsa0 = BBE_MOVVSA32(bexp);

    a0 = BBE_SLLNX16(a0, vsa0);
    a1 = BBE_SLLNX16(a1, vsa0);
    a2 = BBE_SLLNX16(a2, vsa0);
    a3 = BBE_SLLNX16(a3, vsa0);

    // Radix-4 w/ twiddle multiplication
    shiftSum = 3 - bexp;

    BBE_FFTWMODE(0x10 | 3);

    //
    // DFT4 x I8
    //

    BBE_MOVSAV(a2);
    BBE_MOVSBV(a3);

    b0 = BBE_FFTADD4SABNX16(a0, a1, 0, 0);
    b1 = BBE_FFTADD4SABNX16(a0, a1, 1, 0);
    b2 = BBE_FFTADD4SABNX16(a0, a1, 2, 0);
    b3 = BBE_FFTADD4SABNX16(a0, a1, 3, 0);

    //
    // T32_8
    //

    p0 = b0;
    p1 = BBE_MULNX16CPACKQ(b1, tw01);
    p2 = BBE_MULNX16CPACKQ(b2, tw02);
    p3 = BBE_MULNX16CPACKQ(b3, tw03);
  }

  //
  // Stage 2: ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )
  //

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //
    // 2nd scaling stage
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(p0);
    BBE_RANGENX16(p1);
    BBE_RANGENX16(p2);
    BBE_RANGENX16(p3);

    bexp = BBE_RRANGE();

    // Radix-4 w/ twiddle multiplication
    shift = 3 - bexp;

    XT_MOVLTZ(shift, 0, shift);

    shiftSum += shift;

    BBE_FFTWMODE(0x10 | shift);

    // I2 x L16_8
    BBE_DSELNX16I(a1, a0, p1, p0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(a3, a2, p3, p2, BBE_DSELI_INTERLEAVE_2);

    // L16_8 x I2
    BBE_DSELNX16I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4);

    //
    // DFT4 x I8
    //

    BBE_MOVSAV(b2);
    BBE_MOVSBV(b3);

    a0 = BBE_FFTADD4SABNX16(b0, b1, 0, 0);
    a1 = BBE_FFTADD4SABNX16(b0, b1, 1, 0);
    a2 = BBE_FFTADD4SABNX16(b0, b1, 2, 0);
    a3 = BBE_FFTADD4SABNX16(b0, b1, 3, 0);

    //
    // T8_2 x I4
    //

    q0 = a0;
    q1 = BBE_MULNX16CPACKQ(a1, tw11);
    q2 = BBE_MULNX16CPACKQ(a2, tw12);
    q3 = BBE_MULNX16CPACKQ(a3, tw13);
  }

  //
  // Stage 3: ( DFT2 x I16 )*( L8_2 x I4 )
  //

  {
    xb_vecNx16 a0, a1, a2, a3;

    //
    // 3rd scaling stage
    //

    BBE_WRANGE(4);

    BBE_RANGENX16(q0);
    BBE_RANGENX16(q1);
    BBE_RANGENX16(q2);
    BBE_RANGENX16(q3);

    bexp = BBE_RRANGE();

    // Data scaling shift motivation:
    //  1. Radix-2 requires 1 bit position.
    //  2. Forward spectrum conversion doesn't change data magnitude, but includes
    //     twiddle multiplication, which requires additional bit position.
    shift = 2 - bexp;

    XT_MOVLTZ(shift, 0, shift);

    shiftSum += shift;

    BBE_FFTWMODE(shift);

    // L8_2 x I4
    a0 = BBE_SELNX16I(q1, q0, BBE_SELI_EXTRACT_LO_HALVES);
    a1 = BBE_SELNX16I(q3, q2, BBE_SELI_EXTRACT_LO_HALVES);
    a2 = BBE_SELNX16I(q1, q0, BBE_SELI_EXTRACT_HI_HALVES);
    a3 = BBE_SELNX16I(q3, q2, BBE_SELI_EXTRACT_HI_HALVES);

    BBE_MOVSAV(0);
    BBE_MOVSBV(0);

    // DFT2 x I16. We use a radix-4 instruction just because it supports shift
    // amounts exceeding 1.
    r0 = BBE_FFTADD4SABNX16(a0, a2, 0, 0);
    r1 = BBE_FFTADD4SABNX16(a1, a3, 0, 0);
    r2 = BBE_FFTADD4SABNX16(a0, a2, 2, 0);
    r3 = BBE_FFTADD4SABNX16(a1, a3, 2, 0);
  }

  //----------------------------------------------------------------------------
  // Real-to-complex spectrum converter. MATLAB code:
  //  N = 64;
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

  // Transform the central sample:
  //  y(N/4+1) = conj(X(N/4+1));
  s2 = BBE_CONJSNX16C(r2);

  //
  // Manually unrolled loop, iteration n=0.
  //

  {
    xb_vecNx16 a0, a1;
    xb_vecNx40 b0, b1;

    // a0 = X(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1));
    a0 = r1;
    // a1 = X(N/4+1+n*SIMD_WIDTH_2+(SIMD_WIDTH_2:-1:1)); 
    a1 = BBE_SELNX16(r3, r2, sel_rsft);

    // tw = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // twc0 = (1-tw)/2; Q15
    // twc1 = (1+tw)/2; Q15
    // twc2 = (1+conj(tw))/2;
    twc02 = BBE_CONJSNX16C(twc01);

    // b0 = a0*(1-tw)/2; CQ30 <- CQ15*CQ15
    b0 = BBE_MULNX16C(twc00, a0);
    // b1 = conj(a0)*(1+conj(tw))/2; CQ30 <- CQ15*CQ15
    b1 = BBE_MULNX16J(twc02, a0);

    // b0 += conj(a1)*(1+tw)/2;
    BBE_MULANX16J(b0, twc01, a1);
    // b1 += a1*(1-conj(tw))/2;
    BBE_MULANX16J(b1, a1, twc00);

    // CQ15 <- CQ30 - 15 w/ rounding and saturation
    s1 = BBE_PACKQNX40(b0);
    s3 = BBE_PACKQNX40(b1);
  }

  //
  // Iteration n=1.
  //

  {
    xb_vecNx16 a0, a1;
    xb_vecNx40 b0, b1;

    // a0 = X(N/4+1-n*SIMD_WIDTH_2-(SIMD_WIDTH_2:-1:1));
    a0 = r0;
    // a1 = X(N/4+1+n*SIMD_WIDTH_2+(SIMD_WIDTH_2:-1:1)); 
    a1 = BBE_SELNX16(r0, r3, sel_rsft);

    // tw = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2))
    // twc0 = (1-tw)/2; Q15
    // twc1 = (1+tw)/2; Q15
    // twc2 = (1+conj(tw))/2;
    twc12 = BBE_CONJSNX16C(twc11);

    // b0 = a0*(1-tw)/2; CQ30 <- CQ15*CQ15
    b0 = BBE_MULNX16C(twc10, a0);
    // b1 = conj(a0)*(1+conj(tw))/2; CQ30 <- CQ15*CQ15
    b1 = BBE_MULNX16J(twc12, a0);

    // b0 += conj(a1)*(1+tw)/2;
    BBE_MULANX16J(b0, twc11, a1);
    // b1 += a1*(1-conj(tw))/2;
    BBE_MULANX16J(b1, a1, twc10);

    // CQ15 <- CQ30 - 15 w/ rounding and saturation
    s0 = BBE_PACKQNX40(b0);
    s4 = BBE_PACKQNX40(b1);
  }

  //
  // Store results to the output array. The right half of samples
  // is to be stored in reverse order.
  //

  Y = (xb_vecNx16*)y;

  {
    xb_vecNx16 a0, a1, a2, a3, a4;

    a0 = s0;
    a1 = s1;
    // s2 contains only y(N/4+1), at position 0.
    a2 = BBE_SELNX16(s2, s3, sel_rsft);
    a3 = BBE_SELNX16(s3, s4, sel_rsft);
    // Nyquist frequency bin y(N/2+1).
    a4 = s4;

    BBE_SVNX16_I(a0, Y, 0 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_I(a1, Y, 1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_I(a2, Y, 2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_I(a3, Y, 3 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SPNX16_I(a4, Y, 4 * 4 * BBE_SIMD_WIDTH / 2);
  }

  return (shiftSum);
} /* rfft64() */
#endif

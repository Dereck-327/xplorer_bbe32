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
DISCARD_FUN(int, rifft16, (int16_t * restrict y, complex_fract16 * restrict x, int bexp))
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

int rifft16 ( int16_t * restrict y, complex_fract16 * restrict x, int bexp )
{
  //----------------------------------------------------------------------------
  // MATLAB code:
  //  N = 16;
  //  % Complex-to-real spectrum conversion
  //  twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
  //  tw0 = [1-conj(twd);1;wrev(1+twd(2:N/4))];
  //  tw1 = [1+conj(twd);1;wrev(1-twd(2:N/4))];
  //  a0 = [x(1:N/4);wrev(conj(x(2:N/4+1)))];
  //  a1 = [wrev(conj(x(N/4+1+(0:N/4))));x(N/4+1+(1:N/4-1))];
  //  x = a0.*tw0+a1.*tw1;
  //  % Inverse complex-valued DFT of half the original size
  //  DFTN = conj(fft(eye(N/2)));
  //  y = DFTN*x;
  //  y = reshape([real(y(:))';imag(y(:))'],N,1);
  //

  const xb_vecNx16 * TW = (const xb_vecNx16*)fft_tw_tab_rifft_16;

  xb_vecNx16 x0, x1;
  xb_vecNx40 w0;

  xb_vecNx16 tw0, tw1, tw2, tw3;
  xb_vecNx16 tw4, tw5, tw6, tw7;

  xb_vecNx16 p0, p1, p2, p3;
  xb_vecNx16 p4, p5, p6, p7;

  vsaN vsa0;

  const int scale0 = 2; // Spectrum conversion
  const int scale1 = 3; // Radix-8 DFT w/o twiddles

  int shift0, shift1, nsa;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Spectrum conversion

  {
    xb_vecNx16 a0, a1;

    vselN s0, s1, s2;

    // Select pattern to store Nyquist frequency bin value to the imaginary part
    // of 0-th bin.
    BBE_LVNX16_IP(p0, TW, +2 * BBE_SIMD_WIDTH);
    // a0 sequence select pattern
    BBE_LVNX16_IP(p1, TW, +2 * BBE_SIMD_WIDTH);
    // a1 sequence select pattern
    BBE_LVNX16_IP(p2, TW, +2 * BBE_SIMD_WIDTH);

    // tw0 = [1-conj(twd);1;wrev(1+twd(2:N/4))]; CQ14
    BBE_LVNX16_IP(tw0, TW, +2 * BBE_SIMD_WIDTH);
    // tw1 = [1+conj(twd);1;wrev(1-twd(2:N/4))]; CQ14
    BBE_LVNX16_IP(tw1, TW, +2 * BBE_SIMD_WIDTH);

    s0 = BBE_MOVVSV(p0, 0);
    s1 = BBE_MOVVSV(p1, 0);
    s2 = BBE_MOVVSV(p2, 0);

    x0 = BBE_LVNX16_I((const xb_vecNx16*)x, 0 * 4 * BBE_SIMD_WIDTH / 2);
    x1 = BBE_LPNX16_I((const xb_vecNx16*)x, 1 * 4 * BBE_SIMD_WIDTH / 2);

    p0 = BBE_SELNX16(x1, x0, s0);
    p1 = BBE_CONJSNX16C(x0);

    a0 = BBE_SELNX16(p1, x0, s1);
    a1 = BBE_SELNX16(p1, p0, s2);

    shift0 = scale0 - bexp;

    vsa0 = BBE_MOVVSA32(14 + shift0);

    // CQ29 <- CQ14*CQ15
    w0 = BBE_MULRNX16C(a0, tw0, vsa0);

    // CQ29 <- CQ29 + CQ14*CQ15
    BBE_MULANX16C(w0, a1, tw1);

    // CQ(15+bexp-scale0) <- CQ29 - 14 + bexp - scale0 w/ rounding, saturation
    x0 = BBE_PACKVNX40(w0, vsa0);
  }

  //----------------------------------------------------------------------------
  // IDFT8

  {
    xb_vecNx16 y0;

    BBE_LVNX16_IP(tw1, TW, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(tw2, TW, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(tw3, TW, +4 * BBE_SIMD_WIDTH / 2);

    BBE_LVNX16_IP(tw4, TW, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(tw5, TW, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(tw6, TW, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(tw7, TW, +4 * BBE_SIMD_WIDTH / 2);

    p0 = BBE_REPNX16C(x0, 0);
    p1 = BBE_REPNX16C(x0, 1);
    p2 = BBE_REPNX16C(x0, 2);
    p3 = BBE_REPNX16C(x0, 3);

    p4 = BBE_REPNX16C(x0, 4);
    p5 = BBE_REPNX16C(x0, 5);
    p6 = BBE_REPNX16C(x0, 6);
    p7 = BBE_REPNX16C(x0, 7);

    // CQ(30+bexp-scale0) <- CQ(15+bexp-scale0) + 15
    w0 = BBE_UNPKQNX16(p0);

    // CQ(30+bexp-scale0) <- CQ(30+bexp-scale0) + CQ15*CQ(15+bexp-scale0)
    BBE_MULANX16C(w0, tw1, p1);
    BBE_MULANX16C(w0, tw2, p2);
    BBE_MULANX16C(w0, tw3, p3);

    // CQ(30+bexp-scale0) <- CQ(30+bexp-scale0) + CQ15*CQ(15+bexp-scale0)
    BBE_MULANX16C(w0, tw4, p4);
    BBE_MULANX16C(w0, tw5, p5);
    BBE_MULANX16C(w0, tw6, p6);
    BBE_MULANX16C(w0, tw7, p7);

    BBE_WRANGE(4);

    BBE_RANGENX16(x0);

    // Normalization shift amount for the converted spectrum.
    nsa = BBE_RRANGE();

    shift1 = scale1 - nsa;

    vsa0 = BBE_MOVVSA32(15 + shift1);

    // Rounding is applied as late as possible to allow more parallelism for
    // the compiler.
    w0 = BBE_RNDADJNX40(w0, vsa0);

    // CQ(15+bexp+nsa-scale0-scale1) <- CQ(30+bexp-scale0) - 15 + nsa - scale1
    y0 = BBE_PACKVNX40(w0, vsa0);

    BBE_SVNX16_I(y0, (xb_vecNx16*)y, 0);
  }

  return (shift0 + shift1);
} /* rifft16() */
#endif

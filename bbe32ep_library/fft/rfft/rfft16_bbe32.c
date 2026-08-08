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
DISCARD_FUN(int, rfft16, (complex_fract16 * restrict y, int16_t * restrict x, int bexp))
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

int rfft16 ( complex_fract16 * restrict y, int16_t * restrict x, int bexp )
{
  const xb_vecNx16 * TW = (const xb_vecNx16*)fft_tw_tab_rfft_16;

  xb_vecNx16 tw00, tw01, tw02, tw03, tw04, tw05, tw06, tw07;
  xb_vecNx16 tw10, tw11, tw12, tw13, tw14, tw15, tw16, tw17, tw20;

  unsigned int x00, x01, x02, x03;
  unsigned int x10, x11, x12, x13;

  xb_vecNx40 w0, w1;
  xb_vecNx16 x0;
  xb_vecNx16 y0, y1;

  xb_int40 z0;

  int shift;

  vsaN vsa0;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //
  // For spectrum bins 0..7 we use a direct multiplication to the upper half of
  // the DFT16 matrix.
  //

  BBE_LVNX16_IP(tw00, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw01, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw02, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw03, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw04, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw05, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw06, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw07, TW, 4 * BBE_SIMD_WIDTH / 2);

  BBE_LVNX16_IP(tw10, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw11, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw12, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw13, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw14, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw15, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw16, TW, 4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw17, TW, 4 * BBE_SIMD_WIDTH / 2);

  x0 = BBE_LVNX16_I((const xb_vecNx16*)x, 0);

  x00 = BBE_EXTRNX16C(x0, 0);
  x01 = BBE_EXTRNX16C(x0, 1);
  x02 = BBE_EXTRNX16C(x0, 2);
  x03 = BBE_EXTRNX16C(x0, 3);

  x10 = BBE_EXTRNX16C(x0, 4);
  x11 = BBE_EXTRNX16C(x0, 5);
  x12 = BBE_EXTRNX16C(x0, 6);
  x13 = BBE_EXTRNX16C(x0, 7);

  // Select scaling shift for radix-16 w/o twiddles.
  shift = 4 - bexp;

  vsa0 = BBE_MOVVSA32(15 + shift);

  // Each MULNXPR takes a sum of two weighted columns of 8x16. Columns
  // are complex, while weights (signal samples) are real.
  w0 = BBE_MULRNX16PR(tw01, tw00, x00, vsa0);

  BBE_MULANX16PR(w0, tw03, tw02, x01);
  BBE_MULANX16PR(w0, tw05, tw04, x02);
  BBE_MULANX16PR(w0, tw07, tw06, x03);

  BBE_MULANX16PR(w0, tw11, tw10, x10);
  BBE_MULANX16PR(w0, tw13, tw12, x11);
  BBE_MULANX16PR(w0, tw15, tw14, x12);
  BBE_MULANX16PR(w0, tw17, tw16, x13);

  //
  // For the Nyquist frequency bin all we need is to subtract the sum of all
  // odd-numbered samples from the sum of even samples.
  //

  tw20 = BBE_MOVVA16C(0x80007fff);

  w1 = BBE_MULNX16(x0, tw20);
  z0 = BBE_RADDNX40(w1);
  w1 = BBE_MOVNX40_FROM40(z0);
  w1 = BBE_RNDADJNX40(w1, vsa0);

  // CQ(15-(4-bexp)) <- CQ30 + bexp - 4
  y0 = BBE_PACKVNX40(w0, vsa0);
  y1 = BBE_PACKVNX40(w1, vsa0);

  BBE_SVNX16_I(y0, (xb_vecNx16*)y, 0 * 4 * BBE_SIMD_WIDTH / 2);
  BBE_SPNX16_I(y1, (xb_vecNx16*)y, 1 * 4 * BBE_SIMD_WIDTH / 2);

  return (shift);
} /* rfft16() */
#endif

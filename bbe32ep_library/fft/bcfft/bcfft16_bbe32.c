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
    Blockwise radix-2 forward FFT on complex data, no data scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Twiddle tables declarations. */
#include "fft_tw.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"

/*-------------------------------------------------------------------------
Blockwise radix-2 forward FFT on complex data, no data scaling

Description: These functions make forward FFT on L blocks, each of N=2^n
complex samples, where n=4..7. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input signal
  Output:          
    y[L][N]   Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

void bcfft16 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          C;

  xb_vecNx16 tw_T4_2, tw_T16_4_0, tw_T16_4_1;
  xb_vecNx16 t0, t1;

  vselN sel_L16_4_0, sel_L16_4_1;

  int l;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  NASSERT_ALIGN32(fft16_tw);

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file

  C = (const xb_vecNx16*)fft16_tw;

  BBE_LVNX16_IP(tw_T4_2, C, +2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw_T16_4_0, C, +2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw_T16_4_1, C, +2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(t0, C, +2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(t1, C, +2 * BBE_SIMD_WIDTH);

  sel_L16_4_0 = BBE_MOVVSELNX16(t0, 0);
  sel_L16_4_1 = BBE_MOVVSELNX16(t1, 0);

  //----------------------------------------------------------------------------
  // Apply DFT16 to each of L data blocks. DFT16 decomposition (basically, DIF
  // with trivial modifications):
  //   DFT16 -> ( DFT2 x I8 )*( L4_2 x I4 )*( T4_2 x I4 )*( DFT2 x I8 )*L16_4*
  //            T16_4*( DFT2 x I8 )*( L4_2 x I4 )*( T4_2 x I4 )*( DFT2 x I8 )

  Y = (xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;

  for (l = 0; l<L; l++)
  {
    xb_vecNx16 a0, a1;
    xb_vecNx16 b0, b1;

    BBE_LVNX16_IP(a0, X, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(a1, X, +4 * BBE_SIMD_WIDTH / 2);

    // DFT2 x I8
    b0 = BBE_ADDNX16(a0, a1);
    b1 = BBE_SUBNX16(a0, a1);

    // T4_2 x I4
    // tw_T4_2 contains only trivial twiddle factors of 1 and -1j.
    b1 = BBE_MULNX16CPACKQ(b1, tw_T4_2);

    // L4_2 x I4
    a0 = BBE_SELNX16I(b1, b0, BBE_SELI_EXTRACT_LO_HALVES);
    a1 = BBE_SELNX16I(b1, b0, BBE_SELI_EXTRACT_HI_HALVES);

    // DFT2 x I8
    b0 = BBE_ADDNX16(a0, a1);
    b1 = BBE_SUBNX16(a0, a1);

    // T16_4
    b0 = BBE_MULNX16CPACKQ(b0, tw_T16_4_0);
    b1 = BBE_MULNX16CPACKQ(b1, tw_T16_4_1);

    // L16_4
    a0 = BBE_SELNX16(b1, b0, sel_L16_4_0);
    a1 = BBE_SELNX16(b1, b0, sel_L16_4_1);

    // DFT2 x I8
    b0 = BBE_ADDNX16(a0, a1);
    b1 = BBE_SUBNX16(a0, a1);

    // T4_2 x I4
    b1 = BBE_MULNX16CPACKQ(b1, tw_T4_2);

    // L4_2 x I4
    a0 = BBE_SELNX16I(b1, b0, BBE_SELI_EXTRACT_LO_HALVES);
    a1 = BBE_SELNX16I(b1, b0, BBE_SELI_EXTRACT_HI_HALVES);

    // DFT2 x I8
    b0 = BBE_ADDNX16(a0, a1);
    b1 = BBE_SUBNX16(a0, a1);

    BBE_SVNX16_IP(b0, Y, +4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b1, Y, +4 * BBE_SIMD_WIDTH / 2);
  }
} /* bcfft16() */

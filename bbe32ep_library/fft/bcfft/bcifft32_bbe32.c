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
    Blockwise radix-2 inverse FFT on complex data, no data scaling
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
DISCARD_FUN(void, bcifft32, (complex_fract16 * restrict y, complex_fract16 * restrict x, int L))
#else
/*-------------------------------------------------------------------------
Blockwise radix-2 inverse FFT on complex data, no data scaling

Description: These functions make inverse FFT on L blocks, each of N=2^n
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

void bcifft32 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X0;
  const xb_vecNx16 *          X1;
  const xb_vecNx16 *          C;

  xb_vecNx16 tw_T32_8_1, tw_T32_8_2, tw_T32_8_3;
  xb_vecNx16 tw_T8_2_1, tw_T8_2_2, tw_T8_2_3;

  int l;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  NASSERT_ALIGN32(fft32_tw);

  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file

  C = (const xb_vecNx16*)fft32_tw;

  BBE_LVNX16_IP(tw_T32_8_1, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw_T32_8_2, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw_T32_8_3, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw_T8_2_1, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw_T8_2_2, C, +4 * BBE_SIMD_WIDTH / 2);
  BBE_LVNX16_IP(tw_T8_2_3, C, +4 * BBE_SIMD_WIDTH / 2);

  BBE_FFTWMODE(0x10);

  //----------------------------------------------------------------------------
  // Apply DFT32 to each of L data blocks. DFT32 decomposition:
  //   DFT32 -> ( DFT2 x I16 )*( L8_2 x I4 )*
  //            ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )*
  //            T32_8*( DFT4 x I8 )

  Y = (xb_vecNx16*)y;
  X0 = (const xb_vecNx16*)((uintptr_t)x + 0 * 4 * BBE_SIMD_WIDTH / 2);
  X1 = (const xb_vecNx16*)((uintptr_t)x + 2 * 4 * BBE_SIMD_WIDTH / 2);

  for (l = 0; l<L; l++)
  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    //********************************************************************
    // Stage 1: T32_8*( DFT4 x I8 )

    BBE_LVNX16_IP(a0, X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(a1, X0, 3 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_LVA_IP(X1, 1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_LVB_IP(X1, 3 * 4 * BBE_SIMD_WIDTH / 2);

    //
    // DFT4 x I8
    //

    b0 = BBE_FFTADD4SABNX16(a0, a1, 0, 0);
    b3 = BBE_FFTADD4SABNX16(a0, a1, 1, 0);
    b2 = BBE_FFTADD4SABNX16(a0, a1, 2, 0);
    b1 = BBE_FFTADD4SABNX16(a0, a1, 3, 0);

    //
    // T32_8
    //

    b1 = BBE_MULNX16JPACKQ(b1, tw_T32_8_1);
    b2 = BBE_MULNX16JPACKQ(b2, tw_T32_8_2);
    b3 = BBE_MULNX16JPACKQ(b3, tw_T32_8_3);

    //********************************************************************
    // Stage 2: ( T8_2 x I4 )*( DFT4 x I8 )*( L16_8 x I2 )*( I2 x L16_8 )

    // I2 x L16_8
    BBE_DSELNX16I(a1, a0, b1, b0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(a3, a2, b3, b2, BBE_DSELI_INTERLEAVE_2);

    // L16_8 x I2
    BBE_DSELNX16I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4);

    //
    // DFT4 x I8
    //

    BBE_MOVSAV(b2);
    BBE_MOVSBV(b3);

    a0 = BBE_FFTADD4SABNX16(b0, b1, 0, 0);
    a3 = BBE_FFTADD4SABNX16(b0, b1, 1, 0);
    a2 = BBE_FFTADD4SABNX16(b0, b1, 2, 0);
    a1 = BBE_FFTADD4SABNX16(b0, b1, 3, 0);

    //
    // T8_2 x I4
    //

    a1 = BBE_MULNX16JPACKQ(a1, tw_T8_2_1);
    a2 = BBE_MULNX16JPACKQ(a2, tw_T8_2_2);
    a3 = BBE_MULNX16JPACKQ(a3, tw_T8_2_3);

    //********************************************************************
    // Stage 3: ( DFT2 x I16 )*( L8_2 x I4 )

    // L8_2 x I4
    b0 = BBE_SELNX16I(a1, a0, BBE_SELI_EXTRACT_LO_HALVES);
    b1 = BBE_SELNX16I(a3, a2, BBE_SELI_EXTRACT_LO_HALVES);
    b2 = BBE_SELNX16I(a1, a0, BBE_SELI_EXTRACT_HI_HALVES);
    b3 = BBE_SELNX16I(a3, a2, BBE_SELI_EXTRACT_HI_HALVES);

    // DFT2 x I16
    a0 = BBE_ADDNX16(b0, b2);
    a1 = BBE_ADDNX16(b1, b3);
    a2 = BBE_SUBNX16(b0, b2);
    a3 = BBE_SUBNX16(b1, b3);

    BBE_SVNX16_IP(a0, Y, +4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(a1, Y, +4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(a2, Y, +4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(a3, Y, +4 * BBE_SIMD_WIDTH / 2);
  }
} /* bcifft32() */
#endif

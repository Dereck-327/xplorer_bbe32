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
DISCARD_FUN(void, bcifft64, (complex_fract16 * restrict y, complex_fract16 * restrict x, int L))
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

void bcifft64 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          TW0;
  const xb_vecNx16 *          TW1;

  xb_vecNx16 T64_16_10, T64_16_11, T64_16_20;
  xb_vecNx16 T64_16_21, T64_16_30, T64_16_31;

  xb_vecNx16 T16_4_10, T16_4_11, T16_4_20;
  xb_vecNx16 T16_4_21, T16_4_30, T16_4_31;

  int l;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  BBE_FFTWMODE(0x10);

  //----------------------------------------------------------------------------
  // Apply DFT64 to each of L data blocks. DFT64 decomposition:
  //   DFT64 -> ( DFT4 x I16 )*( L16_4 x I4 )*
  //            ( T16_4 x I4 )*( DFT4 x I16 )*
  //            ( I2 x L16_8 x I2 )*( I4 x L16_8 )*( L8_2 x I8 )*
  //            T64_16*( DFT4 x I16 )

  Y = (xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;

  TW0 = (const xb_vecNx16*)fft64_tw1;
  TW1 = (const xb_vecNx16*)fft64_tw2;

  __Pragma("ymemory( X )");
  for (l = 0; l<L; l++)
  {
    xb_vecNx16 a00, a01, a02, a03;
    xb_vecNx16 a10, a11, a12, a13;
    xb_vecNx16 a20, a21, a30, a31;
    xb_vecNx16 b00, b01, b02, b03;
    xb_vecNx16 b10, b11, b12, b13;
    xb_vecNx16 b20, b21, b30, b31;

    xb_vecNx16 tw0, tw1;

    //
    // Reload twiddle factor table on each iteration.
    //

    T64_16_10 = BBE_LVNX16_I(TW0, 0 * 4 * BBE_SIMD_WIDTH / 2);
    T64_16_20 = BBE_LVNX16_I(TW0, 1 * 4 * BBE_SIMD_WIDTH / 2);
    T64_16_30 = BBE_LVNX16_I(TW0, 2 * 4 * BBE_SIMD_WIDTH / 2);
    T64_16_11 = BBE_LVNX16_I(TW0, 3 * 4 * BBE_SIMD_WIDTH / 2);
    T64_16_21 = BBE_LVNX16_I(TW0, 4 * 4 * BBE_SIMD_WIDTH / 2);
    T64_16_31 = BBE_LVNX16_I(TW0, 5 * 4 * BBE_SIMD_WIDTH / 2);

    tw0 = BBE_LVNX16_I(TW1, 0 * 4 * BBE_SIMD_WIDTH / 2);
    tw1 = BBE_LVNX16_I(TW1, 1 * 4 * BBE_SIMD_WIDTH / 2);

    T16_4_10 = BBE_SHFLNX16I(tw0, BBE_SHFLI_REP_2X4_OFFSET_0);
    T16_4_20 = BBE_SHFLNX16I(tw0, BBE_SHFLI_REP_2X4_OFFSET_1);
    T16_4_30 = BBE_SHFLNX16I(tw0, BBE_SHFLI_REP_2X4_OFFSET_2);
    T16_4_11 = BBE_SHFLNX16I(tw0, BBE_SHFLI_REP_2X4_OFFSET_3);
    T16_4_21 = BBE_SHFLNX16I(tw1, BBE_SHFLI_REP_2X4_OFFSET_0);
    T16_4_31 = BBE_SHFLNX16I(tw1, BBE_SHFLI_REP_2X4_OFFSET_1);

    //**************************************************************************
    // Stage 1: T64_16*( DFT4 x I16 )

    BBE_LVNX16_IP(a00, X, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(a01, X, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(a10, X, +4 * BBE_SIMD_WIDTH / 2);
    BBE_LVNX16_IP(a11, X, +4 * BBE_SIMD_WIDTH / 2);

    BBE_LVA_IP(X, +4 * BBE_SIMD_WIDTH / 2); // a20
    BBE_LVC_IP(X, +4 * BBE_SIMD_WIDTH / 2); // a21
    BBE_LVB_IP(X, +4 * BBE_SIMD_WIDTH / 2); // a30
    BBE_LVD_IP(X, +4 * BBE_SIMD_WIDTH / 2); // a31

    //
    // DFT4 x I16
    //

    b00 = BBE_FFTADD4SABNX16(a00, a10, 0, 0);
    b01 = BBE_FFTADD4SCDNX16(a01, a11, 0, 0);
    b30 = BBE_FFTADD4SABNX16(a00, a10, 1, 0);
    b31 = BBE_FFTADD4SCDNX16(a01, a11, 1, 0);
    b20 = BBE_FFTADD4SABNX16(a00, a10, 2, 0);
    b21 = BBE_FFTADD4SCDNX16(a01, a11, 2, 0);
    b10 = BBE_FFTADD4SABNX16(a00, a10, 3, 0);
    b11 = BBE_FFTADD4SCDNX16(a01, a11, 3, 0);

    //
    // T64_16
    //

    b10 = BBE_MULNX16JPACKQ(b10, T64_16_10);
    b11 = BBE_MULNX16JPACKQ(b11, T64_16_11);
    b20 = BBE_MULNX16JPACKQ(b20, T64_16_20);
    b21 = BBE_MULNX16JPACKQ(b21, T64_16_21);
    b30 = BBE_MULNX16JPACKQ(b30, T64_16_30);
    b31 = BBE_MULNX16JPACKQ(b31, T64_16_31);

    //**************************************************************************
    // Stage 2: ( T16_4 x I4 )*( DFT4 x I16 )*( I2 x L16_8 x I2 )*
    //          ( I4 x L16_8 )*( L8_2 x I8 )

    //
    // L8_2 x I8
    //

    a00 = b00; a01 = b10; a02 = b20; a03 = b30;
    a10 = b01; a11 = b11; a12 = b21; a13 = b31;

    //
    // I4 x L16_8
    //

    BBE_DSELNX16I(b01, b00, a01, a00, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(b03, b02, a03, a02, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(b11, b10, a11, a10, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(b13, b12, a13, a12, BBE_DSELI_INTERLEAVE_2);

    //
    // I2 x L16_8 x I2
    //

    BBE_DSELNX16I(a01, a00, b02, b00, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(a11, a10, b03, b01, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(a21, a20, b12, b10, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(a31, a30, b13, b11, BBE_DSELI_INTERLEAVE_4);

    //
    // DFT4 x I16
    //

    BBE_MOVSAV(a20); BBE_MOVSCV(a21);
    BBE_MOVSBV(a30); BBE_MOVSDV(a31);

    b00 = BBE_FFTADD4SABNX16(a00, a10, 0, 0);
    b01 = BBE_FFTADD4SCDNX16(a01, a11, 0, 0);
    b30 = BBE_FFTADD4SABNX16(a00, a10, 1, 0);
    b31 = BBE_FFTADD4SCDNX16(a01, a11, 1, 0);
    b20 = BBE_FFTADD4SABNX16(a00, a10, 2, 0);
    b21 = BBE_FFTADD4SCDNX16(a01, a11, 2, 0);
    b10 = BBE_FFTADD4SABNX16(a00, a10, 3, 0);
    b11 = BBE_FFTADD4SCDNX16(a01, a11, 3, 0);

    //
    // T16_4 x I4
    //

    b10 = BBE_MULNX16JPACKQ(b10, T16_4_10);
    b11 = BBE_MULNX16JPACKQ(b11, T16_4_11);
    b20 = BBE_MULNX16JPACKQ(b20, T16_4_20);
    b21 = BBE_MULNX16JPACKQ(b21, T16_4_21);
    b30 = BBE_MULNX16JPACKQ(b30, T16_4_30);
    b31 = BBE_MULNX16JPACKQ(b31, T16_4_31);

    //**************************************************************************
    // Stage 3: ( DFT4 x I16 )*( L16_4 x I4 )

    //
    // L16_4 x I4
    //

    a00 = BBE_SELNX16I(b10, b00, BBE_SELI_EXTRACT_LO_HALVES);
    a01 = BBE_SELNX16I(b30, b20, BBE_SELI_EXTRACT_LO_HALVES);
    a10 = BBE_SELNX16I(b10, b00, BBE_SELI_EXTRACT_HI_HALVES);
    a11 = BBE_SELNX16I(b30, b20, BBE_SELI_EXTRACT_HI_HALVES);
    a20 = BBE_SELNX16I(b11, b01, BBE_SELI_EXTRACT_LO_HALVES);
    a21 = BBE_SELNX16I(b31, b21, BBE_SELI_EXTRACT_LO_HALVES);
    a30 = BBE_SELNX16I(b11, b01, BBE_SELI_EXTRACT_HI_HALVES);
    a31 = BBE_SELNX16I(b31, b21, BBE_SELI_EXTRACT_HI_HALVES);

    //
    // DFT4 x I16
    //

    BBE_MOVSAV(a20); BBE_MOVSCV(a21);
    BBE_MOVSBV(a30); BBE_MOVSDV(a31);

    b00 = BBE_FFTADD4SABNX16(a00, a10, 0, 0);
    b01 = BBE_FFTADD4SCDNX16(a01, a11, 0, 0);
    b30 = BBE_FFTADD4SABNX16(a00, a10, 1, 0);
    b31 = BBE_FFTADD4SCDNX16(a01, a11, 1, 0);
    b20 = BBE_FFTADD4SABNX16(a00, a10, 2, 0);
    b21 = BBE_FFTADD4SCDNX16(a01, a11, 2, 0);
    b10 = BBE_FFTADD4SABNX16(a00, a10, 3, 0);
    b11 = BBE_FFTADD4SCDNX16(a01, a11, 3, 0);

    BBE_SVNX16_IP(b00, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b01, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b10, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b11, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b20, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b21, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b30, Y, 4 * BBE_SIMD_WIDTH / 2);
    BBE_SVNX16_IP(b31, Y, 4 * BBE_SIMD_WIDTH / 2);
  }
} /* bcifft64() */
#endif

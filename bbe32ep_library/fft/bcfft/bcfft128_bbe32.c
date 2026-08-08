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
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Twiddle tables declarations. */
#include "fft_tw.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(void, bcfft128, (complex_fract16 * restrict y, complex_fract16 * restrict x, int L))
#else

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

void bcfft128 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  xb_vecNx16 * restrict Y;
  xb_vecNx16 * restrict Y0;
  xb_vecNx16 * restrict Y1;
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          X0;
  const xb_vecNx16 *          X1;
  const xb_vecNx16 *          TW;

  valign va0, va1;

  int l;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  BBE_FFTWMODE(0x10);

  //----------------------------------------------------------------------------
  // Apply each stage of DFT128 to L data blocks. DFT128 decomposition:
  //   DFT128 -> ( L8_2 x I16 )*( I4 x DFT2 x I16 )*( T8_2 x I16 )*
  //             ( L8_4 x I16 )*( I2 x DFT4 x I16 )*( L8_2 x I16 )*
  //             ( I4 x L8_2 x I4 )*( L16_4 x I8 )*( T32_8 x I4 )*
  //             ( L16_4 x I8 )*( I4 x DFT4 x I8 )*( L16_4 x I8 )*
  //             ( I4 x L32_16 )*( I4 x L32_16 )*( L16_4 x I8 )*T128_32*
  //             ( L16_4 x I8 )*( I4 x DFT4 x I8 )*( L16_4 x I8 )
  //

  //----------------------------------------------------------------------------
  // Stage 1: ( I4 x L32_16 )*( I4 x L32_16 )*( L16_4 x I8 )*T128_32*
  //          ( L16_4 x I8 )*( I4 x DFT4 x I8 )*( L16_4 x I8 )

  {
    xb_vecNx16 T128_32_0, T128_32_1, T128_32_2, T128_32_3;
    xb_vecNx16 T128_32_4, T128_32_5, T128_32_6, T128_32_7;
    xb_vecNx16 T128_32_8, T128_32_9, T128_32_10, T128_32_11;

    Y = (xb_vecNx16*)y;
    X = (const xb_vecNx16*)x;
    TW = (const xb_vecNx16*)fft128_tw1;

    {
      xb_vecNx16 a00, a01, a02, a03, a10, a11, a12, a13;
      xb_vecNx16 a20, a21, a22, a23, a30, a31, a32, a33;
      xb_vecNx16 b00, b01, b02, b03, b10, b11, b12, b13;
      xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

      T128_32_0 = BBE_LVNX16_I(TW, 0 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_1 = BBE_LVNX16_I(TW, 1 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_2 = BBE_LVNX16_I(TW, 2 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_3 = BBE_LVNX16_I(TW, 3 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_4 = BBE_LVNX16_I(TW, 4 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_5 = BBE_LVNX16_I(TW, 5 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_6 = BBE_LVNX16_I(TW, 6 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_7 = BBE_LVNX16_I(TW, 7 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_8 = BBE_LVNX16_I(TW, 8 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_9 = BBE_LVNX16_I(TW, 9 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_10 = BBE_LVNX16_I(TW, 10 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_11 = BBE_LVNX16_I(TW, 11 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // (I4 x DFT4 x I8)*(L16_4 x I8)
      //

      BBE_LVNX16_IP(a00, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a02
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a12
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a03
      BBE_LVD_XP(X, -11 * 4 * BBE_SIMD_WIDTH / 2); // a13

      b00 = BBE_FFTADD4SABNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SABNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SABNX16(a00, a01, 3, 0);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      BBE_LVNX16_IP(a20, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a30, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a22
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a32
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a23
      BBE_LVD_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a33

      b20 = BBE_FFTADD4SABNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SABNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SABNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SABNX16(a20, a21, 3, 0);

      b30 = BBE_FFTADD4SCDNX16(a30, a31, 0, 0);
      b31 = BBE_FFTADD4SCDNX16(a30, a31, 1, 0);
      b32 = BBE_FFTADD4SCDNX16(a30, a31, 2, 0);
      b33 = BBE_FFTADD4SCDNX16(a30, a31, 3, 0);

      //
      // ( L16_4 x I8 )*T128_32*( L16_4 x I8 )
      //

      b01 = BBE_MULNX16CPACKQ(b01, T128_32_0);
      b02 = BBE_MULNX16CPACKQ(b02, T128_32_1);
      b03 = BBE_MULNX16CPACKQ(b03, T128_32_2);

      b11 = BBE_MULNX16CPACKQ(b11, T128_32_3);
      b12 = BBE_MULNX16CPACKQ(b12, T128_32_4);
      b13 = BBE_MULNX16CPACKQ(b13, T128_32_5);

      b21 = BBE_MULNX16CPACKQ(b21, T128_32_6);
      b22 = BBE_MULNX16CPACKQ(b22, T128_32_7);
      b23 = BBE_MULNX16CPACKQ(b23, T128_32_8);

      b31 = BBE_MULNX16CPACKQ(b31, T128_32_9);
      b32 = BBE_MULNX16CPACKQ(b32, T128_32_10);
      b33 = BBE_MULNX16CPACKQ(b33, T128_32_11);

      //
      // ( I4 x L32_16 )
      //

      BBE_DSELNX16I(a01, a00, b02, b00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a03, a02, b03, b01, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a11, a10, b12, b10, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a13, a12, b13, b11, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a21, a20, b22, b20, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a23, a22, b23, b21, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a31, a30, b32, b30, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a33, a32, b33, b31, BBE_DSELI_INTERLEAVE_2);

      //
      // ( I4 x L32_16 )
      //

      va0 = BBE_MOVUVR(a00);
      va1 = BBE_MOVUVR(a01);

      BBE_SVINTLARNX16_XP(a02, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a03, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(a10, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a11, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a12, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a13, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(a20, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a21, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a22, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a23, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(a30, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a31, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a32, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a33, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    }

    for (l = 1; l<L; l++)
    {
      xb_vecNx16 a00, a01, a02, a03, a10, a11, a12, a13;
      xb_vecNx16 a20, a21, a22, a23, a30, a31, a32, a33;
      xb_vecNx16 b00, b01, b02, b03, b10, b11, b12, b13;
      xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

      T128_32_0 = BBE_LVNX16_I(TW, 0 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_1 = BBE_LVNX16_I(TW, 1 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_2 = BBE_LVNX16_I(TW, 2 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_3 = BBE_LVNX16_I(TW, 3 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_4 = BBE_LVNX16_I(TW, 4 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_5 = BBE_LVNX16_I(TW, 5 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_6 = BBE_LVNX16_I(TW, 6 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_7 = BBE_LVNX16_I(TW, 7 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_8 = BBE_LVNX16_I(TW, 8 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_9 = BBE_LVNX16_I(TW, 9 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_10 = BBE_LVNX16_I(TW, 10 * 4 * BBE_SIMD_WIDTH / 2);
      T128_32_11 = BBE_LVNX16_I(TW, 11 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // (I4 x DFT4 x I8)*(L16_4 x I8)
      //

      BBE_LVNX16_IP(a00, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a02
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a12
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a03
      BBE_LVD_XP(X, -11 * 4 * BBE_SIMD_WIDTH / 2); // a13

      b00 = BBE_FFTADD4SABNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SABNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SABNX16(a00, a01, 3, 0);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      BBE_LVNX16_IP(a20, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a30, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a22
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a32
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a23
      BBE_LVD_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a33

      b20 = BBE_FFTADD4SABNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SABNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SABNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SABNX16(a20, a21, 3, 0);

      b30 = BBE_FFTADD4SCDNX16(a30, a31, 0, 0);
      b31 = BBE_FFTADD4SCDNX16(a30, a31, 1, 0);
      b32 = BBE_FFTADD4SCDNX16(a30, a31, 2, 0);
      b33 = BBE_FFTADD4SCDNX16(a30, a31, 3, 0);

      //
      // ( L16_4 x I8 )*T128_32*( L16_4 x I8 )
      //

      b01 = BBE_MULNX16CPACKQ(b01, T128_32_0);
      b02 = BBE_MULNX16CPACKQ(b02, T128_32_1);
      b03 = BBE_MULNX16CPACKQ(b03, T128_32_2);

      b11 = BBE_MULNX16CPACKQ(b11, T128_32_3);
      b12 = BBE_MULNX16CPACKQ(b12, T128_32_4);
      b13 = BBE_MULNX16CPACKQ(b13, T128_32_5);

      b21 = BBE_MULNX16CPACKQ(b21, T128_32_6);
      b22 = BBE_MULNX16CPACKQ(b22, T128_32_7);
      b23 = BBE_MULNX16CPACKQ(b23, T128_32_8);

      b31 = BBE_MULNX16CPACKQ(b31, T128_32_9);
      b32 = BBE_MULNX16CPACKQ(b32, T128_32_10);
      b33 = BBE_MULNX16CPACKQ(b33, T128_32_11);

      //
      // ( I4 x L32_16 )
      //

      BBE_DSELNX16I(a01, a00, b02, b00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a03, a02, b03, b01, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a11, a10, b12, b10, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a13, a12, b13, b11, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a21, a20, b22, b20, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a23, a22, b23, b21, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a31, a30, b32, b30, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a33, a32, b33, b31, BBE_DSELI_INTERLEAVE_2);

      //
      // ( I4 x L32_16 )
      //

      BBE_SALIGNVRNX16_XP(a00, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a01, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a02, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a03, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);

      BBE_SALIGNVRNX16_XP(a10, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a11, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a12, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a13, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);

      BBE_SALIGNVRNX16_XP(a20, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a21, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a22, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a23, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);

      BBE_SALIGNVRNX16_XP(a30, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a31, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a32, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(a33, va1, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    }

    {
      xb_vecNx16 z = 0;

      BBE_SALIGNVRNX16_XP(z, va0, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(z, va1, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 2: ( I4 x L8_2 x I4 )*( L16_4 x I8 )*( T32_8 x I4 )*
  //          ( L16_4 x I8 )*( I4 x DFT4 x I8 )*( L16_4 x I8 )

  {
    xb_vecNx16 T32_8_0, T32_8_1, T32_8_2;

    xb_vecNx16 tw01, tw02, tw03;
    xb_vecNx16 tw11, tw12, tw13;
    xb_vecNx16 tw21, tw22, tw23;
    xb_vecNx16 tw31, tw32, tw33;

    Y = (xb_vecNx16*)x;
    X = (const xb_vecNx16*)y;
    TW = (const xb_vecNx16*)fft128_tw2;

    {
      xb_vecNx16 a00, a01, a10, a11, a20, a21, a30, a31;
      xb_vecNx16 b00, b01, b02, b03, b10, b11, b12, b13;
      xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

      T32_8_0 = BBE_LVNX16_I(TW, 0 * 4 * BBE_SIMD_WIDTH / 2);
      T32_8_1 = BBE_LVNX16_I(TW, 1 * 4 * BBE_SIMD_WIDTH / 2);
      T32_8_2 = BBE_LVNX16_I(TW, 2 * 4 * BBE_SIMD_WIDTH / 2);

      tw01 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw02 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw03 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw11 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw12 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw13 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw21 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw22 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw23 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw31 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw32 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw33 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_3);

      //
      // (I4 x DFT4 x I8)*(L16_4 x I8)
      //

      BBE_LVNX16_IP(a00, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a02
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a12
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a03
      BBE_LVD_XP(X, -11 * 4 * BBE_SIMD_WIDTH / 2); // a13

      b00 = BBE_FFTADD4SABNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SABNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SABNX16(a00, a01, 3, 0);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      BBE_LVNX16_IP(a20, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a30, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a22
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a32
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a23
      BBE_LVD_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a33

      b20 = BBE_FFTADD4SABNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SABNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SABNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SABNX16(a20, a21, 3, 0);

      b30 = BBE_FFTADD4SCDNX16(a30, a31, 0, 0);
      b31 = BBE_FFTADD4SCDNX16(a30, a31, 1, 0);
      b32 = BBE_FFTADD4SCDNX16(a30, a31, 2, 0);
      b33 = BBE_FFTADD4SCDNX16(a30, a31, 3, 0);

      //
      // ( L16_4 x I8 )*( T32_8 x I4 )*( L16_4 x I8 )
      //

      b01 = BBE_MULNX16CPACKQ(b01, tw01);
      b02 = BBE_MULNX16CPACKQ(b02, tw02);
      b03 = BBE_MULNX16CPACKQ(b03, tw03);

      b11 = BBE_MULNX16CPACKQ(b11, tw11);
      b12 = BBE_MULNX16CPACKQ(b12, tw12);
      b13 = BBE_MULNX16CPACKQ(b13, tw13);

      b21 = BBE_MULNX16CPACKQ(b21, tw21);
      b22 = BBE_MULNX16CPACKQ(b22, tw22);
      b23 = BBE_MULNX16CPACKQ(b23, tw23);

      b31 = BBE_MULNX16CPACKQ(b31, tw31);
      b32 = BBE_MULNX16CPACKQ(b32, tw32);
      b33 = BBE_MULNX16CPACKQ(b33, tw33);

      //
      // I4 x L8_2 x I4
      //

      va0 = BBE_MOVUVR(b00);
      va1 = BBE_MOVUVR(b02);

      BBE_SVINTLARNX16_XP(b01, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b03, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(b10, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b12, va1, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b11, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b13, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(b20, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b22, va1, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b21, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b23, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(b30, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b32, va1, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b31, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b33, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
    }

    __Pragma("ymemory( TW )");
    for (l = 1; l<L; l++)
    {
      xb_vecNx16 a00, a01, a10, a11, a20, a21, a30, a31;
      xb_vecNx16 b00, b01, b02, b03, b10, b11, b12, b13;
      xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

      T32_8_0 = BBE_LVNX16_I(TW, 0 * 4 * BBE_SIMD_WIDTH / 2);
      T32_8_1 = BBE_LVNX16_I(TW, 1 * 4 * BBE_SIMD_WIDTH / 2);
      T32_8_2 = BBE_LVNX16_I(TW, 2 * 4 * BBE_SIMD_WIDTH / 2);

      tw01 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw02 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw03 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw11 = BBE_SHFLNX16I(T32_8_0, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw12 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw13 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw21 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw22 = BBE_SHFLNX16I(T32_8_1, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw23 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw31 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw32 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw33 = BBE_SHFLNX16I(T32_8_2, BBE_SHFLI_REP_2X4_OFFSET_3);

      //
      // (I4 x DFT4 x I8)*(L16_4 x I8)
      //

      BBE_LVNX16_IP(a00, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a02
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a12
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a03
      BBE_LVD_XP(X, -11 * 4 * BBE_SIMD_WIDTH / 2); // a13

      b00 = BBE_FFTADD4SABNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SABNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SABNX16(a00, a01, 3, 0);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      BBE_LVNX16_IP(a20, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a30, X, +3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X, +3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a22
      BBE_LVC_IP(X, +3 * 4 * BBE_SIMD_WIDTH / 2); // a32
      BBE_LVB_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a23
      BBE_LVD_IP(X, +1 * 4 * BBE_SIMD_WIDTH / 2); // a33

      b20 = BBE_FFTADD4SABNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SABNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SABNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SABNX16(a20, a21, 3, 0);

      b30 = BBE_FFTADD4SCDNX16(a30, a31, 0, 0);
      b31 = BBE_FFTADD4SCDNX16(a30, a31, 1, 0);
      b32 = BBE_FFTADD4SCDNX16(a30, a31, 2, 0);
      b33 = BBE_FFTADD4SCDNX16(a30, a31, 3, 0);

      //
      // ( L16_4 x I8 )*( T32_8 x I4 )*( L16_4 x I2 )
      //

      b01 = BBE_MULNX16CPACKQ(b01, tw01);
      b02 = BBE_MULNX16CPACKQ(b02, tw02);
      b03 = BBE_MULNX16CPACKQ(b03, tw03);

      b11 = BBE_MULNX16CPACKQ(b11, tw11);
      b12 = BBE_MULNX16CPACKQ(b12, tw12);
      b13 = BBE_MULNX16CPACKQ(b13, tw13);

      b21 = BBE_MULNX16CPACKQ(b21, tw21);
      b22 = BBE_MULNX16CPACKQ(b22, tw22);
      b23 = BBE_MULNX16CPACKQ(b23, tw23);

      b31 = BBE_MULNX16CPACKQ(b31, tw31);
      b32 = BBE_MULNX16CPACKQ(b32, tw32);
      b33 = BBE_MULNX16CPACKQ(b33, tw33);

      //
      // I4 x L8_2 x I4
      //

      BBE_SALIGNVRNX16_XP(b00, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b02, va1, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(b01, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b03, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);

      BBE_SALIGNVRNX16_XP(b10, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b12, va1, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(b11, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b13, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);

      BBE_SALIGNVRNX16_XP(b20, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b22, va1, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(b21, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b23, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);

      BBE_SALIGNVRNX16_XP(b30, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b32, va1, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(b31, va0, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(b33, va1, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
    }

    {
      xb_vecNx16 z = 0;

      BBE_SALIGNVRNX16_XP(z, va0, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(z, va1, Y, 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 3: ( L8_2 x I16 )*( I4 x DFT2 x I16 )*( T8_2 x I16 )*
  //          ( L8_4 x I16 )*( I2 x DFT4 x I16 )*( L8_2 x I16 )

  {
    xb_vecNx16 T8_2_1, T8_2_2;

    Y0 = (xb_vecNx16*)((uintptr_t)y + 0 * 4 * BBE_SIMD_WIDTH / 2);
    Y1 = (xb_vecNx16*)((uintptr_t)y + 1 * 4 * BBE_SIMD_WIDTH / 2);
    X0 = (const xb_vecNx16*)((uintptr_t)x + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)x + 2 * 4 * BBE_SIMD_WIDTH / 2);

    T8_2_1 = BBE_MOVVA16C(0xa57e5a82);
    T8_2_2 = BBE_MOVVINX16(BBE_MOVVI_CQ15_MI);

    for (l = 0; l<L; l++)
    {
      xb_vecNx16 a00, a01, a02, a03, a10, a11, a12, a13;
      xb_vecNx16 a20, a21, a22, a23, a30, a31, a32, a33;
      xb_vecNx16 b00, b01, b02, b03, b04, b05, b06, b07;
      xb_vecNx16 b10, b11, b12, b13, b14, b15, b16, b17;
      xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

      //
      // ( I2 x DFT4 x I16 )*( L8_2 x I16 )
      //

      BBE_LVNX16_IP(a00, X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X0, 3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a02, X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a03, X0, 3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVA_IP(X0, +4 * 4 * BBE_SIMD_WIDTH / 2); // a04
      BBE_LVB_XP(X0, -3 * 4 * BBE_SIMD_WIDTH / 2); // a06

      b00 = BBE_FFTADD4SABNX16(a00, a02, 0, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a02, 1, 0);
      b04 = BBE_FFTADD4SABNX16(a00, a02, 2, 0);
      b06 = BBE_FFTADD4SABNX16(a00, a02, 3, 0);

      BBE_LVA_IP(X0, +4 * 4 * BBE_SIMD_WIDTH / 2); // a05
      BBE_LVB_IP(X0, +3 * 4 * BBE_SIMD_WIDTH / 2); // a07

      b01 = BBE_FFTADD4SABNX16(a01, a03, 0, 0);
      b03 = BBE_FFTADD4SABNX16(a01, a03, 1, 0);
      b05 = BBE_FFTADD4SABNX16(a01, a03, 2, 0);
      b07 = BBE_FFTADD4SABNX16(a01, a03, 3, 0);

      BBE_LVNX16_IP(a10, X1, 1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X1, 3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a12, X1, 1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a13, X1, 3 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVC_IP(X1, +4 * 4 * BBE_SIMD_WIDTH / 2); // a14
      BBE_LVD_XP(X1, -3 * 4 * BBE_SIMD_WIDTH / 2); // a16

      b10 = BBE_FFTADD4SCDNX16(a10, a12, 0, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a12, 1, 0);
      b14 = BBE_FFTADD4SCDNX16(a10, a12, 2, 0);
      b16 = BBE_FFTADD4SCDNX16(a10, a12, 3, 0);

      BBE_LVC_IP(X1, +4 * 4 * BBE_SIMD_WIDTH / 2); // a15
      BBE_LVD_IP(X1, +3 * 4 * BBE_SIMD_WIDTH / 2); // a17

      b11 = BBE_FFTADD4SCDNX16(a11, a13, 0, 0);
      b13 = BBE_FFTADD4SCDNX16(a11, a13, 1, 0);
      b15 = BBE_FFTADD4SCDNX16(a11, a13, 2, 0);
      b17 = BBE_FFTADD4SCDNX16(a11, a13, 3, 0);

      //
      // ( T8_2 x I16 )*( L8_4 x I16 )
      //

      a00 = b00; a01 = b01; a02 = b10; a03 = b11;
      a10 = b02; a11 = b03; a12 = b12; a13 = b13;
      a20 = b04; a21 = b05; a22 = b14; a23 = b15;
      a30 = b06; a31 = b07; a32 = b16; a33 = b17;

      a12 = BBE_MULNX16CPACKQ(a12, T8_2_1);
      a13 = BBE_MULNX16CPACKQ(a13, T8_2_1);
      a22 = BBE_MULNX16CPACKQ(a22, T8_2_2);
      a23 = BBE_MULNX16CPACKQ(a23, T8_2_2);
      // a32 and a33 are implicitly negated to reuse the T8_2_1 twiddle factor.
      a32 = BBE_MULNX16JPACKQ(a32, T8_2_1);
      a33 = BBE_MULNX16JPACKQ(a33, T8_2_1);

      //
      // I4 x DFT2 x I16
      //

      b00 = BBE_FFTADDSSRNX16(a00, a02);
      b01 = BBE_FFTADDSSRNX16(a01, a03);
      b02 = BBE_FFTSUBSSRNX16(a00, a02);
      b03 = BBE_FFTSUBSSRNX16(a01, a03);

      b10 = BBE_FFTADDSSRNX16(a10, a12);
      b11 = BBE_FFTADDSSRNX16(a11, a13);
      b12 = BBE_FFTSUBSSRNX16(a10, a12);
      b13 = BBE_FFTSUBSSRNX16(a11, a13);

      b20 = BBE_FFTADDSSRNX16(a20, a22);
      b21 = BBE_FFTADDSSRNX16(a21, a23);
      b22 = BBE_FFTSUBSSRNX16(a20, a22);
      b23 = BBE_FFTSUBSSRNX16(a21, a23);

      b30 = BBE_FFTSUBSSRNX16(a30, a32);
      b31 = BBE_FFTSUBSSRNX16(a31, a33);
      b32 = BBE_FFTADDSSRNX16(a30, a32);
      b33 = BBE_FFTADDSSRNX16(a31, a33);

      //
      // L8_2 x I16
      //

      BBE_SVNX16_IP(b00, Y0, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b01, Y1, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b02, Y0, -6 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b03, Y1, -6 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(b10, Y0, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b11, Y1, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b12, Y0, -6 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b13, Y1, -6 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(b20, Y0, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b21, Y1, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b22, Y0, -6 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_XP(b23, Y1, -6 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(b30, Y0, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b31, Y1, +8 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b32, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(b33, Y1, +2 * 4 * BBE_SIMD_WIDTH / 2);
    }
  }
} /* bcfft128() */
#endif

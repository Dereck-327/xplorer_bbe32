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
    Blockwise mixed radix forward FFT on complex data, no data scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Twiddle factor tables. */
#include "fft_tw.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(void, bcnfft108,( void * restrict pScr,
                complex_fract16 * restrict y,
                complex_fract16 * restrict x,
                int L ))
#else
/*-------------------------------------------------------------------------
Blockwise mixed radix forward FFT on complex data, no data scaling
  
Description: These functions make forward FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of BCNFFT_SCRATCH_SIZE(N) bytes
  Input:
    S           Required input/output buffer size may exceed actual data size:
                S >= N. Use BCNFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[L][S]     Complex input signal
  Output:            
    y[L][S]     Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
  L>0           The number of blocks must be positive
-------------------------------------------------------------------------*/

void bcnfft108 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  //
  // DFT108 decomposition:
  //
  //   DFT108 = ( I9 x C16_12 )*( L18_3 x I8 )* ...
  //            ( I6 x DFT3 x I8 )*( L18_6 x I8 )* ...
  //            ...
  //            ( L9_3 x I16 )*( L18_3 x I8 )*( Tp9_3 x I8 )* ...
  //            ( I6 x DFT3 x I8 )*( L18_6 x I8 )* ...
  //            ...
  //            ( L36_9 x I4 )*( E3_4 x I36 )*( I3 x C10_9 x I4 )* ...
  //            ( L15_3 x I8 )*( Tep27_9 x I4 )*( I5 x DFT3 x I8 )* ...
  //            ( L15_5 x I8 )*( I3 x E9_10 x I4 )* ...
  //            ...
  //            ( C32_27 x I4 )*L128_32*( L16_4 x I8 )*Tep108_27* ...
  //            ( I4 x DFT4 x I8 )*( L16_4 x I8 )*( I4 x E27_32 )
  //
  //   Tp9_3 = L18_6*( T9_3 x I2 )*L18_3
  //
  //   Tep27_9 = P*T27_9*P', where P = ( L15_5 x I2 )*( I3 x E9_10 )
  //
  //   Tep108_27 = P*T108_27*P', where P = ( L16_4 x I8 )*( I4 x E27_32 )
  //

  const xb_vecNx16 *          X;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 * restrict Y0;
  xb_vecNx16 * restrict Y1;
  xb_vecNx16 * restrict Y2;
  const xb_vecNx16 *          T;

  int l;

  NASSERT_ALIGN32(pScr);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Stage 1: ( C32_27 x I4 )*L128_32*( L16_4 x I8 )*Tep108_27*
  //          ( I4 x DFT4 x I8 )*( L16_4 x I8 )*( I4 x E27_32 )

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)pScr;
  T = (const xb_vecNx16*)fft108_tw1;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;
    const xb_vecNx16 * X3;

    valign X1_va, X2_va, Y0_va, Y1_va;

    xb_vecNx16 a00, a01, a02, a03;
    xb_vecNx16 a10, a11, a12, a13;
    xb_vecNx16 a20, a21, a22, a23;
    xb_vecNx16 a30, a31, a32, a33;
    xb_vecNx16 b00, b01, b02, b03;
    xb_vecNx16 b10, b11, b12, b13;
    xb_vecNx16 b20, b21, b22, b23;
    xb_vecNx16 b30, b31, b32, b33;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5, tw6, tw7, tw8, tw9, tw10, tw11, tw12;

    BBE_FFTWMODE(0x10);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 27 * 4);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 54 * 4);
    X3 = (const xb_vecNx16*)((uintptr_t)X + 80 * 4);

    for (l = 0; l<L; l++)
    {
      //
      // Load twiddle table Tep108_27
      //

      tw1 = BBE_LVNX16_I(T, 0 * 4 * BBE_SIMD_WIDTH / 2);
      tw2 = BBE_LVNX16_I(T, 1 * 4 * BBE_SIMD_WIDTH / 2);
      tw3 = BBE_LVNX16_I(T, 2 * 4 * BBE_SIMD_WIDTH / 2);
      tw4 = BBE_LVNX16_I(T, 3 * 4 * BBE_SIMD_WIDTH / 2);
      tw5 = BBE_LVNX16_I(T, 4 * 4 * BBE_SIMD_WIDTH / 2);
      tw6 = BBE_LVNX16_I(T, 5 * 4 * BBE_SIMD_WIDTH / 2);
      tw7 = BBE_LVNX16_I(T, 6 * 4 * BBE_SIMD_WIDTH / 2);
      tw8 = BBE_LVNX16_I(T, 7 * 4 * BBE_SIMD_WIDTH / 2);
      tw9 = BBE_LVNX16_I(T, 8 * 4 * BBE_SIMD_WIDTH / 2);
      tw10 = BBE_LVNX16_I(T, 9 * 4 * BBE_SIMD_WIDTH / 2);
      tw11 = BBE_LVNX16_I(T, 10 * 4 * BBE_SIMD_WIDTH / 2);
      tw12 = BBE_LVNX16_I(T, 11 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // ( L16_4 x I8 )*( I4 x E27_32 )
      //

      BBE_LVNX16_IP(a00, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a20, X0, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_XP(a30, X0, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 3 * 4 * BBE_SIMD_WIDTH / 2);

      X1_va = BBE_LANX16_PP(X1);

      BBE_LANX16_IP(a01, X1_va, X1);
      BBE_LANX16_IP(a11, X1_va, X1);
      BBE_LANX16_IP(a21, X1_va, X1);
      BBE_LANX16_IP(a31, X1_va, X1);

      X1 = (const xb_vecNx16*)XT_ADD((uintptr_t)X1, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);

      X2_va = BBE_LANX16_PP(X2);

      BBE_LANX16_IP(a02, X2_va, X2);
      BBE_LANX16_IP(a12, X2_va, X2);
      BBE_LANX16_IP(a22, X2_va, X2);
      BBE_LANX16_IP(a32, X2_va, X2);

      X2 = (const xb_vecNx16*)XT_ADD((uintptr_t)X2, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a03, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a13, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a23, X3, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_XP(a33, X3, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 3 * 4 * BBE_SIMD_WIDTH / 2);

      a03 = BBE_SELNX16I(a13, a03, BBE_SELI_ROTATE_RIGHT_2);
      a13 = BBE_SELNX16I(a23, a13, BBE_SELI_ROTATE_RIGHT_2);
      a23 = BBE_SELNX16I(a33, a23, BBE_SELI_ROTATE_RIGHT_2);
      a33 = BBE_SELNX16I(a33, a33, BBE_SELI_ROTATE_RIGHT_2);

      //
      // Tep108_27*( I4 x DFT4 x I8 )
      //

      BBE_MOVSCV(a02);
      BBE_MOVSDV(a03);

      b00 = BBE_FFTADD4SCDNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SCDNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SCDNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SCDNX16(a00, a01, 3, 0);

      b01 = BBE_MULNX16CPACKQ(b01, tw1);
      b02 = BBE_MULNX16CPACKQ(b02, tw2);
      b03 = BBE_MULNX16CPACKQ(b03, tw3);

      BBE_MOVSCV(a12);
      BBE_MOVSDV(a13);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      b11 = BBE_MULNX16CPACKQ(b11, tw4);
      b12 = BBE_MULNX16CPACKQ(b12, tw5);
      b13 = BBE_MULNX16CPACKQ(b13, tw6);

      BBE_MOVSAV(a22);
      BBE_MOVSBV(a23);

      b20 = BBE_FFTADD4SABNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SABNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SABNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SABNX16(a20, a21, 3, 0);

      b21 = BBE_MULNX16CPACKQ(b21, tw7);
      b22 = BBE_MULNX16CPACKQ(b22, tw8);
      b23 = BBE_MULNX16CPACKQ(b23, tw9);

      BBE_MOVSAV(a32);
      BBE_MOVSBV(a33);

      b30 = BBE_FFTADD4SABNX16(a30, a31, 0, 0);
      b31 = BBE_FFTADD4SABNX16(a30, a31, 1, 0);
      b32 = BBE_FFTADD4SABNX16(a30, a31, 2, 0);
      b33 = BBE_FFTADD4SABNX16(a30, a31, 3, 0);

      b31 = BBE_MULNX16CPACKQ(b31, tw10);
      b32 = BBE_MULNX16CPACKQ(b32, tw11);
      b33 = BBE_MULNX16CPACKQ(b33, tw12);

      //
      // ( C32_27 x I4 )*L128_32*( L16_4 x I8 )
      //
      // L128_32 == L128_64^2
      //

      a00 = b00; a01 = b10; a02 = b20; a03 = b30;
      a10 = b01; a11 = b11; a12 = b21; a13 = b31;
      a20 = b02; a21 = b12; a22 = b22; a23 = b32;
      a30 = b03; a31 = b13; a32 = b23; a33 = b33;

      BBE_DSELNX16I(b01, b00, a20, a00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b03, b02, a21, a01, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b11, b10, a22, a02, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b13, b12, a23, a03, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b21, b20, a30, a10, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b23, b22, a31, a11, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b31, b30, a32, a12, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b33, b32, a33, a13, BBE_DSELI_INTERLEAVE_2);

      b32 = BBE_SELNX16I(a33, a13, BBE_SELI_INTERLEAVE_2_LO);

      Y0_va = BBE_MOVUVR(b00);
      Y1_va = BBE_MOVUVR(b01);

      BBE_SVINTLARNX16_XP(b20, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(b21, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(b02, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b03, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b22, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(b23, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(b10, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b11, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b30, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(b31, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(b12, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b31, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b32, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 0);

      BBE_SALIGNVRNX16_XP(b32, Y0_va, Y, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 13 * 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 2: ( L36_9 x I4 )*( E3_4 x I36 )*( I3 x C10_9 x I4 )*
  //          ( L15_3 x I8 )*( Tep27_9 x I4 )*( I5 x DFT3 x I8 )*
  //          ( L15_5 x I8 )*( I3 x E9_10 x I4 )

  X = (const xb_vecNx16*)pScr;
  Y = (xb_vecNx16*)x;
  T = (const xb_vecNx16*)fft108_tw2;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;
    const xb_vecNx16 * X3;
    const xb_vecNx16 * X4;

    valign Y0_va, Y1_va;

    xb_vecNx16 a00, a01, a02, a03, a04;
    xb_vecNx16 a10, a11, a12, a13, a14;
    xb_vecNx16 a20, a21, a22, a23, a24;
    xb_vecNx16 a30, a31, a32, a41, a42;

    xb_vecNx16 b00, b01, b02;
    xb_vecNx16 b10, b11, b12;
    xb_vecNx16 b20, b21, b22;
    xb_vecNx16 b30, b31, b32;
    xb_vecNx16 b40, b41, b42;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5;
    xb_vecNx16 tw6, tw7, tw8, tw9, tw10;

    xb_vecNx16 r3tw = BBE_MOVVA16C(0x91260000);

    xb_vecNx16 t0, t1, t2;

    t0 = 0;

    BBE_MOVSDV(t0);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X3 = (const xb_vecNx16*)((uintptr_t)X + 32 * 4);
    X4 = (const xb_vecNx16*)((uintptr_t)X + 72 * 4);

    for (l = 0; l<L; l++)
    {
      //
      // Load twiddle factor table Tep27_9
      //

      t0 = BBE_LVNX16_I(T, 0 * 4 * BBE_SIMD_WIDTH / 2);
      t1 = BBE_LVNX16_I(T, 1 * 4 * BBE_SIMD_WIDTH / 2);
      t2 = BBE_LVNX16_I(T, 2 * 4 * BBE_SIMD_WIDTH / 2);

      tw1 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw2 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw3 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw4 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw5 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw6 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw7 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw8 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw9 = BBE_SHFLNX16I(t2, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw10 = BBE_SHFLNX16I(t2, BBE_SHFLI_REP_2X4_OFFSET_1);

      //
      // ( L15_5 x I8 )*( I3 x E9_10 x I4 )
      //

      a10 = BBE_LVNX16_I(X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
      a20 = BBE_LVNX16_I(X0, 2 * 4 * BBE_SIMD_WIDTH / 2);
      a30 = BBE_LVNX16_I(X0, 3 * 4 * BBE_SIMD_WIDTH / 2);

      X0 = (const xb_vecNx16*)XT_ADD((uintptr_t)X0, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16));

      BBE_LVNX16_IP(a01, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X3, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_XP(a41, X3, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);

      a01 = BBE_SELNX16I(a11, a01, BBE_SELI_ROTATE_RIGHT_8);
      a11 = BBE_SELNX16I(a21, a11, BBE_SELI_ROTATE_RIGHT_8);
      a21 = BBE_SELNX16I(a31, a21, BBE_SELI_ROTATE_RIGHT_8);
      a31 = BBE_SELNX16I(a41, a31, BBE_SELI_ROTATE_RIGHT_8);

      a41 = BBE_SHFLNX16I(a41, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

      BBE_LVNX16_IP(a02, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a12, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a22, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a32, X4, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_XP(a42, X4, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // ( Tep27_9 x I4 )*( I5 x DFT3 x I8 )
      //
      // Radix-3 butterfly: b0 = a0 + a1 + a2
      //                    b1 = a1*r3tw - a2*r3tw + a0 - (a1+a2)/2
      //                    b2 = a2*r3tw - a1*r3tw + a0 - (a1+a2)/2
      //
      // where r3tw = sign*1j*(3^0.5)/2
      //

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 4 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a01, a02);

      t0 = BBE_MULNX16CPACKQ(a01, r3tw);
      t1 = BBE_MULNX16CPACKQ(a02, r3tw);

      // b0 <= a1 + a2 + C
      b00 = BBE_FFTADD4SCDNX16(a01, a02, 0, 0);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b01 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b02 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b01 = BBE_MULNX16CPACKQ(b01, tw1);
      b02 = BBE_MULNX16CPACKQ(b02, tw2);

      //--------------------

      // A <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a11, a12);

      t0 = BBE_MULNX16CPACKQ(a11, r3tw);
      t1 = BBE_MULNX16CPACKQ(a12, r3tw);

      // b0 <= a0 + a1 + a2
      b10 = BBE_ADDSNX16(a10, a11);
      b10 = BBE_ADDSNX16(b10, a12);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b11 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b12 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b11 = BBE_MULNX16CPACKQ(b11, tw3);
      b12 = BBE_MULNX16CPACKQ(b12, tw4);

      //--------------------

      // A <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a21, a22);

      t0 = BBE_MULNX16CPACKQ(a21, r3tw);
      t1 = BBE_MULNX16CPACKQ(a22, r3tw);

      // b0 <= a0 + a1 + a2
      b20 = BBE_ADDSNX16(a20, a21);
      b20 = BBE_ADDSNX16(b20, a22);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b21 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b22 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b21 = BBE_MULNX16CPACKQ(b21, tw5);
      b22 = BBE_MULNX16CPACKQ(b22, tw6);

      //--------------------

      // A <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a31, a32);

      t0 = BBE_MULNX16CPACKQ(a31, r3tw);
      t1 = BBE_MULNX16CPACKQ(a32, r3tw);

      // b0 <= a0 + a1 + a2
      b30 = BBE_ADDSNX16(a30, a31);
      b30 = BBE_ADDSNX16(b30, a32);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b31 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b32 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b31 = BBE_MULNX16CPACKQ(b31, tw7);
      b32 = BBE_MULNX16CPACKQ(b32, tw8);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 4 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a41, a42);

      t0 = BBE_MULNX16CPACKQ(a41, r3tw);
      t1 = BBE_MULNX16CPACKQ(a42, r3tw);

      // b0 <= a1 + a2 + C
      b40 = BBE_FFTADD4SCDNX16(a41, a42, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b41 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b42 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b41 = BBE_MULNX16CPACKQ(b41, tw9);
      b42 = BBE_MULNX16CPACKQ(b42, tw10);

      //
      // ( L36_9 x I4 )*( E3_4 x I36 )*( I3 x C10_9 x I4 )*( L15_3 x I8 )
      //

      a00 = b00; a01 = b10; a02 = b20; a03 = b30; a04 = b40;
      a10 = b01; a11 = b11; a12 = b21; a13 = b31; a14 = b41;
      a20 = b02; a21 = b12; a22 = b22; a23 = b32; a24 = b42;

      t0 = 0;

      Y0_va = BBE_MOVUVR(a00);
      Y1_va = BBE_MOVUVR(a20);

      BBE_SVINTLARNX16_XP(a10, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a10, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a01, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a21, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a11, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a11, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a02, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a22, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a12, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a12, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a03, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a23, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a13, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a13, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a04, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a24, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(a14, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);

      BBE_SVINTLARNX16_XP(a14, Y1_va, Y, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 17 * 4 * BBE_SIMD_WIDTH / 2, 1);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 3: ( L9_3 x I16 )*( L18_3 x I8 )*( Tp9_3 x I8 )*
  //          ( I6 x DFT3 x I8 )*( L18_6 x I8 )

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)pScr;
  T = (const xb_vecNx16*)fft108_tw4;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;
    const xb_vecNx16 * X3;
    const xb_vecNx16 * X4;

    xb_vecNx16 a00, a01, a02, a03, a04, a05;
    xb_vecNx16 a10, a11, a12, a13, a14, a15;
    xb_vecNx16 a20, a21, a22, a23, a24, a25;
    xb_vecNx16 a30, a31, a32;
    xb_vecNx16 a41, a42, a51, a52;
    xb_vecNx16 b00, b01, b02, b10, b11, b12;
    xb_vecNx16 b20, b21, b22, b30, b31, b32;
    xb_vecNx16 b40, b41, b42, b50, b51, b52;

    xb_vecNx16 t0, t1;

    xb_vecNx16 tw1, tw2, tw3, tw4;

    xb_vecNx16 r3tw = BBE_MOVVA16C(0x91260000);

    t0 = 0;

    BBE_MOVSDV(t0);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (16 * 108 / 12) / 3 + 3 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (16 * 108 / 12) / 3);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (16 * 108 / 12) / 3);
    X3 = (const xb_vecNx16*)((uintptr_t)X + 1 * 4 * (16 * 108 / 12) / 3);
    X4 = (const xb_vecNx16*)((uintptr_t)X + 2 * 4 * (16 * 108 / 12) / 3);

    for (l = 0; l<L; l++)
    {
      //
      // Load twiddle table Tp9_3
      //

      t0 = BBE_LVNX16_I(T, 0);

      tw1 = BBE_REPNX16C(t0, 2);
      tw2 = BBE_REPNX16C(t0, 3);
      tw3 = BBE_REPNX16C(t0, 4);
      tw4 = BBE_REPNX16C(t0, 5);

      //
      // L18_6 x I8
      //

      BBE_LVNX16_XP(a30, X0, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16));

      BBE_LVNX16_IP(a01, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a41, X3, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a51, X3, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a02, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a12, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a22, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a32, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a42, X4, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a52, X4, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // ( Tp9_3 x I8 )*( I6 x DFT3 x I8 )
      //

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a01, a02);

      t0 = BBE_MULNX16CPACKQ(a01, r3tw);
      t1 = BBE_MULNX16CPACKQ(a02, r3tw);

      // b0 <= a1 + a2 + C
      b00 = BBE_FFTADD4SCDNX16(a01, a02, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b01 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b02 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a11, a12);

      t0 = BBE_MULNX16CPACKQ(a11, r3tw);
      t1 = BBE_MULNX16CPACKQ(a12, r3tw);

      // b0 <= a1 + a2 + C
      b10 = BBE_FFTADD4SCDNX16(a11, a12, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b11 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b12 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 2 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a21, a22);

      t0 = BBE_MULNX16CPACKQ(a21, r3tw);
      t1 = BBE_MULNX16CPACKQ(a22, r3tw);

      // b0 <= a1 + a2 + C
      b20 = BBE_FFTADD4SCDNX16(a21, a22, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b21 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b22 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b21 = BBE_MULNX16CPACKQ(b21, tw1);
      b22 = BBE_MULNX16CPACKQ(b22, tw2);

      //--------------------

      // A <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a31, a32);

      t0 = BBE_MULNX16CPACKQ(a31, r3tw);
      t1 = BBE_MULNX16CPACKQ(a32, r3tw);

      // b0 <= a0 + a1 + a2
      b30 = BBE_FFTADDSSRNX16(a30, a31);
      b30 = BBE_FFTADDSSRNX16(b30, a32);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b31 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b32 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b31 = BBE_MULNX16CPACKQ(b31, tw1);
      b32 = BBE_MULNX16CPACKQ(b32, tw2);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a41, a42);

      t0 = BBE_MULNX16CPACKQ(a41, r3tw);
      t1 = BBE_MULNX16CPACKQ(a42, r3tw);

      // b0 <= a1 + a2 + C
      b40 = BBE_FFTADD4SCDNX16(a41, a42, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b41 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b42 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b41 = BBE_MULNX16CPACKQ(b41, tw3);
      b42 = BBE_MULNX16CPACKQ(b42, tw4);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a51, a52);

      t0 = BBE_MULNX16CPACKQ(a51, r3tw);
      t1 = BBE_MULNX16CPACKQ(a52, r3tw);

      // b0 <= a1 + a2 + C
      b50 = BBE_FFTADD4SCDNX16(a51, a52, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b51 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b52 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b51 = BBE_MULNX16CPACKQ(b51, tw3);
      b52 = BBE_MULNX16CPACKQ(b52, tw4);

      //
      // ( L9_3 x I16 )*( L18_3 x I8 )
      //

      a00 = b00; a01 = b10; a02 = b20; a03 = b30; a04 = b40; a05 = b50;
      a10 = b01; a11 = b11; a12 = b21; a13 = b31; a14 = b41; a15 = b51;
      a20 = b02; a21 = b12; a22 = b22; a23 = b32; a24 = b42; a25 = b52;

      BBE_SVNX16_IP(a00, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a01, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a10, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a11, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a20, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a21, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(a02, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a03, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a12, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a13, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a22, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a23, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(a04, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a05, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a14, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a15, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a24, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVNX16_IP(a25, Y, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 17 * 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 4: ( I9 x C16_12 )*( L18_3 x I8 )*( I6 x DFT3 x I8 )*( L18_6 x I8 )
  //

  X = (const xb_vecNx16*)pScr;
  Y = (xb_vecNx16*)y;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;
    const xb_vecNx16 * X3;
    const xb_vecNx16 * X4;

    valign Y0_va, Y1_va, Y2_va;

    xb_vecNx16 a00, a01, a02, a03, a04, a05;
    xb_vecNx16 a10, a11, a12, a13, a14, a15;
    xb_vecNx16 a20, a21, a22, a23, a24, a25;
    xb_vecNx16 a30, a31, a32;
    xb_vecNx16 a41, a42, a50, a51, a52;
    xb_vecNx16 b00, b01, b02, b10, b11, b12;
    xb_vecNx16 b20, b21, b22, b30, b31, b32;
    xb_vecNx16 b40, b41, b42, b50, b51, b52;

    xb_vecNx16 t0, t1;

    xb_vecNx16 r3tw = BBE_MOVVA16C(0x91260000);

    t0 = 0;

    BBE_MOVSDV(t0);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (108 * 16 / 12) / 3 + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (108 * 16 / 12) / 3 + 1 * 4 * BBE_SIMD_WIDTH / 2);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * (108 * 16 / 12) / 3 + 1 * 4 * BBE_SIMD_WIDTH / 2);
    X3 = (const xb_vecNx16*)((uintptr_t)X + 1 * 4 * (108 * 16 / 12) / 3);
    X4 = (const xb_vecNx16*)((uintptr_t)X + 2 * 4 * (108 * 16 / 12) / 3);

    Y0 = (xb_vecNx16*)((uintptr_t)Y + 0 * 4 * 108 / 3);
    Y1 = (xb_vecNx16*)((uintptr_t)Y + 1 * 4 * 108 / 3);
    Y2 = (xb_vecNx16*)((uintptr_t)Y + 2 * 4 * 108 / 3);

    for (l = 0; l<L; l++)
    {
      //
      // L18_6 x I8
      //

      BBE_LVNX16_IP(a00, X0, 3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a30, X0, 2 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_XP(a50, X0, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a01, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a21, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a31, X3, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a41, X3, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a51, X3, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a02, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a12, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a22, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a32, X4, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a42, X4, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a52, X4, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      //
      // ( Tp9_3 x I8 )*( I6 x DFT3 x I8 )
      //

      // A <= a0
      BBE_MOVSAV(a00);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a01, a02);

      t0 = BBE_MULNX16CPACKQ(a01, r3tw);
      t1 = BBE_MULNX16CPACKQ(a02, r3tw);

      // b0 <= a0 + a1 + a2
      b00 = BBE_FFTADDSSRNX16(a00, a01);
      b00 = BBE_FFTADDSSRNX16(b00, a02);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b01 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b02 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a11, a12);

      t0 = BBE_MULNX16CPACKQ(a11, r3tw);
      t1 = BBE_MULNX16CPACKQ(a12, r3tw);

      // b0 <= a1 + a2 + C
      b10 = BBE_FFTADD4SCDNX16(a11, a12, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b11 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b12 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, 2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, 2 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a21, a22);

      t0 = BBE_MULNX16CPACKQ(a21, r3tw);
      t1 = BBE_MULNX16CPACKQ(a22, r3tw);

      // b0 <= a1 + a2 + C
      b20 = BBE_FFTADD4SCDNX16(a21, a22, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b21 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b22 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A <= a0
      BBE_MOVSAV(a30);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a31, a32);

      t0 = BBE_MULNX16CPACKQ(a31, r3tw);
      t1 = BBE_MULNX16CPACKQ(a32, r3tw);

      // b0 <= a0 + a1 + a2
      b30 = BBE_FFTADDSSRNX16(a30, a31);
      b30 = BBE_FFTADDSSRNX16(b30, a32);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b31 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b32 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 3 * 4 * BBE_SIMD_WIDTH / 2);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a41, a42);

      t0 = BBE_MULNX16CPACKQ(a41, r3tw);
      t1 = BBE_MULNX16CPACKQ(a42, r3tw);

      // b0 <= a1 + a2 + C
      b40 = BBE_FFTADD4SCDNX16(a41, a42, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b41 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b42 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //--------------------

      // A <= a0
      BBE_MOVSAV(a50);

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a51, a52);

      t0 = BBE_MULNX16CPACKQ(a51, r3tw);
      t1 = BBE_MULNX16CPACKQ(a52, r3tw);

      // b0 <= a0 + a1 + a2
      b50 = BBE_FFTADDSSRNX16(a50, a51);
      b50 = BBE_FFTADDSSRNX16(b50, a52);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b51 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b52 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //
      // ( I9 x C16_12 )*( L18_3 x I8 )
      //

      a00 = b00; a01 = b10; a02 = b20; a03 = b30; a04 = b40; a05 = b50;
      a10 = b01; a11 = b11; a12 = b21; a13 = b31; a14 = b41; a15 = b51;
      a20 = b02; a21 = b12; a22 = b22; a23 = b32; a24 = b42; a25 = b52;

      Y0_va = BBE_ZALIGN();

      BBE_SAVNX16_XP(a00, Y0_va, Y0, 4 * 8);
      BBE_SAVNX16_XP(a01, Y0_va, Y0, 4 * 4);
      BBE_SAVNX16_XP(a02, Y0_va, Y0, 4 * 8);
      BBE_SAVNX16_XP(a03, Y0_va, Y0, 4 * 4);
      BBE_SAVNX16_XP(a04, Y0_va, Y0, 4 * 8);
      BBE_SAVNX16_XP(a05, Y0_va, Y0, 4 * 4);

      BBE_SAVNX16POS_FP(Y0_va, Y0);

      Y1_va = BBE_ZALIGN();

      BBE_SAVNX16_XP(a10, Y1_va, Y1, 4 * 8);
      BBE_SAVNX16_XP(a11, Y1_va, Y1, 4 * 4);
      BBE_SAVNX16_XP(a12, Y1_va, Y1, 4 * 8);
      BBE_SAVNX16_XP(a13, Y1_va, Y1, 4 * 4);
      BBE_SAVNX16_XP(a14, Y1_va, Y1, 4 * 8);
      BBE_SAVNX16_XP(a15, Y1_va, Y1, 4 * 4);

      BBE_SAVNX16POS_FP(Y1_va, Y1);

      Y2_va = BBE_ZALIGN();

      BBE_SAVNX16_XP(a20, Y2_va, Y2, 4 * 8);
      BBE_SAVNX16_XP(a21, Y2_va, Y2, 4 * 4);
      BBE_SAVNX16_XP(a22, Y2_va, Y2, 4 * 8);
      BBE_SAVNX16_XP(a23, Y2_va, Y2, 4 * 4);
      BBE_SAVNX16_XP(a24, Y2_va, Y2, 4 * 8);
      BBE_SAVNX16_XP(a25, Y2_va, Y2, 4 * 4);

      BBE_SAVNX16POS_FP(Y2_va, Y2);

      Y0 = (xb_vecNx16*)((uintptr_t)Y0 + BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 36 * 4);
      Y1 = (xb_vecNx16*)((uintptr_t)Y1 + BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 36 * 4);
      Y2 = (xb_vecNx16*)((uintptr_t)Y2 + BCNFFT_BUF_SIZE(108)*sizeof(complex_fract16) - 36 * 4);
    }
  }
} /* bcnfft108() */
#endif

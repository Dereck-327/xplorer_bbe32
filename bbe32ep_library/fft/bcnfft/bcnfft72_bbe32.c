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
DISCARD_FUN(void, bcnfft72,( void * restrict pScr,
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

void bcnfft72 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  //
  // DFT72 decomposition:
  //
  //   DFT72 = ( I6 x C16_12 )*( DFT6 x I16 )*  ...
  //           ...
  //           ( L24_6 x I4 )*( E12_16 x I6 )*( T18_6 x I4 )*( DFT3 x I24 )* ...
  //           ...
  //           ( C24_18 x I4 )*L96_48*L96_48*( L12_4 x I8 )*Tep72_18* ...
  //           ( I3 x DFT4 x I8 )*( L12_3 x I8 )*( I4 x E18_24 )
  //
  //   Tep72_18 = P*T72_18*P', where P = ( L12_3 x I8 )*( I4 x E18_24 )
  //

  const xb_vecNx16 *          X;
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          T;

  int l;

  NASSERT_ALIGN32(pScr);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Stage 1: ( C24_18 x I4 )*L96_48*L96_48*( L12_4 x I8 )*Tep72_18*
  //          ( I3 x DFT4 x I8 )*( L12_3 x I8 )*( I4 x E18_24 )

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)y;
  T = (const xb_vecNx16*)fft72_tw1;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;

    valign X1_va, Y0_va, Y1_va;

    xb_vecNx16 a00, a01, a02, a03;
    xb_vecNx16 a10, a11, a12, a13;
    xb_vecNx16 a20, a21, a22, a23;
    xb_vecNx16 a30, a31, a32;
    xb_vecNx16 b00, b01, b02, b03, b04;
    xb_vecNx16 b10, b11, b12, b13, b14;
    xb_vecNx16 b20, b21, b22, b23;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5, tw6, tw7, tw8, tw9;

    BBE_FFTWMODE(0x10);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 18 * 4);

  //  __Pragma("super_swp ii=22, unroll=1");
    for (l = 0; l<L; l++)
    {
      // 
      // Load twiddle factor table Tep72_18
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

      //
      // ( L12_3 x I8 )*( I4 x E18_24 )
      //

      a00 = BBE_LVNX16_I(X0, 0 * 4 * BBE_SIMD_WIDTH / 2);
      a10 = BBE_LVNX16_I(X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
      a20 = BBE_LVNX16_I(X0, 2 * 4 * BBE_SIMD_WIDTH / 2);

      X0 = (const xb_vecNx16*)((uintptr_t)X0 + BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16));

      X1_va = BBE_LAVNX16_PP(X1);

      BBE_LAVNX16_XP(a01, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a11, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a21, X1_va, X1, 2 * 4);
      BBE_LAVNX16_XP(a02, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a12, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a22, X1_va, X1, 2 * 4);
      BBE_LAVNX16_XP(a03, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a13, X1_va, X1, 8 * 4);
      BBE_LAVNX16_XP(a23, X1_va, X1, 2 * 4);

      X1 = (const xb_vecNx16*)((uintptr_t)X1 + BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 54 * 4);

      //
      // Tep72_18*( I3 x DFT4 x I8 )
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

      BBE_MOVSAV(a12);
      BBE_MOVSBV(a13);

      b10 = BBE_FFTADD4SABNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SABNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SABNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SABNX16(a10, a11, 3, 0);

      b11 = BBE_MULNX16CPACKQ(b11, tw4);
      b12 = BBE_MULNX16CPACKQ(b12, tw5);
      b13 = BBE_MULNX16CPACKQ(b13, tw6);

      BBE_MOVSCV(a22);
      BBE_MOVSDV(a23);

      b20 = BBE_FFTADD4SCDNX16(a20, a21, 0, 0);
      b21 = BBE_FFTADD4SCDNX16(a20, a21, 1, 0);
      b22 = BBE_FFTADD4SCDNX16(a20, a21, 2, 0);
      b23 = BBE_FFTADD4SCDNX16(a20, a21, 3, 0);

      b21 = BBE_MULNX16CPACKQ(b21, tw7);
      b22 = BBE_MULNX16CPACKQ(b22, tw8);
      b23 = BBE_MULNX16CPACKQ(b23, tw9);

      //
      // ( C24_18 x I4 )*L96_48*L96_48*( L12_4 x I8 )
      //

      // L12_4 x I8
      a00 = b00; a01 = b10; a02 = b20;
      a10 = b01; a11 = b11; a12 = b21;
      a20 = b02; a21 = b12; a22 = b22;
      a30 = b03; a31 = b13; a32 = b23;

      // L96_48
      BBE_DSELNX16I(b01, b00, a20, a00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b03, b02, a21, a01, BBE_DSELI_INTERLEAVE_2);

      b04 = BBE_SELNX16I(a22, a02, BBE_SELI_INTERLEAVE_2_LO);

      BBE_DSELNX16I(b11, b10, a30, a10, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b13, b12, a31, a11, BBE_DSELI_INTERLEAVE_2);

      b14 = BBE_SELNX16I(a32, a12, BBE_SELI_INTERLEAVE_2_LO);

      Y0_va = BBE_MOVUVR(b00);
      Y1_va = BBE_MOVUVR(b01);

      // ( C24_18 x I4 )*L96_48
      BBE_SVINTLARNX16_XP(b10, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(b11, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(b02, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b03, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(b12, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SVINTLARNX16_XP(b13, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
      BBE_SALIGNVRNX16_XP(b04, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(b13, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_SVINTLARNX16_XP(b14, Y0_va, Y, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 8 * 4 * BBE_SIMD_WIDTH / 2, 0);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 2: ( L24_6 x I4 )*( E12_16 x I6 )*( T18_6 x I4 )*( DFT3 x I24 )
  //
  // Radix-3 butterfly: b0 = a0 + a1 + a2
  //                    b1 = a1*r3tw - a2*r3tw + a0 - (a1+a2)/2
  //                    b2 = a2*r3tw - a1*r3tw + a0 - (a1+a2)/2
  //
  // where r3tw = sign*1j*(3^0.5)/2
  //

  X = (const xb_vecNx16*)y;
  Y = (xb_vecNx16*)x;
  T = (const xb_vecNx16*)bcnfft72_T18_6;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;

    valign Y0_va, Y1_va;

    xb_vecNx16 a00, a01, a02;
    xb_vecNx16 a10, a11, a12;
    xb_vecNx16 a20, a21, a22;
    xb_vecNx16 a30, a31, a32;
    xb_vecNx16 b00, b01, b02;
    xb_vecNx16 b10, b11, b12;
    xb_vecNx16 b20, b21, b22;

    xb_vecNx16 p0, p1, t0, t1;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5, tw6;

    xb_vecNx16 r3tw = BBE_MOVVA16C(0x91260000);

    t0 = 0;

    BBE_MOVSDV(t0);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 2 * 4 * BBE_SIMD_WIDTH / 2);

    for (l = 0; l<L; l++)
    {
      //
      // Load twiddle table Tp18_6
      //

      p0 = BBE_LVNX16_I(T, 0 * 4 * BBE_SIMD_WIDTH / 2);
      p1 = BBE_LVNX16_I(T, 1 * 4 * BBE_SIMD_WIDTH / 2);

      tw1 = BBE_SHFLNX16I(p0, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw2 = BBE_SHFLNX16I(p0, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw3 = BBE_SHFLNX16I(p0, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw4 = BBE_SHFLNX16I(p0, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw5 = BBE_SHFLNX16I(p1, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw6 = BBE_SHFLNX16I(p1, BBE_SHFLI_REP_2X4_OFFSET_1);

      //
      // ( T18_6 x I4 )*( DFT3 x I24 ) ==
      // ( L9_3 x I8 )*( Tp18_6 x I4 )*( I3 x DFT3 x I8 )*( L9_3 x I8 )
      //

      a00 = BBE_LVNX16_I(X0, 0 * 4 * BBE_SIMD_WIDTH / 2);
      a10 = BBE_LVNX16_I(X0, 1 * 4 * BBE_SIMD_WIDTH / 2);

      a01 = BBE_LVNX16_I(X0, 3 * 4 * BBE_SIMD_WIDTH / 2);
      a11 = BBE_LVNX16_I(X0, 4 * 4 * BBE_SIMD_WIDTH / 2);
      a21 = BBE_LVNX16_I(X0, 5 * 4 * BBE_SIMD_WIDTH / 2);
      a02 = BBE_LVNX16_I(X0, 6 * 4 * BBE_SIMD_WIDTH / 2);
      a12 = BBE_LVNX16_I(X0, 7 * 4 * BBE_SIMD_WIDTH / 2);
      a22 = BBE_LVNX16_I(X0, 8 * 4 * BBE_SIMD_WIDTH / 2);

      X0 = (const xb_vecNx16*)((uintptr_t)X0 + BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16));

      // A <= a0
      BBE_LVA_IP(X1, 4 * BBE_SIMD_WIDTH / 2);
      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a01, a02);

      t0 = BBE_MULNX16CPACKQ(a01, r3tw);
      t1 = BBE_MULNX16CPACKQ(a02, r3tw);

      // b0 <= a0 + a1 + a2
      b00 = BBE_ADDNX16(a00, a01);
      b00 = BBE_ADDNX16(b00, a02);

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
      b10 = BBE_ADDNX16(a10, a11);
      b10 = BBE_ADDNX16(b10, a12);

      // b1 <= a1*r3tw - a2*r3tw + A - B
      b11 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b12 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      b11 = BBE_MULNX16CPACKQ(b11, tw3);
      b12 = BBE_MULNX16CPACKQ(b12, tw4);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 0 * 4 * BBE_SIMD_WIDTH / 2);

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

      b21 = BBE_MULNX16CPACKQ(b21, tw5);
      b22 = BBE_MULNX16CPACKQ(b22, tw6);

      // L9_3 x I8
      a00 = b00; a01 = b10; a02 = b20;
      a10 = b01; a11 = b11; a12 = b21;
      a20 = b02; a21 = b12; a22 = b22;

      //
      // ( L24_6 x I4 )*( E12_16 x I6 )
      //
      a30 = 0; a31 = 0; a32 = 0;

      Y0_va = BBE_MOVUVR(a00);
      Y1_va = BBE_MOVUVR(a20);

      BBE_SVINTLARNX16_XP(a10, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a30, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a01, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a21, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a11, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a31, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a02, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a22, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a12, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(a32, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a12, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);

      BBE_SALIGNVRNX16_XP(a32, Y1_va, Y, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 11 * 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 3: ( I6 x C16_12 )*( DFT6 x I16 )
  //
  // We use the Prime-Factor FFT algorithm for the radix-6 butterfly:
  //
  //         |1 0 0 0 0 0|                                 |1 0 0 0 0 0|
  //         |0 0 0 0 1 0|                                 |0 0 1 0 0 0|
  // DFT6 =  |0 0 1 0 0 0| * ( DFT2 x I3 )*( I2 x DFT3 ) * |0 0 0 0 1 0|
  //         |0 0 0 1 0 0|                                 |0 0 0 1 0 0|
  //         |0 1 0 0 0 0|                                 |0 0 0 0 0 1|
  //         |0 0 0 0 0 1|                                 |0 1 0 0 0 0|
  //               P2                                            P1
  //
  // After substituting DFT6 we obtain:
  //   ( I6 x C16_12 )*( DFT6 x I16 ) ->
  //   ( I6 x C16_12 )*( P2 x I16 )*( DFT2 x I48 )*( I2 x L6_3 x I8 )*
  //   ( I4 x DFT3 x I8 )*( I2 x L6_2 x I8 )*( P1 x I16 )
  //

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)y;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;

    valign Y_va;

    xb_vecNx16 a00, a01, a02, a03, a04, a05;
    xb_vecNx16 a10, a11, a12, a13, a14, a15;
    xb_vecNx16 a20, a21, a22;
    xb_vecNx16 a31, a32;
    xb_vecNx16 b00, b01, b02, b03, b04, b05;
    xb_vecNx16 b10, b11, b12, b13, b14, b15;
    xb_vecNx16 b20, b21, b22;
    xb_vecNx16 b30, b31, b32;

    xb_vecNx16 t0, t1;

    xb_vecNx16 r3tw;

    t0 = 0;

    BBE_MOVSDV(t0);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 7 * 4 * BBE_SIMD_WIDTH / 2);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 7 * 4 * BBE_SIMD_WIDTH / 2);

    Y_va = BBE_ZALIGN();

    for (l = 0; l<L; l++)
    {
      r3tw = BBE_LVNX16_I((const xb_vecNx16*)bcnfft72_r3tw, 0);

      //
      // ( I2 x L6_2 x I8 )*( P1 x I16 )
      //

      a00 = BBE_LVNX16_I(X0, 0 * 4 * BBE_SIMD_WIDTH / 2);
      a01 = BBE_LVNX16_I(X0, 4 * 4 * BBE_SIMD_WIDTH / 2);
      a02 = BBE_LVNX16_I(X0, 8 * 4 * BBE_SIMD_WIDTH / 2);
      a10 = BBE_LVNX16_I(X0, 1 * 4 * BBE_SIMD_WIDTH / 2);
      a11 = BBE_LVNX16_I(X0, 5 * 4 * BBE_SIMD_WIDTH / 2);
      a12 = BBE_LVNX16_I(X0, 9 * 4 * BBE_SIMD_WIDTH / 2);
      a20 = BBE_LVNX16_I(X0, 6 * 4 * BBE_SIMD_WIDTH / 2);
      a21 = BBE_LVNX16_I(X0, 10 * 4 * BBE_SIMD_WIDTH / 2);
      a22 = BBE_LVNX16_I(X0, 2 * 4 * BBE_SIMD_WIDTH / 2);
      a31 = BBE_LVNX16_I(X0, 11 * 4 * BBE_SIMD_WIDTH / 2);
      a32 = BBE_LVNX16_I(X0, 3 * 4 * BBE_SIMD_WIDTH / 2);

      X0 = (const xb_vecNx16*)XT_ADD((uintptr_t)X0, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16));

      //
      // I4 x DFT3 x I8
      //
      // Radix-3 butterfly: b0 = a0 + a1 + a2
      //                    b1 = a1*r3tw - a2*r3tw + a0 - (a1+a2)/2
      //                    b2 = a2*r3tw - a1*r3tw + a0 - (a1+a2)/2
      //
      // where r3tw = sign*1j*(3^0.5)/2
      //

      // A,C <= a0
      BBE_MOVSAV(a00);
      BBE_MOVSCV(a00);

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
      BBE_MOVSAV(a10);
      BBE_MOVSCV(a10);

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
      BBE_MOVSAV(a20);
      BBE_MOVSCV(a20);

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

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16));
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16));

      // B <= (a1+a2)/2
      BBE_FFTAVGNX16SB(a31, a32);

      t0 = BBE_MULNX16CPACKQ(a31, r3tw);
      t1 = BBE_MULNX16CPACKQ(a32, r3tw);

      // b0 <= a1 + a2 + C
      b30 = BBE_FFTADD4SCDNX16(a31, a32, 0, 0);
      // b1 <= a1*r3tw - a2*r3tw + A - B
      b31 = BBE_FFTADD4SABNX16(t0, t1, 2, 0);
      // b2 <- a2*r3tw - a1*r3tw + A - B
      b32 = BBE_FFTADD4SABNX16(t1, t0, 2, 0);

      //
      // ( DFT2 x I48 )*( I2 x L6_3 x I8 )
      //

      a00 = b00; a01 = b10; a02 = b01; a03 = b11; a04 = b02; a05 = b12;
      a10 = b20; a11 = b30; a12 = b21; a13 = b31; a14 = b22; a15 = b32;

      b00 = BBE_FFTADDSSRNX16(a00, a10); // 0
      b01 = BBE_FFTADDSSRNX16(a01, a11);
      b02 = BBE_FFTADDSSRNX16(a02, a12); // 1
      b03 = BBE_FFTADDSSRNX16(a03, a13);
      b04 = BBE_FFTADDSSRNX16(a04, a14); // 2
      b05 = BBE_FFTADDSSRNX16(a05, a15);

      b10 = BBE_FFTSUBSSRNX16(a00, a10); // 3
      b11 = BBE_FFTSUBSSRNX16(a01, a11);
      b12 = BBE_FFTSUBSSRNX16(a02, a12); // 4
      b13 = BBE_FFTSUBSSRNX16(a03, a13);
      b14 = BBE_FFTSUBSSRNX16(a04, a14); // 5
      b15 = BBE_FFTSUBSSRNX16(a05, a15);

      //
      // ( I6 x C16_12 )*( P2 x I16 )
      //

      BBE_SAVNX16_XP(b00, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b01, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b12, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b13, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b04, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b05, Y_va, Y, 4 * 4);

      BBE_SAVNX16_XP(b10, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b11, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b02, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b03, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b14, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b15, Y_va, Y, 4 * 4);

      Y = (xb_vecNx16*)XT_ADD((uintptr_t)Y, BCNFFT_BUF_SIZE(72)*sizeof(complex_fract16) - 72 * 4);
    }
  }
} /* bcnfft72() */
#endif

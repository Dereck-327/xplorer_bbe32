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
DISCARD_FUN(void, bcnfft60,( void * restrict pScr,
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

void bcnfft60 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  //
  // DFT60 decomposition:
  //
  //   DFT60 = ( I5 x C4_3 x I4 )*( DFT5 x I16 )* ...
  //           ...
  //           ( I5 x E3_4 x I4 )*( C6_5 x I12 )*( L18_6 x I4 )* ...
  //           ( L9_3 x I8 )*( Tep15_5 x I4 )*( I3 x DFT3 x I8 )* ...
  //           ( L9_3 x I8 )*( I3 x E5_6 x I4 )* ...
  //           ...
  //           ( C16_15 x I4 )*L64_32*L64_32*( L8_4 x I8 )*Tep60_15* ...
  //           ( I2 x DFT4 x I8 )*( L8_2 x I8 )*( I4 x E15_16 )
  //
  //   Tep15_5 = P*T15_5*P', where P = ( L9_3 x I2 )*( I3 x E5_6 )
  //
  //   Tep60_15 = P*T60_15*P', where P = ( L8_2 x I8 )*( I4 x E15_16 )
  //

  const xb_vecNx16 *          X;
  xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          T;

  int l;

  NASSERT_ALIGN32(pScr);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  //----------------------------------------------------------------------------
  // Stage 1: ( C16_15 x I4 )*L64_32*L64_32*( L8_4 x I8 )*Tep60_15*
  //          ( I2 x DFT4 x I8 )*( L8_2 x I8 )*( I4 x E15_16 )
  //

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)y;
  T = (const xb_vecNx16*)fft60_tw1;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;

    valign X_va;

    xb_vecNx16 a00, a01, a02, a03;
    xb_vecNx16 a10, a11, a12, a13;
    xb_vecNx16 b00, b01, b02, b03;
    xb_vecNx16 b10, b11, b12, b13;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5, tw6;

    BBE_FFTWMODE(0x10);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 15 * 4);

    tw1 = BBE_LVNX16_I(T, 0 * 4 * BBE_SIMD_WIDTH / 2);
    tw4 = BBE_LVNX16_I(T, 3 * 4 * BBE_SIMD_WIDTH / 2);

    for (l = 0; l<L; l++)
    {
      //
      // ( L8_2 x I8 )*( I4 x E15_16 )
      //

      BBE_LVNX16_IP(a00, X0, +4 * BBE_SIMD_WIDTH / 2); //  0
      BBE_LVNX16_IP(a10, X0, -4 * BBE_SIMD_WIDTH / 2 + BCNFFT_BUF_SIZE(60)*(int)sizeof(complex_fract16)); //  8

      X_va = BBE_LAVNX16_PP(X1);

      BBE_LAVNX16_XP(a01, X_va, X1, 8 * 4); // 15
      BBE_LAVNX16_XP(a11, X_va, X1, 7 * 4); // 23
      BBE_LAVNX16_XP(a02, X_va, X1, 8 * 4); // 30
      BBE_LAVNX16_XP(a12, X_va, X1, 7 * 4); // 38
      BBE_LAVNX16_XP(a03, X_va, X1, 8 * 4); // 45
      BBE_LAVNX16_XP(a13, X_va, X1, 7 * 4); // 53

      X1 = (const xb_vecNx16*)((uintptr_t)X1 + BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 45 * 4);

      BBE_MOVSAV(a02);
      BBE_MOVSBV(a03);

      BBE_MOVSCV(a12);
      BBE_MOVSDV(a13);

      //
      // Tep60_15*( I2 x DFT4 x I8 )
      //

      b00 = BBE_FFTADD4SABNX16(a00, a01, 0, 0);
      b01 = BBE_FFTADD4SABNX16(a00, a01, 1, 0);
      b02 = BBE_FFTADD4SABNX16(a00, a01, 2, 0);
      b03 = BBE_FFTADD4SABNX16(a00, a01, 3, 0);

      b10 = BBE_FFTADD4SCDNX16(a10, a11, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a10, a11, 1, 0);
      b12 = BBE_FFTADD4SCDNX16(a10, a11, 2, 0);
      b13 = BBE_FFTADD4SCDNX16(a10, a11, 3, 0);

      tw2 = BBE_LVNX16_I(T, 1 * 4 * BBE_SIMD_WIDTH / 2);
      tw3 = BBE_LVNX16_I(T, 2 * 4 * BBE_SIMD_WIDTH / 2);

      tw5 = BBE_LVNX16_I(T, 4 * 4 * BBE_SIMD_WIDTH / 2);
      tw6 = BBE_LVNX16_I(T, 5 * 4 * BBE_SIMD_WIDTH / 2);

      b01 = BBE_MULNX16CPACKQ(b01, tw1);
      b02 = BBE_MULNX16CPACKQ(b02, tw2);
      b03 = BBE_MULNX16CPACKQ(b03, tw3);

      b11 = BBE_MULNX16CPACKQ(b11, tw4);
      b12 = BBE_MULNX16CPACKQ(b12, tw5);
      b13 = BBE_MULNX16CPACKQ(b13, tw6);

      //
      // ( C16_15 x I4 )*L64_32*L64_32*( L8_4 x I8 )
      //

      a00 = b00; a01 = b10;
      a02 = b01; a03 = b11;
      a10 = b02; a11 = b12;
      a12 = b03; a13 = b13;

      BBE_DSELNX16I(b01, b00, a10, a00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b03, b02, a11, a01, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b11, b10, a12, a02, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(b13, b12, a13, a03, BBE_DSELI_INTERLEAVE_2);

      BBE_DSELNX16I(a01, a00, b10, b00, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a03, a02, b11, b01, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a11, a10, b12, b02, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(a13, a12, b13, b03, BBE_DSELI_INTERLEAVE_2);

      BBE_SVNX16_I(a00, Y, 0 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a01, Y, 1 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a02, Y, 2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a03, Y, 3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a10, Y, 4 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a11, Y, 5 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a12, Y, 6 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_I(a13, Y, 7 * 4 * BBE_SIMD_WIDTH / 2);

      Y = (xb_vecNx16*)((uintptr_t)Y + BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16));
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 2: ( I5 x E3_4 x I4 )*( C6_5 x I12 )*( L18_6 x I4 )*
  //          ( L9_3 x I8 )*( Tep15_5 x I4 )*( I3 x DFT3 x I8 )*
  //          ( L9_3 x I8 )*( I3 x E5_6 x I4 )

  X = (const xb_vecNx16*)y;
  Y = (xb_vecNx16*)x;
  T = (const xb_vecNx16*)fft60_tw2;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    const xb_vecNx16 * X2;

    valign Y0_va, Y1_va;

    xb_vecNx16 a00, a01, a02, a10, a11, a12, a20, a21, a22;
    xb_vecNx16 b00, b01, b02, b10, b11, b12, b20, b21, b22;

    xb_vecNx16 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecNx16 tw1, tw2, tw3, tw4, tw5, tw6;

    xb_vecNx16 r3tw = BBE_MOVVA16C(0x91260000);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X2 = (const xb_vecNx16*)((uintptr_t)X + 1 * 4 * BBE_SIMD_WIDTH / 2);

    t0 = 0;

    BBE_MOVSDV(t0);

    for (l = 0; l<L; l++)
    {
      //
      // Load twiddle factor table Tep15_5
      //

      t0 = BBE_LVNX16_I(T, 0 * 4 * BBE_SIMD_WIDTH / 2);
      t1 = BBE_LVNX16_I(T, 1 * 4 * BBE_SIMD_WIDTH / 2);

      tw1 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw2 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_1);
      tw3 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_2);
      tw4 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4_OFFSET_3);
      tw5 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_0);
      tw6 = BBE_SHFLNX16I(t1, BBE_SHFLI_REP_2X4_OFFSET_1);

      //
      // ( L9_3 x I8 )*( I3 x E5_6 x I4 )
      //

      BBE_LVNX16_IP(t0, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t1, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t2, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t3, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t4, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t5, X0, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(t6, X0, 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(t7, X0, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 7 * 4 * BBE_SIMD_WIDTH / 2);

      a00 = t0;
      a10 = t1;
      a20 = t2;

      a01 = BBE_SELNX16I(t3, t2, BBE_SELI_ROTATE_RIGHT_8);
      a11 = BBE_SELNX16I(t4, t3, BBE_SELI_ROTATE_RIGHT_8);
      a21 = BBE_SHFLNX16I(t4, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

      a02 = t5;
      a12 = t6;
      a22 = t7;

      //
      // ( Tep15_5 x I4 )*( I3 x DFT3 x I8 )
      //
      // Radix-3 butterfly: b0 = a0 + a1 + a2
      //                    b1 = a1*r3tw - a2*r3tw + a0 - (a1+a2)/2
      //                    b2 = a2*r3tw - a1*r3tw + a0 - (a1+a2)/2
      //
      // where r3tw = sign*1j*(3^0.5)/2
      //

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

      b11 = BBE_MULNX16CPACKQ(b11, tw3);
      b12 = BBE_MULNX16CPACKQ(b12, tw4);

      //--------------------

      // A,C <= a0
      BBE_LVA_IP(X1, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVC_IP(X2, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 1 * 4 * BBE_SIMD_WIDTH / 2);
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

      //
      // ( I5 x E3_4 x I4 )*( C6_5 x I12 )*( L18_6 x I4 )*( L9_3 x I8 )
      //

      // L9_3 x I8
      a00 = b00; a01 = b10; a02 = b20;
      a10 = b01; a11 = b11; a12 = b21;
      a20 = b02; a21 = b12; a22 = b22;

      t0 = 0;

      Y0_va = BBE_MOVUVR(a00);
      Y1_va = BBE_MOVUVR(a20);

      // ( I5 x E3_4 x I4 )*( C6_5 x I12 )*( L18_6 x I4 )
      BBE_SVINTLARNX16_XP(a10, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(t0, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a01, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a21, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a11, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(t0, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SALIGNVRNX16_XP(a02, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SALIGNVRNX16_XP(a22, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVINTLARNX16_XP(a12, Y0_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
      BBE_SVINTLARNX16_XP(t0, Y1_va, Y, 4 * BBE_SIMD_WIDTH / 2, 1);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Stage 3: ( I5 x C4_3 x I4 )*( DFT5 x I16 )
  //
  // We use Rader decomposition for the radix-5 butterfly:
  //
  //                [1        ]   [1  0  0  0  0]   [1     1       ]   [1  0  0  0  0]   [1        ]
  //                [  1      ]   [0            ]   [1 -0.25       ]   [0            ]   [  1      ]
  // DFT(5, sign) = [    1    ] o [0 DFT(4, -1) ] o [        w1    ] o [0 DFT(4, -1) ] o [    1    ]
  //                [        1]   [0            ]   [          w2  ]   [0            ]   [        1]
  //                [      1  ]   [0            ]   [            w3]   [0            ]   [      1  ]
  //
  // Where w1 = sign * (-0.29389262614623656 + j*0.47552825814757676)
  //       w2 = 0.55901699437494751;
  //       w3 = -conj(w1);
  //
  // Note that two internal DFT4 are both forward (-1) regardless of transform sign!
  //

  X = (const xb_vecNx16*)x;
  Y = (xb_vecNx16*)y;

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;

    valign Y_va;

    xb_vecNx16 a00, a01, a10, a11, a20, a21, a30, a31, a40, a41;
    xb_vecNx16 b00, b01, b10, b11, b20, b21, b30, b31, b40, b41;

    xb_vecNx16 r5tw1 = BBE_MOVVA16C(0xc322259e); //  0.29388 - j*0.47552
    xb_vecNx16 r5tw2 = BBE_MOVVA16C(0x0000478e); //  0.55902
    xb_vecNx16 r5tw3 = BBE_MOVVA16C(0xc322da62); // -0.29388 - j*0.47552

    BBE_FFTWMODE(0x10);

    X0 = (const xb_vecNx16*)((uintptr_t)X + 0 * 4 * BBE_SIMD_WIDTH / 2);
    X1 = (const xb_vecNx16*)((uintptr_t)X + 6 * 4 * BBE_SIMD_WIDTH / 2);

    for (l = 0; l<L; l++)
    {
      //
      // DFT5 x I16
      //

      BBE_LVNX16_IP(a00, X0, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a10, X0, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_XP(a20, X0, -3 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a01, X0, +2 * 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a11, X0, +2 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(a21, X0, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 5 * 4 * BBE_SIMD_WIDTH / 2);

      BBE_LVB_IP(X1, +2 * 4 * BBE_SIMD_WIDTH / 2); // a30
      BBE_LVA_IP(X1, -1 * 4 * BBE_SIMD_WIDTH / 2); // a40
      BBE_LVD_IP(X1, +2 * 4 * BBE_SIMD_WIDTH / 2); // a31

      BBE_LVC_IP(X1, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 3 * 4 * BBE_SIMD_WIDTH / 2); // a41

      b00 = a00; b01 = a01;

      b10 = BBE_FFTADD4SABNX16(a10, a20, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a11, a21, 0, 0);
      b20 = BBE_FFTADD4SABNX16(a10, a20, 1, 0);
      b21 = BBE_FFTADD4SCDNX16(a11, a21, 1, 0);
      b30 = BBE_FFTADD4SABNX16(a10, a20, 2, 0);
      b31 = BBE_FFTADD4SCDNX16(a11, a21, 2, 0);
      b40 = BBE_FFTADD4SABNX16(a10, a20, 3, 0);
      b41 = BBE_FFTADD4SCDNX16(a11, a21, 3, 0);

      a00 = BBE_ADDSNX16(b00, b10);
      a01 = BBE_ADDSNX16(b01, b11);

      a10 = BBE_SRAINX16(b10, 2);
      a10 = BBE_SUBSNX16(b00, a10);
      a11 = BBE_SRAINX16(b11, 2);
      a11 = BBE_SUBSNX16(b01, a11);

      a20 = BBE_MULNX16CPACKQ(b20, r5tw1);
      a21 = BBE_MULNX16CPACKQ(b21, r5tw1);
      a30 = BBE_MULNX16CPACKQ(b30, r5tw2);
      a31 = BBE_MULNX16CPACKQ(b31, r5tw2);
      a40 = BBE_MULNX16CPACKQ(b40, r5tw3);
      a41 = BBE_MULNX16CPACKQ(b41, r5tw3);

      BBE_MOVSAV(a30);
      BBE_MOVSCV(a31);
      BBE_MOVSBV(a40);
      BBE_MOVSDV(a41);

      b00 = a00;
      b01 = a01;

      b10 = BBE_FFTADD4SABNX16(a10, a20, 0, 0);
      b11 = BBE_FFTADD4SCDNX16(a11, a21, 0, 0);
      b20 = BBE_FFTADD4SABNX16(a10, a20, 1, 0);
      b21 = BBE_FFTADD4SCDNX16(a11, a21, 1, 0);
      b40 = BBE_FFTADD4SABNX16(a10, a20, 2, 0);
      b41 = BBE_FFTADD4SCDNX16(a11, a21, 2, 0);
      b30 = BBE_FFTADD4SABNX16(a10, a20, 3, 0);
      b31 = BBE_FFTADD4SCDNX16(a11, a21, 3, 0);

      // 
      // I5 x C4_3 x I4
      //

      Y_va = BBE_ZALIGN();

      BBE_SAVNX16_XP(b00, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b01, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b10, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b11, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b20, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b21, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b30, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b31, Y_va, Y, 4 * 4);
      BBE_SAVNX16_XP(b40, Y_va, Y, 8 * 4);
      BBE_SAVNX16_XP(b41, Y_va, Y, 4 * 4);

      BBE_SAVNX16POS_FP(Y_va, Y);

      Y = (xb_vecNx16*)XT_ADD((uintptr_t)Y, BCNFFT_BUF_SIZE(60)*sizeof(complex_fract16) - 60 * 4);
    }
  }
} /* bcnfft60() */
#endif

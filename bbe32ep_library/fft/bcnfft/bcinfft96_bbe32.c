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
    Blockwise mixed radix inverse FFT on complex data, no data scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "fft_tw.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#if !(HAVE_FFT && 1)
DISCARD_FUN(void, bcinfft96,( void * restrict pScr,
                 complex_fract16 * restrict y,
                 complex_fract16 * restrict x,
                 int L ))
#else
/*-------------------------------------------------------------------------
Blockwise mixed radix inverse FFT on complex data, no data scaling
  
Description: These functions make inverse FFT on complex data of the following
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
    pScr[]      Scratch memory area of BCINFFT_SCRATCH_SIZE(N) bytes
  Input:
    S           Required input/output buffer size may exceed actual data size:
                S >= 2*N. Use BCINFFT_BUF_SIZE(N) macro to determine the minimum
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

void bcinfft96 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  int l;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 *          X0;
  xb_vecNx16 *          X1;
  xb_vecNx16 *          X2;

  const xb_vecNx16 *          C0 = (xb_vecNx16 *)fft96_tw1;
  const xb_vecNx16 *          C1 = (xb_vecNx16 *)fft96_tw2;
  xb_vecNx16 x0, x1, x2, x3, y0, y1, y2, y3;
  xb_vecNx16 x4, x5, x6, x7, x8, x9, x10, x11;
  xb_vecNx16 y4, y5, y6, y7, y8, y9, y10, y11;
  xb_vecNx16 tw96_24_0, tw96_24_1, tw96_24_2, tw96_24_3, tw96_24_4, tw96_24_5;
  xb_vecNx16 tw96_24_6, tw96_24_7, tw96_24_8;
  xb_vecNx16 tw2_0, tw2_1, tw2_2, tw3;
  xb_vecNx16 tw24_6_0, tw24_6_1, tw24_6_2, tw24_6_3, tw24_6_4, tw24_6_5, tw24_6_6, tw24_6_7, tw24_6_8;
  valign x_align;
  valign Y0_va, Y1_va;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  //DFT96 = (DFT6 x I16)*(L24_6 x I4)*(T24_6*(DFT4 x I6) x I4)*L96_24*T96_24*(DFT4 x I24)
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  X1 = (xb_vecNx16*)(x + BBE_SIMD_WIDTH/2);
  X2 = (xb_vecNx16*)(x + BBE_SIMD_WIDTH);
  x_align = BBE_LAVNX16_PP(X0);
  //----------------------------------------------------------------------------

  BBE_FFTWMODE(0x00);
  // DFT4 x I24  
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(x0, X0, 3 * 2 * BBE_SIMD_WIDTH); // 0
    BBE_LVNX16_IP(x1, X0, 3 * 2 * BBE_SIMD_WIDTH); // 3
    BBE_LVNX16_IP(x2, X0, 3 * 2 * BBE_SIMD_WIDTH); // 6
    BBE_LVNX16_IP(x3, X0, 3 * 2 * BBE_SIMD_WIDTH); // 9
    BBE_LVNX16_IP(x4, X1, 3 * 2 * BBE_SIMD_WIDTH); // 1
    BBE_LVNX16_IP(x5, X1, 3 * 2 * BBE_SIMD_WIDTH); // 4
    BBE_LVNX16_IP(x6, X1, 3 * 2 * BBE_SIMD_WIDTH); // 7
    BBE_LVNX16_IP(x7, X1, 3 * 2 * BBE_SIMD_WIDTH); // 10
    BBE_LVNX16_IP(x8, X2, 3 * 2 * BBE_SIMD_WIDTH); // 2
    BBE_LVNX16_IP(x9, X2, 3 * 2 * BBE_SIMD_WIDTH); // 5
    BBE_LVNX16_IP(x10, X2, 3 * 2 * BBE_SIMD_WIDTH); // 8
    BBE_LVNX16_IP(x11, X2, 3 * 2 * BBE_SIMD_WIDTH); // 11

    // Load twiddles and constants into the registers file
    tw96_24_0 = BBE_LVNX16_I(C0, 0 * 2 * BBE_SIMD_WIDTH);
    tw96_24_1 = BBE_LVNX16_I(C0, 1 * 2 * BBE_SIMD_WIDTH);
    tw96_24_2 = BBE_LVNX16_I(C0, 2 * 2 * BBE_SIMD_WIDTH);
    tw96_24_3 = BBE_LVNX16_I(C0, 3 * 2 * BBE_SIMD_WIDTH);
    tw96_24_4 = BBE_LVNX16_I(C0, 4 * 2 * BBE_SIMD_WIDTH);
    tw96_24_5 = BBE_LVNX16_I(C0, 5 * 2 * BBE_SIMD_WIDTH);
    tw96_24_6 = BBE_LVNX16_I(C0, 6 * 2 * BBE_SIMD_WIDTH);
    tw96_24_7 = BBE_LVNX16_I(C0, 7 * 2 * BBE_SIMD_WIDTH);
    tw96_24_8 = BBE_LVNX16_I(C0, 8 * 2 * BBE_SIMD_WIDTH);
    BBE_MOVSAV(x2);
    BBE_MOVSBV(x3);
    BBE_MOVSCV(x6);
    BBE_MOVSDV(x7);
    y0 = BBE_FFTADD4SABNX16(x0, x1, 0, 0);
    y1 = BBE_FFTADD4SABNX16(x0, x1, 1, 0);
    y2 = BBE_FFTADD4SABNX16(x0, x1, 2, 0);
    y3 = BBE_FFTADD4SABNX16(x0, x1, 3, 0);

    y4 = BBE_FFTADD4SCDNX16(x4, x5, 0, 0);
    y5 = BBE_FFTADD4SCDNX16(x4, x5, 1, 0);
    y6 = BBE_FFTADD4SCDNX16(x4, x5, 2, 0);
    y7 = BBE_FFTADD4SCDNX16(x4, x5, 3, 0);

    BBE_MOVSAV(x10);
    BBE_MOVSBV(x11);

    y8 = BBE_FFTADD4SABNX16(x8, x9, 0, 0);
    y9 = BBE_FFTADD4SABNX16(x8, x9, 1, 0);
    y10 = BBE_FFTADD4SABNX16(x8, x9, 2, 0);
    y11 = BBE_FFTADD4SABNX16(x8, x9, 3, 0);

    y1 = BBE_MULNX16JPACKQ(y1, tw96_24_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw96_24_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw96_24_2);

    y5 = BBE_MULNX16JPACKQ(y5, tw96_24_3);
    y6 = BBE_MULNX16JPACKQ(y6, tw96_24_4);
    y7 = BBE_MULNX16JPACKQ(y7, tw96_24_5);

    y9 = BBE_MULNX16JPACKQ(y9, tw96_24_6);
    y10 = BBE_MULNX16JPACKQ(y10, tw96_24_7);
    y11 = BBE_MULNX16JPACKQ(y11, tw96_24_8);

    BBE_DSELNX16I(x1, x0, y2, y0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x3, x2, y3, y1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x5, x4, y6, y4, BBE_DSELI_INTERLEAVE_2);
    x6 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_LO);
    x7 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_HI);
    x8 = BBE_SELNX16I(y10, y8, BBE_SELI_INTERLEAVE_2_LO);
    x9 = BBE_SELNX16I(y10, y8, BBE_SELI_INTERLEAVE_2_HI);
    x10 = BBE_SELNX16I(y11, y9, BBE_SELI_INTERLEAVE_2_LO);
    x11 = BBE_SELNX16I(y11, y9, BBE_SELI_INTERLEAVE_2_HI);

    Y0_va = BBE_MOVUVR(x0);
    Y1_va = BBE_MOVUVR(x1);

    BBE_SVINTLARNX16_XP(x2, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(x3, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(x4, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(x5, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(x6, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(x7, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(x8, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(x9, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(x10, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(x11, Y1_va, Y, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(x10, Y0_va, Y, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(x11, Y1_va, Y, BCNFFT_BUF_SIZE(96)*sizeof(complex_fract16) - 11 * 4 * BBE_SIMD_WIDTH / 2);
  }
  //********************************************************************
  __Pragma("no_reorder");
  Y = (xb_vecNx16*)x;
  X0 = (xb_vecNx16*)y;
  X1 = (xb_vecNx16*)(y + BBE_SIMD_WIDTH/2);
  X2 = (xb_vecNx16*)(y + BBE_SIMD_WIDTH);
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(x0, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x6, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x7, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x8, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x9, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x10, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x11, X2, 3 * 2 * BBE_SIMD_WIDTH);

    // Load twiddles and constants into the registers file
    BBE_LVNX16_IP(tw2_0, C1, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw2_1, C1, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw2_2, C1, -2 * 2 * BBE_SIMD_WIDTH);
    tw24_6_0 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw24_6_1 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw24_6_2 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_2);


    tw24_6_3 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_3);
    tw24_6_4 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw24_6_5 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_1);

    tw24_6_6 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_2);
    tw24_6_7 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_3);
    tw24_6_8 = BBE_SHFLNX16I(tw2_2, BBE_SHFLI_REP_2X4_OFFSET_0);

    BBE_MOVSAV(x2);
    BBE_MOVSBV(x3);
    BBE_MOVSCV(x6);
    BBE_MOVSDV(x7);

    y0 = BBE_FFTADD4SABNX16(x0, x1, 0, 0);
    y1 = BBE_FFTADD4SABNX16(x0, x1, 1, 0);
    y2 = BBE_FFTADD4SABNX16(x0, x1, 2, 0);
    y3 = BBE_FFTADD4SABNX16(x0, x1, 3, 0);

    y4 = BBE_FFTADD4SCDNX16(x4, x5, 0, 0);
    y5 = BBE_FFTADD4SCDNX16(x4, x5, 1, 0);
    y6 = BBE_FFTADD4SCDNX16(x4, x5, 2, 0);
    y7 = BBE_FFTADD4SCDNX16(x4, x5, 3, 0);

    BBE_MOVSAV(x10);
    BBE_MOVSBV(x11);

    y8 = BBE_FFTADD4SABNX16(x8, x9, 0, 0);
    y9 = BBE_FFTADD4SABNX16(x8, x9, 1, 0);
    y10 = BBE_FFTADD4SABNX16(x8, x9, 2, 0);
    y11 = BBE_FFTADD4SABNX16(x8, x9, 3, 0);

    y1 = BBE_MULNX16JPACKQ(y1, tw24_6_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw24_6_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw24_6_2);

    y5 = BBE_MULNX16JPACKQ(y5, tw24_6_3);
    y6 = BBE_MULNX16JPACKQ(y6, tw24_6_4);
    y7 = BBE_MULNX16JPACKQ(y7, tw24_6_5);

    y9 = BBE_MULNX16JPACKQ(y9, tw24_6_6);
    y10 = BBE_MULNX16JPACKQ(y10, tw24_6_7);
    y11 = BBE_MULNX16JPACKQ(y11, tw24_6_8);

    Y0_va = BBE_MOVUVR(y0);
    Y1_va = BBE_MOVUVR(y2);

    BBE_SVINTLARNX16_XP(y1, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SVINTLARNX16_XP(y3, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SALIGNVRNX16_XP(y4, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y6, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(y5, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SVINTLARNX16_XP(y7, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SALIGNVRNX16_XP(y8, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y10, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(y9, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SVINTLARNX16_XP(y11, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2, 1);
    BBE_SALIGNVRNX16_XP(y9, Y0_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y11, Y1_va, Y, +1 * 4 * BBE_SIMD_WIDTH / 2);
  }
  //********************************************************************
  __Pragma("no_reorder");
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  X1 = (xb_vecNx16*)(x + BBE_SIMD_WIDTH/2);
  X2 = (xb_vecNx16*)(x + 3 * BBE_SIMD_WIDTH);
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);

  for (l = 0; l<L; l++)
  {
    tw3 = BBE_MOVVA16C(0x91260000);

    BBE_LVNX16_IP(y0, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y3, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y5, X0, -9 * 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(y6, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y7, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y8, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y9, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y10, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y11, X0, 2 * BBE_SIMD_WIDTH);

    x0 = BBE_FFTADDSSRNX16(y0, y3);
    x1 = BBE_FFTSUBSSRNX16(y0, y3);
    x2 = BBE_FFTADDSSRNX16(y4, y1);
    x3 = BBE_FFTSUBSSRNX16(y4, y1);
    x4 = BBE_FFTADDSSRNX16(y2, y5);
    x5 = BBE_FFTSUBSSRNX16(y2, y5);

    x6 = BBE_FFTADDSSRNX16(y6, y9);
    x7 = BBE_FFTSUBSSRNX16(y6, y9);
    x8 = BBE_FFTADDSSRNX16(y10, y7);
    x9 = BBE_FFTSUBSSRNX16(y10, y7);
    x10 = BBE_FFTADDSSRNX16(y8, y11);
    x11 = BBE_FFTSUBSSRNX16(y8, y11);

    BBE_MOVSAV(x0);
    BBE_MOVSCV(x0);

    /* t0 = t1 + t2 + C + D; D==0; */
    BBE_FFTAVGNX16SB(x2, x4);
    y0 = BBE_FFTADD4SCDNX16(x2, x4, 0, 0);

    x2 = BBE_MULNX16JPACKQ(x2, tw3);
    x4 = BBE_MULNX16JPACKQ(x4, tw3);
    y1 = BBE_FFTADD4SABNX16(x2, x4, 2, 0);
    y2 = BBE_FFTADD4SABNX16(x4, x2, 2, 0);

    BBE_MOVSAV(x1);
    BBE_MOVSCV(x1);

    BBE_FFTAVGNX16SB(x3, x5);
    y3 = BBE_FFTADD4SCDNX16(x3, x5, 0, 0);

    x3 = BBE_MULNX16JPACKQ(x3, tw3);
    x5 = BBE_MULNX16JPACKQ(x5, tw3);
    y4 = BBE_FFTADD4SABNX16(x3, x5, 2, 0);
    y5 = BBE_FFTADD4SABNX16(x5, x3, 2, 0);

    BBE_MOVSAV(x6);
    BBE_MOVSCV(x6);

    /* t0 = t1 + t2 + C + D; D==0; */
    BBE_FFTAVGNX16SB(x8, x10);
    y6 = BBE_FFTADD4SCDNX16(x8, x10, 0, 0);

    x8 = BBE_MULNX16JPACKQ(x8, tw3);
    x10 = BBE_MULNX16JPACKQ(x10, tw3);
    y7 = BBE_FFTADD4SABNX16(x8, x10, 2, 0);
    y8 = BBE_FFTADD4SABNX16(x10, x8, 2, 0);

    BBE_MOVSAV(x7);
    BBE_MOVSCV(x7);

    BBE_FFTAVGNX16SB(x9, x11);
    y9 = BBE_FFTADD4SCDNX16(x9, x11, 0, 0);

    x9 = BBE_MULNX16JPACKQ(x9, tw3);
    x11 = BBE_MULNX16JPACKQ(x11, tw3);
    y10 = BBE_FFTADD4SABNX16(x9, x11, 2, 0);
    y11 = BBE_FFTADD4SABNX16(x11, x9, 2, 0);

    BBE_SVNX16_IP(y0, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y5, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y1, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y3, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y2, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(y4, Y, -9 * 2 * BBE_SIMD_WIDTH);

    BBE_SVNX16_IP(y6, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y11, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y7, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y9, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y8, Y, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y10, Y, 1 * 2 * BBE_SIMD_WIDTH);
  }
} /* bcinfft96() */
#endif

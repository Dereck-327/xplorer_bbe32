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
DISCARD_FUN(void, bcinfft120,( void * restrict pScr,
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

void bcinfft120 ( void * restrict pScr, complex_fract16 * restrict _y, complex_fract16 * restrict _x, int L )
{
  int16_t * restrict y=(int16_t *)_y;
  int16_t * restrict x=(int16_t *)_x;
  int l;
  xb_vecNx16 * restrict Y0;
  xb_vecNx16 * restrict Y1;
  xb_vecNx16 * restrict Y2;
  xb_vecNx16 *          X0;
  xb_vecNx16 *          X1;
  xb_vecNx16 *          X2;
  xb_vecNx16 *          X3;

  const xb_vecNx16 *          C0 = (xb_vecNx16 *)fft120_tw1;
  const xb_vecNx16 *          C1 = (xb_vecNx16 *)fft120_tw2;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
  xb_vecNx16 x8, x9, x10, x11, x12, x13, x14, x15;
  xb_vecNx16 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecNx16 y8, y9, y10, y11, y12, y13, y14, y15;
  xb_vecNx16 tw120_30_0, tw120_30_1, tw120_30_2, tw120_30_3, tw120_30_4, tw120_30_5;
  xb_vecNx16 tw120_30_6, tw120_30_7, tw120_30_8, tw120_30_9, tw120_30_10, tw120_30_11;
  xb_vecNx16 tw2_0, tw2_1, tw2_2, tw2_3;
  xb_vecNx16 tw30_5_0, tw30_5_1, tw30_5_2, tw30_5_3, tw30_5_4;
  xb_vecNx16 tw3;
  xb_vecNx16 tw5_0, tw5_1, tw5_2;
  valign Y0_va, Y1_va;
  valign x1_align, x2_align, x3_align;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  //DFT120 = (DFT5 x I24)*(L30_5 x I4)*(T30_5*(DFT6 x I5) x I4)*L120_30*T120_30*(DFT4 x I30)
  X0 = (xb_vecNx16*)x;
  X1 = (xb_vecNx16*)(x + 3 * BBE_SIMD_WIDTH + 12);
  X2 = (xb_vecNx16*)(x + 7 * BBE_SIMD_WIDTH + 8);
  X3 = (xb_vecNx16*)(x + 11 * BBE_SIMD_WIDTH);
  Y0 = (xb_vecNx16*)y;
  BBE_FFTWMODE(0x00);
  
  x2_align = BBE_LANX16_PP(X2);
  x3_align = BBE_LANX16_PP(X3);
  for (l = 0; l<L; l++)
  {
    x1_align = BBE_LANX16_PP(X1);
    BBE_LVNX16_IP(y0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y8, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y12, X0, 12 * 2 * BBE_SIMD_WIDTH);

    BBE_LANX16_IP(y1, x1_align, X1);
    BBE_LANX16_IP(y5, x1_align, X1);
    BBE_LANX16_IP(y9, x1_align, X1);
    BBE_LANX16_IP(y13, x1_align, X1);

    BBE_LANX16_IP(y2, x2_align, X2);
    BBE_LANX16_IP(y6, x2_align, X2);
    BBE_LANX16_IP(y10, x2_align, X2);
    BBE_LANX16_IP(y14, x2_align, X2);

    x0 = BBE_LVNX16_I(X3, 0 * 2 * BBE_SIMD_WIDTH);
    x1 = BBE_LVNX16_I(X3, 1 * 2 * BBE_SIMD_WIDTH);
    x2 = BBE_LVNX16_I(X3, 2 * 2 * BBE_SIMD_WIDTH);
    x3 = BBE_LVNX16_I(X3, 3 * 2 * BBE_SIMD_WIDTH);

    X3 = (xb_vecNx16*)XT_ADD((uintptr_t)X3, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16));
    y3 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_4);
    y7 = BBE_SELNX16I(x2, x1, BBE_SELI_ROTATE_RIGHT_4);
    y11 = BBE_SELNX16I(x3, x2, BBE_SELI_ROTATE_RIGHT_4);
    y15 = BBE_SELNX16I(x3, x3, BBE_SELI_ROTATE_RIGHT_4);

    X1 = (xb_vecNx16*)XT_ADD((uintptr_t)X1, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 4 * 2 * BBE_SIMD_WIDTH);
    X2 = (xb_vecNx16*)XT_ADD((uintptr_t)X2, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 4 * 2 * BBE_SIMD_WIDTH);
    
    x2_align = BBE_LANX16_PP(X2);
    // Load twiddles and constants into the registers file
    tw120_30_0 = BBE_LVNX16_I(C0, 0 * 2 * BBE_SIMD_WIDTH);
    tw120_30_1 = BBE_LVNX16_I(C0, 1 * 2 * BBE_SIMD_WIDTH);
    tw120_30_2 = BBE_LVNX16_I(C0, 2 * 2 * BBE_SIMD_WIDTH);
    tw120_30_3 = BBE_LVNX16_I(C0, 3 * 2 * BBE_SIMD_WIDTH);
    tw120_30_4 = BBE_LVNX16_I(C0, 4 * 2 * BBE_SIMD_WIDTH);
    tw120_30_5 = BBE_LVNX16_I(C0, 5 * 2 * BBE_SIMD_WIDTH);
    tw120_30_6 = BBE_LVNX16_I(C0, 6 * 2 * BBE_SIMD_WIDTH);
    tw120_30_7 = BBE_LVNX16_I(C0, 7 * 2 * BBE_SIMD_WIDTH);
    tw120_30_8 = BBE_LVNX16_I(C0, 8 * 2 * BBE_SIMD_WIDTH);
    tw120_30_9 = BBE_LVNX16_I(C0, 9 * 2 * BBE_SIMD_WIDTH);
    tw120_30_10 = BBE_LVNX16_I(C0, 10 * 2 * BBE_SIMD_WIDTH);
    tw120_30_11 = BBE_LVNX16_I(C0, 11 * 2 * BBE_SIMD_WIDTH);

    BBE_MOVSAV(y2);
    BBE_MOVSBV(y3);
    BBE_MOVSCV(y6);
    BBE_MOVSDV(y7);
    x0 = BBE_FFTADD4SABNX16(y0, y1, 0, 0);
    x1 = BBE_FFTADD4SABNX16(y0, y1, 1, 0);
    x2 = BBE_FFTADD4SABNX16(y0, y1, 2, 0);
    x3 = BBE_FFTADD4SABNX16(y0, y1, 3, 0);

    x4 = BBE_FFTADD4SCDNX16(y4, y5, 0, 0);
    x5 = BBE_FFTADD4SCDNX16(y4, y5, 1, 0);
    x6 = BBE_FFTADD4SCDNX16(y4, y5, 2, 0);
    x7 = BBE_FFTADD4SCDNX16(y4, y5, 3, 0);

    BBE_MOVSAV(y10);
    BBE_MOVSBV(y11);
    BBE_MOVSCV(y14);
    BBE_MOVSDV(y15);
    x8 = BBE_FFTADD4SABNX16(y8, y9, 0, 0);
    x9 = BBE_FFTADD4SABNX16(y8, y9, 1, 0);
    x10 = BBE_FFTADD4SABNX16(y8, y9, 2, 0);
    x11 = BBE_FFTADD4SABNX16(y8, y9, 3, 0);

    x12 = BBE_FFTADD4SCDNX16(y12, y13, 0, 0);
    x13 = BBE_FFTADD4SCDNX16(y12, y13, 1, 0);
    x14 = BBE_FFTADD4SCDNX16(y12, y13, 2, 0);
    x15 = BBE_FFTADD4SCDNX16(y12, y13, 3, 0);

    x1 = BBE_MULNX16JPACKQ(x1, tw120_30_0);
    x2 = BBE_MULNX16JPACKQ(x2, tw120_30_1);
    x3 = BBE_MULNX16JPACKQ(x3, tw120_30_2);

    x5 = BBE_MULNX16JPACKQ(x5, tw120_30_3);
    x6 = BBE_MULNX16JPACKQ(x6, tw120_30_4);
    x7 = BBE_MULNX16JPACKQ(x7, tw120_30_5);

    x9 = BBE_MULNX16JPACKQ(x9, tw120_30_6);
    x10 = BBE_MULNX16JPACKQ(x10, tw120_30_7);
    x11 = BBE_MULNX16JPACKQ(x11, tw120_30_8);

    x13 = BBE_MULNX16JPACKQ(x13, tw120_30_9);
    x14 = BBE_MULNX16JPACKQ(x14, tw120_30_10);
    x15 = BBE_MULNX16JPACKQ(x15, tw120_30_11);

    BBE_DSELNX16I(y2, y0, x2, x0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y3, y1, x3, x1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y6, y4, x6, x4, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y7, y5, x7, x5, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y10, y8, x10, x8, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y11, y9, x11, x9, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y14, y12, x14, x12, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y15, y13, x15, x13, BBE_DSELI_INTERLEAVE_2);

    Y0_va = BBE_MOVUVR(y0);
    Y1_va = BBE_MOVUVR(y2);

    BBE_SVINTLARNX16_XP(y1, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(y3, Y1_va, Y0, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(y4, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y6, Y1_va, Y0, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(y5, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(y7, Y1_va, Y0, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(y8, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y10, Y1_va, Y0, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(y9, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(y11, Y1_va, Y0, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(y12, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2);
    BBE_SALIGNVRNX16_XP(y14, Y1_va, Y0, +1 * 4 * BBE_SIMD_WIDTH / 2);

    BBE_SVINTLARNX16_XP(y13, Y0_va, Y0, +2 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SVINTLARNX16_XP(y15, Y1_va, Y0, -1 * 4 * BBE_SIMD_WIDTH / 2, 0);
    BBE_SALIGNVRNX16_XP(y13, Y0_va, Y0, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 13 * 4 * BBE_SIMD_WIDTH / 2);
  }
  //********************************************************************
  __Pragma("no_reorder");
  Y0 = (xb_vecNx16*)x;
  X0 = (xb_vecNx16*)(y);
  X1 = (xb_vecNx16*)(y + 2 * (BBE_SIMD_WIDTH + 4));
  X2 = (xb_vecNx16*)(y + 7 * BBE_SIMD_WIDTH);
  X3 = (xb_vecNx16*)(y + 12 * BBE_SIMD_WIDTH);
  
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  for (l = 0; l<L; l++)
  {
    tw3 = BBE_MOVVA16C(0x91260000);
    x1_align = BBE_LANX16_PP(X1);
    BBE_LVNX16_IP(y0, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LANX16_IP(y1, x1_align, X1);

    X1 = (xb_vecNx16*)XT_ADD((uintptr_t)X1, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);
    //x1_align = BBE_LANX16_PP(X1); //260:

    BBE_LVNX16_IP(x0, X2, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, X2, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(x2, X3, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x3, X3, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);

    y3 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_8);
    y5 = BBE_SELNX16I(x3, x2, BBE_SELI_ROTATE_RIGHT_8);

    x0 = BBE_FFTADDSSRNX16(y0, y3);
    x1 = BBE_FFTSUBSSRNX16(y0, y3);
    x2 = BBE_FFTADDSSRNX16(y4, y1);
    x3 = BBE_FFTSUBSSRNX16(y4, y1);
    x4 = BBE_FFTADDSSRNX16(y2, y5);
    x5 = BBE_FFTSUBSSRNX16(y2, y5);

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

    // Load twiddles and constants into the registers file
    tw2_0 = BBE_LVNX16_I(C1, 0 * 2 * BBE_SIMD_WIDTH);
    tw2_1 = BBE_LVNX16_I(C1, 1 * 2 * BBE_SIMD_WIDTH);

    tw30_5_0 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw30_5_1 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw30_5_2 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_2);
    tw30_5_3 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_3);
    tw30_5_4 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_0);

    x0 = y0;
    x1 = BBE_MULNX16JPACKQ(y5, tw30_5_0);
    x2 = BBE_MULNX16JPACKQ(y1, tw30_5_1);
    x3 = BBE_MULNX16JPACKQ(y3, tw30_5_2);
    x4 = BBE_MULNX16JPACKQ(y2, tw30_5_3);
    x5 = BBE_MULNX16JPACKQ(y4, tw30_5_4);
    Y0_va = BBE_MOVUVR(x0);
    Y1_va = BBE_MOVUVR(x2);

    BBE_SVINTLARNX16_XP(x1, Y0_va, Y0, +1 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SVINTLARNX16_XP(x3, Y1_va, Y0, +2 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SALIGNVRNX16_XP(x4, Y0_va, Y0, +1 * 2 * BBE_SIMD_WIDTH);
    BBE_SALIGNVRNX16_XP(x4, Y1_va, Y0, -2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(x5, Y0_va, Y0, +3 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SALIGNVRNX16_XP(x5, Y0_va, Y0, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 5 * 2 * BBE_SIMD_WIDTH);
  }

  Y0 = (xb_vecNx16*)(x + 6 * BBE_SIMD_WIDTH);
  X0 = (xb_vecNx16*)(y + BBE_SIMD_WIDTH);
  X1 = (xb_vecNx16*)(y + 3 * BBE_SIMD_WIDTH + 8);
  X2 = (xb_vecNx16*)(y + 8 * BBE_SIMD_WIDTH);
  X3 = (xb_vecNx16*)(y + 13 * BBE_SIMD_WIDTH);
  
  for (l = 0; l<L; l++)
  {
    tw3 = BBE_MOVVA16C(0x91260000);
    x1_align = BBE_LANX16_PP(X1);
    BBE_LVNX16_IP(y0, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LANX16_IP(y1, x1_align, X1);
    BBE_LVNX16_IP(x0, X2, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, X2, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X3, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x3, X3, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);

    y3 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_8);
    y5 = BBE_SELNX16I(x3, x2, BBE_SELI_ROTATE_RIGHT_8);

    X1 = (xb_vecNx16*)XT_ADD((uintptr_t)X1, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);
   // x1_align = BBE_LANX16_PP(X1);  //352: 

    x0 = BBE_FFTADDSSRNX16(y0, y3);
    x1 = BBE_FFTSUBSSRNX16(y0, y3);
    x2 = BBE_FFTADDSSRNX16(y4, y1);
    x3 = BBE_FFTSUBSSRNX16(y4, y1);
    x4 = BBE_FFTADDSSRNX16(y2, y5);
    x5 = BBE_FFTSUBSSRNX16(y2, y5);

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

    // Load twiddles and constants into the registers file
    tw2_2 = BBE_LVNX16_I(C1, 2 * 2 * BBE_SIMD_WIDTH);
    tw2_3 = BBE_LVNX16_I(C1, 3 * 2 * BBE_SIMD_WIDTH);

    tw30_5_0 = BBE_SHFLNX16I(tw2_2, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw30_5_1 = BBE_SHFLNX16I(tw2_2, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw30_5_2 = BBE_SHFLNX16I(tw2_2, BBE_SHFLI_REP_2X4_OFFSET_2);
    tw30_5_3 = BBE_SHFLNX16I(tw2_2, BBE_SHFLI_REP_2X4_OFFSET_3);
    tw30_5_4 = BBE_SHFLNX16I(tw2_3, BBE_SHFLI_REP_2X4_OFFSET_0);

    x0 = y0;
    x1 = BBE_MULNX16JPACKQ(y5, tw30_5_0);
    x2 = BBE_MULNX16JPACKQ(y1, tw30_5_1);
    x3 = BBE_MULNX16JPACKQ(y3, tw30_5_2);
    x4 = BBE_MULNX16JPACKQ(y2, tw30_5_3);
    x5 = BBE_MULNX16JPACKQ(y4, tw30_5_4);

    Y0_va = BBE_MOVUVR(x0);
    Y1_va = BBE_MOVUVR(x2);

    BBE_SVINTLARNX16_XP(x1, Y0_va, Y0, +1 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SVINTLARNX16_XP(x3, Y1_va, Y0, +2 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SALIGNVRNX16_XP(x4, Y0_va, Y0, +1 * 2 * BBE_SIMD_WIDTH);
    BBE_SALIGNVRNX16_XP(x4, Y1_va, Y0, -2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(x5, Y0_va, Y0, +3 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SALIGNVRNX16_XP(x5, Y0_va, Y0, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 5 * 2 * BBE_SIMD_WIDTH);
  }

  Y0 = (xb_vecNx16*)(x + 12 * BBE_SIMD_WIDTH);
  X0 = (xb_vecNx16*)(y + 2 * BBE_SIMD_WIDTH);
  X1 = (xb_vecNx16*)(y + 4 * BBE_SIMD_WIDTH + 8);
  X2 = (xb_vecNx16*)(y + 9 * BBE_SIMD_WIDTH + 8);
  X3 = (xb_vecNx16*)(y + 14 * BBE_SIMD_WIDTH + 8);
  x1_align = BBE_LANX16_PP(X1);
  x2_align = BBE_LANX16_PP(X2);
  //x3_align = BBE_LANX16_PP(X3);
  for (l = 0; l<L; l++)
  {
    x1_align = BBE_LANX16_PP(X1);
    x2_align = BBE_LANX16_PP(X2);
    x3_align = BBE_LAVNX16_PP(X3);
    tw3 = BBE_MOVVA16C(0x91260000);

    BBE_LVNX16_IP(y0, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X0, 5 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 5 * 2 * BBE_SIMD_WIDTH);

    BBE_LANX16_IP(y1, x1_align, X1);
    BBE_LANX16_IP(y3, x2_align, X2);
    BBE_LAVNX16_XP(y5, x3_align, X3, BBE_SIMD_WIDTH);   //433: 

    X1 = (xb_vecNx16*)XT_ADD((uintptr_t)X1, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);
    X2 = (xb_vecNx16*)XT_ADD((uintptr_t)X2, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 1 * 2 * BBE_SIMD_WIDTH);
    X3 = (xb_vecNx16*)XT_ADD((uintptr_t)X3, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - BBE_SIMD_WIDTH);
    x0 = BBE_FFTADDSSRNX16(y0, y3);
    x1 = BBE_FFTSUBSSRNX16(y0, y3);
    x2 = BBE_FFTADDSSRNX16(y4, y1);
    x3 = BBE_FFTSUBSSRNX16(y4, y1);
    x4 = BBE_FFTADDSSRNX16(y2, y5);
    x5 = BBE_FFTSUBSSRNX16(y2, y5);

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

    // Load twiddles and constants into the registers file
    tw2_0 = BBE_LVNX16_I(C1, 4 * 2 * BBE_SIMD_WIDTH);
    tw2_1 = BBE_LVNX16_I(C1, 5 * 2 * BBE_SIMD_WIDTH);
    tw30_5_0 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw30_5_1 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw30_5_2 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_2);
    tw30_5_3 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_3);
    tw30_5_4 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_0);

    x1 = BBE_MULNX16JPACKQ(y5, tw30_5_0);
    x2 = BBE_MULNX16JPACKQ(y1, tw30_5_1);
    x3 = BBE_MULNX16JPACKQ(y3, tw30_5_2);
    x4 = BBE_MULNX16JPACKQ(y2, tw30_5_3);
    x5 = BBE_MULNX16JPACKQ(y4, tw30_5_4);

    Y0_va = BBE_MOVUVR(y0);
    Y1_va = BBE_MOVUVR(x2);
    x5 = BBE_SELNX16I(x5, x4, BBE_SELI_EXTRACT_LO_HALVES);
    BBE_SVINTLARNX16_XP(x1, Y0_va, Y0, +1 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SVINTLARNX16_XP(x3, Y1_va, Y0, +1 * 2 * BBE_SIMD_WIDTH, 1);
    BBE_SVNX16_IP(x5, Y0, BCNFFT_BUF_SIZE(120)*sizeof(complex_fract16) - 2 * 2 * BBE_SIMD_WIDTH);
  }
  //********************************************************************
  __Pragma("no_reorder");
  X0 = (xb_vecNx16*)x;
  X1 = (xb_vecNx16*)(x + BBE_SIMD_WIDTH);
  X2 = (xb_vecNx16*)(x + 2 * BBE_SIMD_WIDTH);
  X3 = (xb_vecNx16*)(x + 9 * BBE_SIMD_WIDTH);
  Y0 = (xb_vecNx16*)(y);
  Y1 = (xb_vecNx16*)(y + BBE_SIMD_WIDTH);
  Y2 = (xb_vecNx16*)(y + 2 * BBE_SIMD_WIDTH);
  BBE_FFTWMODE(0x10);
  for (l = 0; l<L; l++)
  {

    BBE_LVNX16_IP(x0, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVB_IP(X0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVA_IP(X0, 3 * 2 * BBE_SIMD_WIDTH);

    y0 = BBE_FFTSRANX16(x0);
    y3 = BBE_FFTADD4SABNX16(x1, x2, 0, 1);
    y4 = BBE_FFTADD4SABNX16(x1, x2, 1, 1);
    y5 = BBE_FFTADD4SABNX16(x1, x2, 2, 1);
    y6 = BBE_FFTADD4SABNX16(x1, x2, 3, 1);

    y1 = BBE_FFTADDSSRNX16(y0, y3);

    tw5_0 = BBE_MOVVA16C(0x3cdeda62);
    tw5_1 = BBE_MOVVA16C(0x0000478e);
    tw5_2 = BBE_MOVVA16C(0x3cde259e);

    {
      xb_vecNx16 tmpc;
      tmpc = BBE_SRAINX16(y3, 2);
      y2 = BBE_SUBSNX16(y0, tmpc);   /* t2 = t0 - (t3 >> 2);*/
    }
    y4 = BBE_MULNX16CPACKQ(y4, tw5_0);
    y5 = BBE_MULNX16CPACKQ(y5, tw5_1);
    y6 = BBE_MULNX16CPACKQ(y6, tw5_2);

    BBE_MOVSCV(y5);
    BBE_MOVSDV(y6);

    x0 = BBE_FFTADD4SCDNX16(y2, y4, 0, 0);
    x1 = BBE_FFTADD4SCDNX16(y2, y4, 1, 0);
    x2 = BBE_FFTADD4SCDNX16(y2, y4, 2, 0);
    x3 = BBE_FFTADD4SCDNX16(y2, y4, 3, 0);

    BBE_SVNX16_IP(y1, Y0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x0, Y0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, Y0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, Y0, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, Y0, 3 * 2 * BBE_SIMD_WIDTH);


    BBE_LVNX16_IP(x5, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x6, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x7, X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVB_IP(X1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVA_IP(X1, 3 * 2 * BBE_SIMD_WIDTH);

    y0 = BBE_FFTSRANX16(x5);
    y3 = BBE_FFTADD4SABNX16(x6, x7, 0, 1);
    y4 = BBE_FFTADD4SABNX16(x6, x7, 1, 1);
    y5 = BBE_FFTADD4SABNX16(x6, x7, 2, 1);
    y6 = BBE_FFTADD4SABNX16(x6, x7, 3, 1);

    y1 = BBE_FFTADDSSRNX16(y0, y3);

    {
      xb_vecNx16 tmpc;
      tmpc = BBE_SRAINX16(y3, 2);
      y2 = BBE_SUBSNX16(y0, tmpc);   /* t2 = t0 - (t3 >> 2);*/
    }

    tw5_0 = BBE_MOVVA16C(0x3cdeda62);
    tw5_1 = BBE_MOVVA16C(0x0000478e);
    tw5_2 = BBE_MOVVA16C(0x3cde259e);

    y4 = BBE_MULNX16CPACKQ(y4, tw5_0);
    y5 = BBE_MULNX16CPACKQ(y5, tw5_1);
    y6 = BBE_MULNX16CPACKQ(y6, tw5_2);

    BBE_MOVSCV(y5);
    BBE_MOVSDV(y6);

    x0 = BBE_FFTADD4SCDNX16(y2, y4, 0, 0);
    x1 = BBE_FFTADD4SCDNX16(y2, y4, 1, 0);
    x2 = BBE_FFTADD4SCDNX16(y2, y4, 2, 0);
    x3 = BBE_FFTADD4SCDNX16(y2, y4, 3, 0);

    BBE_SVNX16_IP(y1, Y1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x0, Y1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, Y1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, Y1, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, Y1, 3 * 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(x10, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x11, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x12, X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVB_IP(X2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_LVA_IP(X2, 3 * 2 * BBE_SIMD_WIDTH);

    y0 = BBE_FFTSRANX16(x10);
    y3 = BBE_FFTADD4SABNX16(x11, x12, 0, 1);
    y4 = BBE_FFTADD4SABNX16(x11, x12, 1, 1);
    y5 = BBE_FFTADD4SABNX16(x11, x12, 2, 1);
    y6 = BBE_FFTADD4SABNX16(x11, x12, 3, 1);

    y1 = BBE_FFTADDSSRNX16(y0, y3);

    {
      xb_vecNx16 tmpc;
      tmpc = BBE_SRAINX16(y3, 2);
      y2 = BBE_SUBSNX16(y0, tmpc);   /* t2 = t0 - (t3 >> 2);*/
    }
    y4 = BBE_MULNX16CPACKQ(y4, tw5_0);
    y5 = BBE_MULNX16CPACKQ(y5, tw5_1);
    y6 = BBE_MULNX16CPACKQ(y6, tw5_2);

    BBE_MOVSCV(y5);
    BBE_MOVSDV(y6);

    x0 = BBE_FFTADD4SCDNX16(y2, y4, 0, 0);
    x1 = BBE_FFTADD4SCDNX16(y2, y4, 1, 0);
    x2 = BBE_FFTADD4SCDNX16(y2, y4, 2, 0);
    x3 = BBE_FFTADD4SCDNX16(y2, y4, 3, 0);

    BBE_SVNX16_IP(y1, Y2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x0, Y2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, Y2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, Y2, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, Y2, 3 * 2 * BBE_SIMD_WIDTH);
  }
} /* bcinfft120() */
#endif

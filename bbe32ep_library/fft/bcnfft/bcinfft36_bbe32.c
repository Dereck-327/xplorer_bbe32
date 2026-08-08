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
DISCARD_FUN(void, bcinfft36,( void * restrict pScr,
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

void bcinfft36 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  int l;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 *          X0;

  const xb_vecNx16 *          C0 = (xb_vecNx16 *)fft36_tw1;
  const xb_vecNx16 *          C1 = (xb_vecNx16 *)fft36_tw2;
  xb_vecNx16 x0, x1, x2, x3, y0, y1, y2, y3;
  xb_vecNx16 x4, x5, x6, x7, y4, y5, y6, y7;
  xb_vecNx16 tw36_9_0, tw36_9_1, tw36_9_2, tw36_9_3, tw36_9_4, tw36_9_5;
  xb_vecNx16 tw9_3, tw9_3_0, tw9_3_1, tw9_3_2, tw9_3_3, tw3;
  valign x_align;
  valign uu0, uu1;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  x_align = BBE_LAVNX16_PP(X0);
  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file
  BBE_LVNX16_IP(tw9_3, C1, 2 * BBE_SIMD_WIDTH);
  tw9_3_0 = BBE_SHFLNX16I(tw9_3, BBE_SHFLI_REP_2X4_OFFSET_0);
  tw9_3_1 = BBE_SHFLNX16I(tw9_3, BBE_SHFLI_REP_2X4_OFFSET_1);
  tw9_3_2 = BBE_SHFLNX16I(tw9_3, BBE_SHFLI_REP_2X4_OFFSET_2);
  tw9_3_3 = BBE_SHFLNX16I(tw9_3, BBE_SHFLI_REP_2X4_OFFSET_3);


  //DFT36 = (DFT3 x I12)*(L9_3 x I4)*(T9_3*(DFT3 x I3) x I4)*L36_9*T36_9*(DFT4 x I9)
  BBE_FFTWMODE(0x00);
  BBE_LVNX16_IP(y0, X0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y1, X0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y2, X0, 2 * BBE_SIMD_WIDTH);
  for (l = 0; l<L - 1; l++)
  {
    // Load twiddles and constants into the registers file
    BBE_LVNX16_IP(tw36_9_0, C0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw36_9_1, C0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw36_9_2, C0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw36_9_3, C0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(tw36_9_4, C0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(tw36_9_5, C0, -5 * 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(y3, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, X0, 2 * BBE_SIMD_WIDTH);

    x0 = y0;
    x1 = BBE_SELNX16I(y2, y1, BBE_SELI_ROTATE_RIGHT_2);
    x2 = BBE_SELNX16I(y3, y2, BBE_SELI_ROTATE_RIGHT_4);
    x3 = BBE_SELNX16I(y4, y3, BBE_SELI_ROTATE_RIGHT_6);

    x4 = y1;
    x5 = BBE_SELNX16I(y3, y2, BBE_SELI_ROTATE_RIGHT_2);
    x6 = BBE_SELNX16I(y4, y3, BBE_SELI_ROTATE_RIGHT_4);
    x7 = BBE_SELNX16I(y4, y4, BBE_SELI_ROTATE_RIGHT_6);

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

    y1 = BBE_MULNX16JPACKQ(y1, tw36_9_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw36_9_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw36_9_2);

    y5 = BBE_MULNX16JPACKQ(y5, tw36_9_3);
    y6 = BBE_MULNX16JPACKQ(y6, tw36_9_4);
    y7 = BBE_MULNX16JPACKQ(y7, tw36_9_5);

    x0 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_LO);
    x2 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_HI);
    x1 = BBE_SELNX16I(y3, y1, BBE_SELI_INTERLEAVE_2_LO);
    x3 = BBE_SELNX16I(y3, y1, BBE_SELI_INTERLEAVE_2_HI);

    x4 = BBE_SELNX16I(y6, y4, BBE_SELI_INTERLEAVE_2_LO);
    x6 = BBE_SELNX16I(y6, y4, BBE_SELI_INTERLEAVE_2_HI);
    x5 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_LO);
    x7 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_HI);

    uu0 = BBE_MOVUVR(x0);
    uu1 = BBE_MOVUVR(x2);
    BBE_SVINTLARNX16_XP(x1, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);
    BBE_SALIGNVRNX16_XP(x4, uu0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(x3, uu1, Y, 2 * BBE_SIMD_WIDTH, 0);
    BBE_SALIGNVRNX16_XP(x5, uu1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(x5, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);

    BBE_LVNX16_IP(y0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X0, 2 * BBE_SIMD_WIDTH);
  }
  // Load twiddles and constants into the registers file
  BBE_LVNX16_IP(tw36_9_0, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw36_9_1, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw36_9_2, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw36_9_3, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw36_9_4, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_XP(tw36_9_5, C0, -5 * 2 * BBE_SIMD_WIDTH);

  BBE_LVNX16_IP(y3, X0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y4, X0, 2 * BBE_SIMD_WIDTH);

  x0 = y0;
  x1 = BBE_SELNX16I(y2, y1, BBE_SELI_ROTATE_RIGHT_2);
  x2 = BBE_SELNX16I(y3, y2, BBE_SELI_ROTATE_RIGHT_4);
  x3 = BBE_SELNX16I(y4, y3, BBE_SELI_ROTATE_RIGHT_6);

  x4 = y1;
  x5 = BBE_SELNX16I(y3, y2, BBE_SELI_ROTATE_RIGHT_2);
  x6 = BBE_SELNX16I(y4, y3, BBE_SELI_ROTATE_RIGHT_4);
  x7 = BBE_SELNX16I(y4, y4, BBE_SELI_ROTATE_RIGHT_6);

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

  y1 = BBE_MULNX16JPACKQ(y1, tw36_9_0);
  y2 = BBE_MULNX16JPACKQ(y2, tw36_9_1);
  y3 = BBE_MULNX16JPACKQ(y3, tw36_9_2);

  y5 = BBE_MULNX16JPACKQ(y5, tw36_9_3);
  y6 = BBE_MULNX16JPACKQ(y6, tw36_9_4);
  y7 = BBE_MULNX16JPACKQ(y7, tw36_9_5);

  x0 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_LO);
  x2 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_HI);
  x1 = BBE_SELNX16I(y3, y1, BBE_SELI_INTERLEAVE_2_LO);
  x3 = BBE_SELNX16I(y3, y1, BBE_SELI_INTERLEAVE_2_HI);

  x4 = BBE_SELNX16I(y6, y4, BBE_SELI_INTERLEAVE_2_LO);
  x6 = BBE_SELNX16I(y6, y4, BBE_SELI_INTERLEAVE_2_HI);
  x5 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_LO);
  x7 = BBE_SELNX16I(y7, y5, BBE_SELI_INTERLEAVE_2_HI);

  uu0 = BBE_MOVUVR(x0);
  uu1 = BBE_MOVUVR(x2);
  BBE_SVINTLARNX16_XP(x1, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);
  BBE_SALIGNVRNX16_XP(x4, uu0, Y, 2 * BBE_SIMD_WIDTH);
  BBE_SVINTLARNX16_XP(x3, uu1, Y, 2 * BBE_SIMD_WIDTH, 0);
  BBE_SALIGNVRNX16_XP(x5, uu1, Y, 2 * BBE_SIMD_WIDTH);
  BBE_SVINTLARNX16_XP(x5, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);
  //********************************************************************
  __Pragma("no_reorder");
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  tw3 = BBE_MOVVA16C(0x91260000);
  BBE_LVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(x3, Y, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(x4, Y, 2 * BBE_SIMD_WIDTH);
  for (l = 0; l<L - 1; l++)
  {
    BBE_LVNX16_IP(x2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, Y, 2 * BBE_SIMD_WIDTH);

    x1 = BBE_SELNX16I(x4, x3, BBE_SELI_ROTATE_RIGHT_8);
    x4 = BBE_SELNX16I(x5, x4, BBE_SELI_ROTATE_RIGHT_8);

    BBE_MOVSAV(x0);
    BBE_MOVSCV(x0);

    /* t0 = t1 + t2 + C + D; D==0; */
    BBE_FFTAVGNX16SB(x1, x2);
    y0 = BBE_FFTADD4SCDNX16(x1, x2, 0, 0);

    x1 = BBE_MULNX16JPACKQ(x1, tw3);
    x2 = BBE_MULNX16JPACKQ(x2, tw3);
    y1 = BBE_FFTADD4SABNX16(x1, x2, 2, 0);
    y2 = BBE_FFTADD4SABNX16(x2, x1, 2, 0);

    BBE_MOVSAV(x3);
    BBE_MOVSCV(x3);

    BBE_FFTAVGNX16SB(x4, x5);
    y3 = BBE_FFTADD4SCDNX16(x4, x5, 0, 0);

    x4 = BBE_MULNX16JPACKQ(x4, tw3);
    x5 = BBE_MULNX16JPACKQ(x5, tw3);
    y4 = BBE_FFTADD4SABNX16(x4, x5, 2, 0);
    y5 = BBE_FFTADD4SABNX16(x5, x4, 2, 0);

    y1 = BBE_MULNX16JPACKQ(y1, tw9_3_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw9_3_1);

    y4 = BBE_MULNX16JPACKQ(y4, tw9_3_2);
    y5 = BBE_MULNX16JPACKQ(y5, tw9_3_3);

    uu0 = BBE_MOVUVR(y0);
    uu1 = BBE_MOVUVR(y2);
    BBE_SVINTLARNX16_XP(y1, uu0, X0, 2 * BBE_SIMD_WIDTH, 1);
    BBE_SVINTLARNX16_XP(y3, uu1, X0, 2 * BBE_SIMD_WIDTH, 1);

    BBE_SALIGNVRNX16_XP(y5, uu0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_SALIGNVRNX16_XP(y4, uu1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(y4, uu0, X0, 2 * BBE_SIMD_WIDTH, 1);

    BBE_LVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, Y, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, Y, 2 * BBE_SIMD_WIDTH);
  }
  BBE_LVNX16_IP(x2, Y, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(x5, Y, 2 * BBE_SIMD_WIDTH);

  x1 = BBE_SELNX16I(x4, x3, BBE_SELI_ROTATE_RIGHT_8);
  x4 = BBE_SELNX16I(x5, x4, BBE_SELI_ROTATE_RIGHT_8);

  BBE_MOVSAV(x0);
  BBE_MOVSCV(x0);

  /* t0 = t1 + t2 + C + D; D==0; */
  BBE_FFTAVGNX16SB(x1, x2);
  y0 = BBE_FFTADD4SCDNX16(x1, x2, 0, 0);

  x1 = BBE_MULNX16JPACKQ(x1, tw3);
  x2 = BBE_MULNX16JPACKQ(x2, tw3);
  y1 = BBE_FFTADD4SABNX16(x1, x2, 2, 0);
  y2 = BBE_FFTADD4SABNX16(x2, x1, 2, 0);

  BBE_MOVSAV(x3);
  BBE_MOVSCV(x3);

  BBE_FFTAVGNX16SB(x4, x5);
  y3 = BBE_FFTADD4SCDNX16(x4, x5, 0, 0);

  x4 = BBE_MULNX16JPACKQ(x4, tw3);
  x5 = BBE_MULNX16JPACKQ(x5, tw3);
  y4 = BBE_FFTADD4SABNX16(x4, x5, 2, 0);
  y5 = BBE_FFTADD4SABNX16(x5, x4, 2, 0);

  y1 = BBE_MULNX16JPACKQ(y1, tw9_3_0);
  y2 = BBE_MULNX16JPACKQ(y2, tw9_3_1);

  y4 = BBE_MULNX16JPACKQ(y4, tw9_3_2);
  y5 = BBE_MULNX16JPACKQ(y5, tw9_3_3);

  uu0 = BBE_MOVUVR(y0);
  uu1 = BBE_MOVUVR(y2);
  BBE_SVINTLARNX16_XP(y1, uu0, X0, 2 * BBE_SIMD_WIDTH, 1);
  BBE_SVINTLARNX16_XP(y3, uu1, X0, 2 * BBE_SIMD_WIDTH, 1);

  BBE_SALIGNVRNX16_XP(y5, uu0, X0, 2 * BBE_SIMD_WIDTH);
  BBE_SALIGNVRNX16_XP(y4, uu1, X0, 2 * BBE_SIMD_WIDTH);
  BBE_SVINTLARNX16_XP(y4, uu0, X0, 2 * BBE_SIMD_WIDTH, 1);
  //********************************************************************
  __Pragma("no_reorder");
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  x_align = BBE_ZALIGN();
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, X0, 2 * BBE_SIMD_WIDTH);
    x2 = BBE_SELNX16I(x5, x3, BBE_SELI_EXTRACT_HI_HALVES);

    BBE_MOVSAV(x0);
    BBE_MOVSCV(x0);

    /* t0 = t1 + t2 + C + D; D==0; */
    BBE_FFTAVGNX16SB(x1, x2);
    y0 = BBE_FFTADD4SCDNX16(x1, x2, 0, 0);

    x1 = BBE_MULNX16JPACKQ(x1, tw3);
    x2 = BBE_MULNX16JPACKQ(x2, tw3);
    y1 = BBE_FFTADD4SABNX16(x1, x2, 2, 0);
    y2 = BBE_FFTADD4SABNX16(x2, x1, 2, 0);

    BBE_MOVSAV(x3);
    BBE_MOVSCV(x3);

    BBE_FFTAVGNX16SB(x4, x5);
    y3 = BBE_FFTADD4SCDNX16(x4, x5, 0, 0);

    x4 = BBE_MULNX16JPACKQ(x4, tw3);
    x5 = BBE_MULNX16JPACKQ(x5, tw3);
    y4 = BBE_FFTADD4SABNX16(x4, x5, 2, 0);
    y5 = BBE_FFTADD4SABNX16(x5, x4, 2, 0);

    BBE_SAVNX16_XP(y0, x_align, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SAVNX16_XP(y3, x_align, Y, BBE_SIMD_WIDTH);
    BBE_SAVNX16_XP(y1, x_align, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SAVNX16_XP(y4, x_align, Y, BBE_SIMD_WIDTH);
    BBE_SAVNX16_XP(y2, x_align, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SAVNX16_XP(y5, x_align, Y, 2 * BBE_SIMD_WIDTH);
  }
  BBE_SAVNX16POS_FP(x_align, Y);
} /* bcinfft36() */
#endif

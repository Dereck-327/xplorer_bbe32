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
#include "fft_tw.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#if !(HAVE_FFT && 1)
DISCARD_FUN(void, bcnfft24,( void * restrict pScr,
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

void bcnfft24 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  int l;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 *          X0;
  xb_vecNx16 *          X1;
  const xb_vecNx16 *          C0 = (xb_vecNx16 *)fft24_tw1;
  const xb_vecNx16 *          C1 = (xb_vecNx16 *)fft24_tw2;
  xb_vecNx16 x0, x1, x2, x3, y0, y1, y2, y3, z0;
  xb_vecNx16 tw24_6_0, tw24_6_1, tw24_6_2, tw6_3, tw6_3_0, tw6_3_1, tw3;
  valign c_align;
  valign uu0;
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  X1 = (xb_vecNx16*)(x + BBE_SIMD_WIDTH);
  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file
  BBE_LVNX16_IP(tw24_6_0, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw24_6_1, C0, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(tw24_6_2, C0, 2 * BBE_SIMD_WIDTH);

  c_align = BBE_LAVNX16_PP(C1);
  BBE_LAVNX16_XP(tw6_3, c_align, C1, 2 * 8);
  tw6_3_0 = BBE_SHFLNX16I(tw6_3, BBE_SHFLI_REP_2X4_OFFSET_0);
  tw6_3_1 = BBE_SHFLNX16I(tw6_3, BBE_SHFLI_REP_2X4_OFFSET_1);
  tw3 = BBE_MOVVA16C(0x91260000);
  // DFT24 = (DFT3 x I8)*(L6_3 x I4)*(T6_3*(DFT2 x I3) x I4)*L24_6*T24_6*(DFT4 x I6)
  BBE_FFTWMODE(0x10);
  z0 = 0;

  for (l = 0; l<L; l++)
  {
    // first stage : L24_6*T24_6*(DFT4 x I6)

    BBE_LVNX16_IP(y0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, X0, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, X1, 3 * 2 * BBE_SIMD_WIDTH);

    x0 = y0;
    x1 = BBE_SELNX16I(y1, y0, BBE_SELI_ROTATE_LEFT_4);
    x2 = BBE_SELNX16I(y2, y1, BBE_SELI_ROTATE_RIGHT_8);
    x3 = BBE_SELNX16I(z0, y2, BBE_SELI_ROTATE_RIGHT_4);
    BBE_MOVSAV(x2);
    BBE_MOVSBV(x3);

    //
    // DFT4 x I6
    //

    y0 = BBE_FFTADD4SABNX16(x0, x1, 0, 0);
    y1 = BBE_FFTADD4SABNX16(x0, x1, 1, 0);
    y2 = BBE_FFTADD4SABNX16(x0, x1, 2, 0);
    y3 = BBE_FFTADD4SABNX16(x0, x1, 3, 0);

    //
    // T24_6
    //

    y1 = BBE_MULNX16CPACKQ(y1, tw24_6_0);
    y2 = BBE_MULNX16CPACKQ(y2, tw24_6_1);
    y3 = BBE_MULNX16CPACKQ(y3, tw24_6_2);

    BBE_DSELNX16I(x3, x1, y3, y1, BBE_DSELI_INTERLEAVE_2);
    x0 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_LO);
    x2 = BBE_SELNX16I(y2, y0, BBE_SELI_INTERLEAVE_2_HI);

    uu0 = BBE_MOVUVR(x0);
    BBE_SVINTLARNX16_XP(x1, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);
    BBE_SALIGNVRNX16_XP(x2, uu0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVINTLARNX16_XP(x3, uu0, Y, 2 * BBE_SIMD_WIDTH, 0);
  }
  __Pragma("no_reorder");
  //********************************************************************
  // Stage 2: (T6_3*(DFT2 x I3) x I4)
  BBE_MOVSAV(0);
  BBE_MOVSBV(0);
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, Y, 2 * BBE_SIMD_WIDTH);
    x0 = y0;
    x2 = y1;
    x1 = BBE_SELNX16I(y2, y1, BBE_SELI_ROTATE_RIGHT_8);
    x3 = BBE_SHFLNX16I(y2, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

    y0 = BBE_FFTADD4SABNX16(x0, x1, 0, 0);
    y1 = BBE_FFTADD4SABNX16(x0, x1, 2, 0);
    y2 = BBE_FFTADD4SABNX16(x2, x3, 0, 0);
    y3 = BBE_FFTADD4SABNX16(x2, x3, 2, 0);

    y1 = BBE_MULNX16CPACKQ(y1, tw6_3_0);
    y3 = BBE_MULNX16CPACKQ(y3, tw6_3_1);

    x0 = BBE_SELNX16I(y1, y0, BBE_SELI_EXTRACT_LO_HALVES);
    x1 = BBE_SELNX16I(y1, y0, BBE_SELI_EXTRACT_HI_HALVES);
    x2 = BBE_SELNX16I(y3, y2, BBE_SELI_EXTRACT_LO_HALVES);

    BBE_SVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, X0, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  //********************************************************************
  // Stage 3: (DFT3 x I8)
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X0, 2 * BBE_SIMD_WIDTH);
    BBE_MOVSAV(x0);
    BBE_MOVSCV(x0);
    BBE_FFTAVGNX16SB(x1, x2);
    /* t0 = t1 + t2 + C + D; D==0; */
    y0 = BBE_FFTADD4SCDNX16(x1, x2, 0, 0);

    x1 = BBE_MULNX16CPACKQ(x1, tw3);
    x2 = BBE_MULNX16CPACKQ(x2, tw3);
    y1 = BBE_FFTADD4SABNX16(x1, x2, 2, 0);
    y2 = BBE_FFTADD4SABNX16(x2, x1, 2, 0);

    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y2, Y, 2 * BBE_SIMD_WIDTH);
  }
} /* bcnfft24() */
#endif

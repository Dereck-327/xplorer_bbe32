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
DISCARD_FUN(void, bcinfft48,( void * restrict pScr,
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

void bcinfft48 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  int l;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 *          X0;

  const xb_vecNx16 *          C0 = (xb_vecNx16 *)fft48_tw1;
  const xb_vecNx16 *          C1 = (xb_vecNx16 *)fft48_tw2;
  xb_vecNx16 x0, x1, x2, x3, y0, y1, y2, y3, z0;
  xb_vecNx16 x4, x5, x6, x7, y4, y5, y6, y7;
  xb_vecNx16 tw48_12_0, tw48_12_1, tw48_12_2, tw48_12_3, tw48_12_4, tw48_12_5;
  xb_vecNx16 tw2_0, tw2_1, tw12_3_0, tw12_3_1, tw12_3_2, tw12_3_3, tw12_3_4, tw12_3_5, tw3;
  xb_vecNx16 * px = (xb_vecNx16*)x;
  valign vx; 
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  vx = BBE_ZALIGN(); 
  //----------------------------------------------------------------------------

  tw3 = BBE_MOVVA16C(0x91260000);

  //DFT48 = (DFT3 x I16)*(L12_3 x I4)*(T12_3*(DFT4 x I3) x I4)*L48_12*T48_12*(DFT4 x I12)
  BBE_FFTWMODE(0x00);
  z0 = 0;
  for (l = 0; l < L ; l++)
  {

    BBE_LAVNX16_XP(y0, vx, px, 2 * BBE_SIMD_WIDTH);
    BBE_LAVNX16_XP(y1, vx, px,     BBE_SIMD_WIDTH);

    BBE_LAVNX16_XP(x1, vx, px, 2 * BBE_SIMD_WIDTH);
    BBE_LAVNX16_XP(x5, vx, px,     BBE_SIMD_WIDTH);

    BBE_LAVNX16_XP(y3, vx, px, 2 * BBE_SIMD_WIDTH);
    BBE_LAVNX16_XP(y4, vx, px,     BBE_SIMD_WIDTH);

    BBE_LAVNX16_XP(x3, vx, px, 2 * BBE_SIMD_WIDTH);
    BBE_LAVNX16_XP(x7, vx, px,     BBE_SIMD_WIDTH);

    // Load twiddles and constants into the registers file
    tw48_12_0 = BBE_LVNX16_I(C0, 0 * 2 * BBE_SIMD_WIDTH);
    tw48_12_1 = BBE_LVNX16_I(C0, 1 * 2 * BBE_SIMD_WIDTH);
    tw48_12_2 = BBE_LVNX16_I(C0, 2 * 2 * BBE_SIMD_WIDTH);
    tw48_12_3 = BBE_LVNX16_I(C0, 3 * 2 * BBE_SIMD_WIDTH);
    tw48_12_4 = BBE_LVNX16_I(C0, 4 * 2 * BBE_SIMD_WIDTH);
    tw48_12_5 = BBE_LVNX16_I(C0, 5 * 2 * BBE_SIMD_WIDTH);

    x0 = y0;
    x4 = y1;
    x2 = y3;
    x6 = y4;

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

    y1 = BBE_MULNX16JPACKQ(y1, tw48_12_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw48_12_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw48_12_2);

    y5 = BBE_MULNX16JPACKQ(y5, tw48_12_3);
    y6 = BBE_MULNX16JPACKQ(y6, tw48_12_4);
    y7 = BBE_MULNX16JPACKQ(y7, tw48_12_5);

    BBE_DSELNX16I(x1, x0, y2, y0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x3, x2, y3, y1, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(y1, y0, x2, x0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y3, y2, x3, x1, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(x5, x4, y6, y4, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x7, x6, y7, y5, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(y5, y4, x6, x4, BBE_DSELI_INTERLEAVE_2);

    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y3, Y, 2 * BBE_SIMD_WIDTH);

    BBE_SVNX16_IP(y4, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y5, Y, 2 * BBE_SIMD_WIDTH);

  }

  //********************************************************************
  __Pragma("no_reorder");
  px = (xb_vecNx16*)y;
  vx = BBE_ZALIGN(); 
  Y = (xb_vecNx16*)x;
  for (l = 0; l<L; l++)
  {
        BBE_LAVNX16_XP(y0, vx, px, 2 * BBE_SIMD_WIDTH);
        BBE_LAVNX16_XP(y1, vx, px,     BBE_SIMD_WIDTH);

        BBE_LAVNX16_XP(x1, vx, px, 2 * BBE_SIMD_WIDTH);
        BBE_LAVNX16_XP(x5, vx, px, BBE_SIMD_WIDTH);

        BBE_LAVNX16_XP(y3, vx, px, 2 * BBE_SIMD_WIDTH);
        BBE_LAVNX16_XP(y4, vx, px, BBE_SIMD_WIDTH);

        BBE_LAVNX16_XP(x3, vx, px, 2 * BBE_SIMD_WIDTH);
        BBE_LAVNX16_XP(x7, vx, px, BBE_SIMD_WIDTH);

    x0 = y0;
    x4 = y1;
    x2 = y3;
    x6 = y4;

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

    // Load twiddles and constants into the registers file
    BBE_LVNX16_IP(tw2_0, C1, 2 * BBE_SIMD_WIDTH);
    tw12_3_0 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw12_3_1 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw12_3_2 = BBE_SHFLNX16I(tw2_0, BBE_SHFLI_REP_2X4_OFFSET_2);

    BBE_LVNX16_IP(tw2_1, C1, -2 * BBE_SIMD_WIDTH);
    tw12_3_3 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_0);
    tw12_3_4 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_1);
    tw12_3_5 = BBE_SHFLNX16I(tw2_1, BBE_SHFLI_REP_2X4_OFFSET_2);

    y1 = BBE_MULNX16JPACKQ(y1, tw12_3_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw12_3_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw12_3_2);

    y5 = BBE_MULNX16JPACKQ(y5, tw12_3_3);
    y6 = BBE_MULNX16JPACKQ(y6, tw12_3_4);
    y7 = BBE_MULNX16JPACKQ(y7, tw12_3_5);

    x0 = BBE_SELNX16I(y1, y0, BBE_SELI_EXTRACT_LO_HALVES);
    x2 = BBE_SELNX16I(y1, y0, BBE_SELI_EXTRACT_HI_HALVES);
    x1 = BBE_SELNX16I(y3, y2, BBE_SELI_EXTRACT_LO_HALVES);
    x3 = BBE_SELNX16I(y3, y2, BBE_SELI_EXTRACT_HI_HALVES);

    x4 = BBE_SELNX16I(y5, y4, BBE_SELI_EXTRACT_LO_HALVES);
    x5 = BBE_SELNX16I(y7, y6, BBE_SELI_EXTRACT_LO_HALVES);

    BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x4, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x5, Y, 2 * BBE_SIMD_WIDTH);
  }
 
  __Pragma("no_reorder");
  //********************************************************************
  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  for (l = 0; l<L; l++)
  {
    BBE_LVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, X0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, X0, 2 * BBE_SIMD_WIDTH);

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

    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y3, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y4, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y5, Y, 2 * BBE_SIMD_WIDTH);

  }
} /* bcinfft48() */
#endif

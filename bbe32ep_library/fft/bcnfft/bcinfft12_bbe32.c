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
DISCARD_FUN(void, bcinfft12,( void * restrict pScr,
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

void bcinfft12 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L )
{
  int l;
  xb_vecNx16 * restrict Y;
  xb_vecNx16 *          X0;
  const xb_vecNx16 *          C0 = (xb_vecNx16 *)bcnfft12_tw;
  xb_vecNx16 x0, x1, x2, x3, y0, y1, y2, y3;
  xb_vecNx16 tw12_3_0, tw12_3_1, tw12_3_2, tw3;
  valign x_align, y_align;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(pScr);
  NASSERT(L > 0);

  Y = (xb_vecNx16*)y;
  X0 = (xb_vecNx16*)x;
  //----------------------------------------------------------------------------
  // Load twiddles and constants into the registers file
  y_align = BBE_LAVNX16_PP(C0);
  BBE_LAVNX16_XP(tw12_3_0, y_align, C0, 2 * 6);
  BBE_LAVNX16_XP(tw12_3_1, y_align, C0, 2 * 6);
  BBE_LAVNX16_XP(tw12_3_2, y_align, C0, 2 * 6);
  tw3 = BBE_MOVVA16C(0x91260000);
  // DFT12 = (DFT3 x I4)*(L12_3)*(T12_3)*(DFT4 x I3)
  BBE_FFTWMODE(0x00);
  BBE_MOVSBV(0);
  BBE_MOVSDV(0);
  x_align = BBE_LAVNX16_PP(X0);
  for (l = 0; l<L; l++)
  {
    // first stage : (L12_3)*(T12_3)*(DFT4 x I3)
    BBE_LAVNX16_XP(x0, x_align, X0, 2 * 6);
    BBE_LAVNX16_XP(x1, x_align, X0, 2 * 6);
    BBE_LAVNX16_XP(x2, x_align, X0, 2 * 6);
    BBE_LAVNX16_XP(x3, x_align, X0, 2 * 14);

    BBE_MOVSAV(x2);
    BBE_MOVSBV(x3);

    //
    // DFT4 x I3
    //

    y0 = BBE_FFTADD4SABNX16(x0, x1, 0, 0);
    y1 = BBE_FFTADD4SABNX16(x0, x1, 1, 0);
    y2 = BBE_FFTADD4SABNX16(x0, x1, 2, 0);
    y3 = BBE_FFTADD4SABNX16(x0, x1, 3, 0);

    //
    // T12_3
    //

    y1 = BBE_MULNX16JPACKQ(y1, tw12_3_0);
    y2 = BBE_MULNX16JPACKQ(y2, tw12_3_1);
    y3 = BBE_MULNX16JPACKQ(y3, tw12_3_2);

    x0 = BBE_SELNX16I(y1, y0, BBE_SELI_INTERLEAVE_2_LO);
    x1 = BBE_SELNX16I(y3, y2, BBE_SELI_INTERLEAVE_2_LO);
    y0 = BBE_SELNX16I(x1, x0, BBE_SELI_INTERLEAVE_4_LO);
    y1 = BBE_SELNX16I(x1, x0, BBE_SELI_INTERLEAVE_4_ODD);
    y2 = BBE_SELNX16I(x1, x0, BBE_SELI_INTERLEAVE_4_HI);

    // DFT3 x I4
    BBE_MOVSAV(y0);
    BBE_MOVSCV(y0);

    BBE_FFTAVGNX16SB(y1, y2);
    /* t0 = t1 + t2 + C + D; D==0; */
    x0 = BBE_FFTADD4SCDNX16(y1, y2, 0, 0);

    y1 = BBE_MULNX16JPACKQ(y1, tw3);
    y2 = BBE_MULNX16JPACKQ(y2, tw3);
    x1 = BBE_FFTADD4SABNX16(y1, y2, 2, 0);
    x2 = BBE_FFTADD4SABNX16(y2, y1, 2, 0);

    y_align = BBE_ZALIGN();
    BBE_SAVNX16_XP(x0, y_align, Y, 2 * 8);
    BBE_SAVNX16_XP(x1, y_align, Y, 2 * 8);
    BBE_SAVNX16_XP(x2, y_align, Y, 2 * 16);
    BBE_SAVNX16POS_FP(y_align, Y);
  }
} /* bcinfft12() */
#endif

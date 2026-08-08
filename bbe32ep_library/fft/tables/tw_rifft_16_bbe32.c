/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
	NatureDSP_Baseband library. FFT part.
    Radix-2 forward/inverse FFT on real data, auto scaling
    Twiddle factor tables 16-point real-valued inverse DFT
	IntegrIT, 2006-2016
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Twiddle factor table declarations. */
#include "fft_tw.h"
/* Common utility declarations. */
#include "common.h"

const int16_t ALIGN(32) fft_tw_tab_rifft_16[3*16+9*16] = {
  // Store Nyquist frequency bin value to the imaginary part of 0-th bin.
  0,16,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  // a0 sequence select pattern
  0,1,2,3,4,5,6,7,24,25,22,23,20,21,18,19,
  // a1 sequence select pattern
  1,17,30,31,28,29,26,27,24,25,10,11,12,13,14,15,
  // N = 16;
  // twd = 1j*exp(-2*pi*1j*1/2*(0:N/4-1)'/(N/2));
  // tw0 = [1-conj(twd);1;wrev(1+twd(2:N/4))];
  // tw1 = [1+conj(twd);1;wrev(1-twd(2:N/4))];
  // IDFTN = conj(fft(eye(N/2)));
  // tw0, CQ14
  (int16_t)0x4000,(int16_t)0x4000,(int16_t)0x2782,(int16_t)0x3b21,
  (int16_t)0x12bf,(int16_t)0x2d41,(int16_t)0x04df,(int16_t)0x187e,
  (int16_t)0x4000,(int16_t)0x0000,(int16_t)0x7b21,(int16_t)0x187e,
  (int16_t)0x6d41,(int16_t)0x2d41,(int16_t)0x587e,(int16_t)0x3b21,
  // tw1, CQ14
  (int16_t)0x4000,(int16_t)0xc000,(int16_t)0x587e,(int16_t)0xc4df,
  (int16_t)0x6d41,(int16_t)0xd2bf,(int16_t)0x7b21,(int16_t)0xe782,
  (int16_t)0x4000,(int16_t)0x0000,(int16_t)0x04df,(int16_t)0xe782,
  (int16_t)0x12bf,(int16_t)0xd2bf,(int16_t)0x2782,(int16_t)0xc4df,
  // IDFT8(:,2:N), CQ15
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0x5a82,
  (int16_t)0x0000,(int16_t)0x7fff,(int16_t)0xa57e,(int16_t)0x5a82,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0xa57e,(int16_t)0xa57e,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x7fff,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x7fff,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0xa57e,(int16_t)0x5a82,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x5a82,(int16_t)0x5a82,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x0000,(int16_t)0x7fff,(int16_t)0xa57e,(int16_t)0xa57e,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0xa57e,(int16_t)0xa57e,
  (int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0x5a82,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0x5a82,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x7fff,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x7fff,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xa57e,
  (int16_t)0x8000,(int16_t)0x0000,(int16_t)0xa57e,(int16_t)0x5a82,
  (int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x5a82,(int16_t)0x5a82
};

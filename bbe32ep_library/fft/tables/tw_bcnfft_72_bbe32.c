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
    Blockwise mixed radix forward/inverse FFT on complex data, no scaling
    Twiddle factor tables for 72-point transform
	IntegrIT, 2006-2016
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Twiddle factor table declarations. */
#include "fft_tw.h"
/* Common utility declarations. */
#include "common.h"

const int16_t ALIGN(32) bcnfft72_T18_6[] = {
  // Tp18_6 = ( L6_3 x I2 )*T18_6(7:18), CQ15
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7848,(int16_t)0xd439,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x620e,(int16_t)0xadb9,
  (int16_t)0x620e,(int16_t)0xadb9,(int16_t)0x4000,(int16_t)0x9126,
  (int16_t)0x163a,(int16_t)0x81f2,(int16_t)0xc000,(int16_t)0x9126,
  (int16_t)0x163a,(int16_t)0x81f2,(int16_t)0xe9c6,(int16_t)0x81f2,
  (int16_t)0x87b8,(int16_t)0xd439,(int16_t)0x87b8,(int16_t)0x2bc7,
  (int16_t)0x0000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x0000,
  (int16_t)0x0000,(int16_t)0x0000,(int16_t)0x0000,(int16_t)0x0000,
};

const int16_t ALIGN(32) bcnfft72_r3tw[] = {
  // -1j*(3^0.5)/2, CQ15
  (int16_t)0x0000,(int16_t)0x9126,(int16_t)0x0000,(int16_t)0x9126,
  (int16_t)0x0000,(int16_t)0x9126,(int16_t)0x0000,(int16_t)0x9126,
  (int16_t)0x0000,(int16_t)0x9126,(int16_t)0x0000,(int16_t)0x9126,
  (int16_t)0x0000,(int16_t)0x9126,(int16_t)0x0000,(int16_t)0x9126
};

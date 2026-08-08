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
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_fft.h"
#include "common.h"
#include "fft_tw.h"



ALIGN(32) const int16_t fft16_tw[16+2*16+2*16] = {
  // +1 replicated 4 times
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000, 
  // -1j replicated 4 times
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,(int16_t)0x8000, 
  // T16_4
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7642,(int16_t)0xcf04,
  (int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x30fc,(int16_t)0x89be, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xa57e, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x30fc,(int16_t)0x89be,
  (int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0x89be,(int16_t)0x30fc,
  // Select patterns for 4x4 matrix transpose.
  0,1,8,9,16,17,24,25,2,3,10,11,18,19,26,27,
  4,5,12,13,20,21,28,29,6,7,14,15,22,23,30,31
};


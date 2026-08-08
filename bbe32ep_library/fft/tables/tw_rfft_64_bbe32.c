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
    Twiddle factor tables for real-to-complex spectrum conversion routines
	IntegrIT, 2006-2016
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Twiddle factor table declarations. */
#include "fft_tw.h"
/* Common utility declarations. */
#include "common.h"

const int16_t ALIGN(32) fft_tw_tab_rfft_64[6*16+5*16] = {
  //----------------------------------------------------------------------------
  // Complex-valued 32-point FFT twiddle table
  // T32_8, 8..31
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7d8a,(int16_t)0xe707,
  (int16_t)0x7642,(int16_t)0xcf04,(int16_t)0x6a6e,(int16_t)0xb8e3, 
  (int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x471d,(int16_t)0x9592,
  (int16_t)0x30fc,(int16_t)0x89be,(int16_t)0x18f9,(int16_t)0x8276, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7642,(int16_t)0xcf04,
  (int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x30fc,(int16_t)0x89be, 
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xcf04,(int16_t)0x89be,
  (int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0x89be,(int16_t)0xcf04, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x6a6e,(int16_t)0xb8e3,
  (int16_t)0x30fc,(int16_t)0x89be,(int16_t)0xe707,(int16_t)0x8276, 
  (int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0x8276,(int16_t)0xe707,
  (int16_t)0x89be,(int16_t)0x30fc,(int16_t)0xb8e3,(int16_t)0x6a6e, 
  // T8_2 x I4, 8..31
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000, 
  (int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x5a82,(int16_t)0xa57e,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000, 
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,(int16_t)0x8000,
  (int16_t)0x0000,(int16_t)0x8000,(int16_t)0x0000,(int16_t)0x8000, 
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
  (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000, 
  (int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0xa57e,
  (int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0xa57e,
  //----------------------------------------------------------------------------
  // Real-to-complex spectrum converter twiddle table and constants
  // N = 64;
  // twd = reshape(1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2)),8,N/32).';
  // twd = [(1-twd)/2,(1+twd)/2].';
  // twd = reshape([real(twd(:))';imag(twd(:))'],4,N/4)';
  // twd = max(-2^15,min(2^15-1,floor(twd*2^15+0.5)));
  // twd = mod(2^16+twd,2^16);
  (int16_t)0x4000,(int16_t)0xc000,(int16_t)0x39ba,(int16_t)0xc04f,
  (int16_t)0x3384,(int16_t)0xc13b,(int16_t)0x2d6c,(int16_t)0xc2c1,
  (int16_t)0x2782,(int16_t)0xc4df,(int16_t)0x21d5,(int16_t)0xc78f,
  (int16_t)0x1c72,(int16_t)0xcac9,(int16_t)0x1766,(int16_t)0xce87,
  (int16_t)0x4000,(int16_t)0x4000,(int16_t)0x4646,(int16_t)0x3fb1,
  (int16_t)0x4c7c,(int16_t)0x3ec5,(int16_t)0x5294,(int16_t)0x3d3f,
  (int16_t)0x587e,(int16_t)0x3b21,(int16_t)0x5e2b,(int16_t)0x3871,
  (int16_t)0x638e,(int16_t)0x3537,(int16_t)0x689a,(int16_t)0x3179,
  (int16_t)0x12bf,(int16_t)0xd2bf,(int16_t)0x0e87,(int16_t)0xd766,
  (int16_t)0x0ac9,(int16_t)0xdc72,(int16_t)0x078f,(int16_t)0xe1d5,
  (int16_t)0x04df,(int16_t)0xe782,(int16_t)0x02c1,(int16_t)0xed6c,
  (int16_t)0x013b,(int16_t)0xf384,(int16_t)0x004f,(int16_t)0xf9ba,
  (int16_t)0x6d41,(int16_t)0x2d41,(int16_t)0x7179,(int16_t)0x289a,
  (int16_t)0x7537,(int16_t)0x238e,(int16_t)0x7871,(int16_t)0x1e2b,
  (int16_t)0x7b21,(int16_t)0x187e,(int16_t)0x7d3f,(int16_t)0x1294,
  (int16_t)0x7ec5,(int16_t)0x0c7c,(int16_t)0x7fb1,(int16_t)0x0646,
  // Shift-and-reverse pattern for select.
  16,17,14,15,12,13,10,11,8,9,6,7,4,5,2,3
};

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

// ****************** stage 1 ****************** 

ALIGN(32) const int16_t fft96_tw1[] = 
{
	32767,	0,	32698,	-2143,	32488,	-4277,	32138,	-6393,
	31651,	-8481,	31029,	-10533,	30274,	-12540,	29389,	-14493,
	32767,	0,	32488,	-4277,	31651,	-8481,	30274,	-12540,
	28378,	-16384,	25997,	-19948,	23170,	-23170,	19948,	-25997,
	32767,	0,	32138,	-6393,	30274,	-12540,	27246,	-18205,
	23170,	-23170,	18205,	-27246,	12540,	-30274,	6393,	-32138,
	28378,	-16384,	27246,	-18205,	25997,	-19948,	24636,	-21605,
	23170,	-23170,	21605,	-24636,	19948,	-25997,	18205,	-27246,
	16384,	-28378,	12540,	-30274,	8481,	-31651,	4277,	-32488,
	0,	-32768,	-4277,	-32488,	-8481,	-31651,	-12540,	-30274,
	0,	-32768,	-6393,	-32138,	-12540,	-30274,	-18205,	-27246,
	-23170,	-23170,	-27246,	-18205,	-30274,	-12540,	-32138,	-6393,
	16384,	-28378,	14493,	-29389,	12540,	-30274,	10533,	-31029,
	8481,	-31651,	6393,	-32138,	4277,	-32488,	2143,	-32698,
	-16384,	-28378,	-19948,	-25997,	-23170,	-23170,	-25997,	-19948,
	-28378,	-16384,	-30274,	-12540,	-31651,	-8481,	-32488,	-4277,
	-32768,	0,	-32138,	6393,	-30274,	12540,	-27246,	18205,
	-23170,	23170,	-18205,	27246,	-12540,	30274,	-6393,	32138,

};

// ****************** stage 2 ****************** 

ALIGN(32) const int16_t fft96_tw2[] = 
{
	32767,	0,	31651,	-8481,	32767,	0,	28378,	-16384,
	32767,	0,	23170,	-23170,	28378,	-16384,	23170,	-23170,
	16384,	-28378,	0,	-32768,	0,	-32768,	-23170,	-23170,
	16384,	-28378,	8481,	-31651,	-16384,	-28378,	-28378,	-16384,
	-32768,	0,	-23170,	23170,	0,	0,	0,	0,
  0,	0,	0,	0, 0,	0,	0,	0,

};

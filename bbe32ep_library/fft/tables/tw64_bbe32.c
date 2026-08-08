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

ALIGN(32) const int16_t fft64_tw1[] = 
{
	32767,	0,	32610,	-3212,	32138,	-6393,	31357,	-9512,
	30274,	-12540,	28899,	-15447,	27246,	-18205,	25330,	-20788,
	32767,	0,	32138,	-6393,	30274,	-12540,	27246,	-18205,
	23170,	-23170,	18205,	-27246,	12540,	-30274,	6393,	-32138,
	32767,	0,	31357,	-9512,	27246,	-18205,	20788,	-25330,
	12540,	-30274,	3212,	-32610,	-6393,	-32138,	-15447,	-28899,
	23170,	-23170,	20788,	-25330,	18205,	-27246,	15447,	-28899,
	12540,	-30274,	9512,	-31357,	6393,	-32138,	3212,	-32610,
	0,	-32768,	-6393,	-32138,	-12540,	-30274,	-18205,	-27246,
	-23170,	-23170,	-27246,	-18205,	-30274,	-12540,	-32138,	-6393,
	-23170,	-23170,	-28899,	-15447,	-32138,	-6393,	-32610,	3212,
	-30274,	12540,	-25330,	20788,	-18205,	27246,	-9512,	31357,

};

// ****************** stage 2 ****************** 

ALIGN(32) const int16_t fft64_tw2[] = 
{
	32767,	0,	30274,	-12540,	32767,	0,	23170,	-23170,
	32767,	0,	12540,	-30274,	23170,	-23170,	12540,	-30274,
	0,	-32768,	-23170,	-23170,	-23170,	-23170,	-30274,	12540,

};

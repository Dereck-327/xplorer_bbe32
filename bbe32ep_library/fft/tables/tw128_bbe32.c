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

ALIGN(32) const int16_t fft128_tw1[] = 
{
	32767,	0,	32729,	-1608,	32610,	-3212,	32413,	-4808,
	32138,	-6393,	31786,	-7962,	31357,	-9512,	30853,	-11039,
	32767,	0,	32610,	-3212,	32138,	-6393,	31357,	-9512,
	30274,	-12540,	28899,	-15447,	27246,	-18205,	25330,	-20788,
	32767,	0,	32413,	-4808,	31357,	-9512,	29622,	-14010,
	27246,	-18205,	24279,	-22006,	20788,	-25330,	16846,	-28106,
	30274,	-12540,	29622,	-14010,	28899,	-15447,	28106,	-16846,
	27246,	-18205,	26320,	-19520,	25330,	-20788,	24279,	-22006,
	23170,	-23170,	20788,	-25330,	18205,	-27246,	15447,	-28899,
	12540,	-30274,	9512,	-31357,	6393,	-32138,	3212,	-32610,
	12540,	-30274,	7962,	-31786,	3212,	-32610,	-1608,	-32729,
	-6393,	-32138,	-11039,	-30853,	-15447,	-28899,	-19520,	-26320,
	23170,	-23170,	22006,	-24279,	20788,	-25330,	19520,	-26320,
	18205,	-27246,	16846,	-28106,	15447,	-28899,	14010,	-29622,
	0,	-32768,	-3212,	-32610,	-6393,	-32138,	-9512,	-31357,
	-12540,	-30274,	-15447,	-28899,	-18205,	-27246,	-20788,	-25330,
	-23170,	-23170,	-26320,	-19520,	-28899,	-15447,	-30853,	-11039,
	-32138,	-6393,	-32729,	-1608,	-32610,	3212,	-31786,	7962,
	12540,	-30274,	11039,	-30853,	9512,	-31357,	7962,	-31786,
	6393,	-32138,	4808,	-32413,	3212,	-32610,	1608,	-32729,
	-23170,	-23170,	-25330,	-20788,	-27246,	-18205,	-28899,	-15447,
	-30274,	-12540,	-31357,	-9512,	-32138,	-6393,	-32610,	-3212,
	-30274,	12540,	-28106,	16846,	-25330,	20788,	-22006,	24279,
	-18205,	27246,	-14010,	29622,	-9512,	31357,	-4808,	32413,

};

// ****************** stage 2 ****************** 

ALIGN(32) const int16_t fft128_tw2[] = 
{
	32767,	0,	32138,	-6393,	32767,	0,	30274,	-12540,
	32767,	0,	27246,	-18205,	30274,	-12540,	27246,	-18205,
	23170,	-23170,	12540,	-30274,	12540,	-30274,	-6393,	-32138,
	23170,	-23170,	18205,	-27246,	0,	-32768,	-12540,	-30274,
	-23170,	-23170,	-32138,	-6393,	12540,	-30274,	6393,	-32138,
	-23170,	-23170,	-30274,	-12540,	-30274,	12540,	-18205,	27246,

};

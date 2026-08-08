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

ALIGN(64) const int16_t fft60_tw1[] = 
{
	32767,	0,	32588,	-3425,	32052,	-6813,	31164,	-10126,
	29935,	-13328,	28378,	-16384,	26510,	-19261,	24351,	-21926,
	32767,	0,	32052,	-6813,	29935,	-13328,	26510,	-19261,
	21926,	-24351,	16384,	-28378,	10126,	-31164,	3425,	-32588,
	32767,	0,	31164,	-10126,	26510,	-19261,	19261,	-26510,
	10126,	-31164,	0,	-32768,	-10126,	-31164,	-19261,	-26510,
	21926,	-24351,	19261,	-26510,	16384,	-28378,	13328,	-29935,
	10126,	-31164,	6813,	-32052,	3425,	-32588,	0,	0,
	-3425,	-32588,	-10126,	-31164,	-16384,	-28378,	-21926,	-24351,
	-26510,	-19261,	-29935,	-13328,	-32052,	-6813,	0,	0,
	-26510,	-19261,	-31164,	-10126,	-32768,	0,	-31164,	10126,
	-26510,	19261,	-19261,	26510,	-10126,	31164,	0,	0,

};

// ****************** stage 2 ****************** 

ALIGN(64) const int16_t fft60_tw2[] = 
{
	32767,	0,	29935,	-13328,	32767,	0,	21926,	-24351,
	21926,	-24351,	10126,	-31164,	-3425,	-32588,	-26510,	-19261,
	-3425,	-32588,	0,	0,	-32052,	6813,	0,	0,

};

// ****************** stage 3 ****************** 

// ****************** stage 4 ****************** 

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

ALIGN(64) const int16_t fft72_tw1[] = 
{
	32767,	0,	32643,	-2856,	32270,	-5690,	31651,	-8481,
	30792,	-11207,	29698,	-13848,	28378,	-16384,	26842,	-18795,
	32767,	0,	32270,	-5690,	30792,	-11207,	28378,	-16384,
	25102,	-21063,	21063,	-25102,	16384,	-28378,	11207,	-30792,
	32767,	0,	31651,	-8481,	28378,	-16384,	23170,	-23170,
	16384,	-28378,	8481,	-31651,	0,	-32768,	-8481,	-31651,
	25102,	-21063,	23170,	-23170,	21063,	-25102,	18795,	-26842,
	16384,	-28378,	13848,	-29698,	11207,	-30792,	8481,	-31651,
	5690,	-32270,	0,	-32768,	-5690,	-32270,	-11207,	-30792,
	-16384,	-28378,	-21063,	-25102,	-25102,	-21063,	-28378,	-16384,
	-16384,	-28378,	-23170,	-23170,	-28378,	-16384,	-31651,	-8481,
	-32768,	0,	-31651,	8481,	-28378,	16384,	-23170,	23170,
	5690,	-32270,	2856,	-32643,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	-30792,	-11207,	-32270,	-5690,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	-16384,	28378,	-8481,	31651,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,

};

// ****************** stage 2 ****************** 

ALIGN(64) const int16_t fft72_tw2[] = 
{
	32767,	0,	30792,	-11207,	32767,	0,	25102,	-21063,
	32767,	0,	16384,	-28378,	32767,	0,	5690,	-32270,
	32767,	0,	-5690,	-32270,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	25102,	-21063,	0,	0,	5690,	-32270,	0,	0,
	-16384,	-28378,	0,	0,	-30792,	-11207,	0,	0,
	-30792,	11207,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,

};

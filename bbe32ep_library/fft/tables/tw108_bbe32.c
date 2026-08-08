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

ALIGN(64) const int16_t fft108_tw1[] = 
{
	32767,	0,	32713,	-1905,	32546,	-3804,	32270,	-5690,
	31885,	-7557,	31391,	-9398,	30792,	-11207,	30088,	-12979,
	32767,	0,	32546,	-3804,	31885,	-7557,	30792,	-11207,
	29283,	-14706,	27377,	-18006,	25102,	-21063,	22487,	-23835,
	32767,	0,	32270,	-5690,	30792,	-11207,	28378,	-16384,
	25102,	-21063,	21063,	-25102,	16384,	-28378,	11207,	-30792,
	29283,	-14706,	28378,	-16384,	27377,	-18006,	26284,	-19568,
	25102,	-21063,	23835,	-22487,	22487,	-23835,	21063,	-25102,
	19568,	-26284,	16384,	-28378,	12979,	-30088,	9398,	-31391,
	5690,	-32270,	1905,	-32713,	-1905,	-32713,	-5690,	-32270,
	5690,	-32270,	0,	-32768,	-5690,	-32270,	-11207,	-30792,
	-16384,	-28378,	-21063,	-25102,	-25102,	-21063,	-28378,	-16384,
	19568,	-26284,	18006,	-27377,	16384,	-28378,	14706,	-29283,
	12979,	-30088,	11207,	-30792,	9398,	-31391,	7557,	-31885,
	-9398,	-31391,	-12979,	-30088,	-16384,	-28378,	-19568,	-26284,
	-22487,	-23835,	-25102,	-21063,	-27377,	-18006,	-29283,	-14706,
	-30792,	-11207,	-32270,	-5690,	-32768,	0,	-32270,	5690,
	-30792,	11207,	-28378,	16384,	-25102,	21063,	-21063,	25102,
	5690,	-32270,	3804,	-32546,	1905,	-32713,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	-30792,	-11207,	-31885,	-7557,	-32546,	-3804,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	-16384,	28378,	-11207,	30792,	-5690,	32270,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,

};

// ****************** stage 2 ****************** 

ALIGN(64) const int16_t fft108_tw2[] = 
{
	32767,	0,	31885,	-7557,	32767,	0,	29283,	-14706,
	29283,	-14706,	25102,	-21063,	19568,	-26284,	5690,	-32270,
	19568,	-26284,	12979,	-30088,	-9398,	-31391,	-22487,	-23835,
	5690,	-32270,	-1905,	-32713,	-30792,	-11207,	-32546,	3804,
	-9398,	-31391,	0,	0,	-27377,	18006,	0,	0,
  0,0,0,0,0,0,0,0,

};

// ****************** stage 3 ****************** 

// ****************** stage 4 ****************** 

ALIGN(64) const int16_t fft108_tw4[] = 
{
	32767,	0,	32767,	0,	25102,	-21063,	5690,	-32270,
	5690,	-32270,	-30792,	-11207,	0,	0,	0,	0,

};

// ****************** stage 5 ****************** 

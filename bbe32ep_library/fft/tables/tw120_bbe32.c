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
                  
		       		  
				  
				  
ALIGN(64) const int16_t fft120_tw1[] = 
{
	32767,	0,	32723,	-1715,	32588,	-3425,	32365,	-5126,
	32052,	-6813,	31651,	-8481,	31164,	-10126,	30592,	-11743,
	32767,	0,	32588,	-3425,	32052,	-6813,	31164,	-10126,
	29935,	-13328,	28378,	-16384,	26510,	-19261,	24351,	-21926,
	32767,	0,	32365,	-5126,	31164,	-10126,	29197,	-14876,
	26510,	-19261,	23170,	-23170,	19261,	-26510,	14876,	-29197,
	29935,	-13328,	29197,	-14876,	28378,	-16384,	27482,	-17847,
	26510,	-19261,	25466,	-20622,	24351,	-21926,	23170,	-23170,
	21926,	-24351,	19261,	-26510,	16384,	-28378,	13328,	-29935,
	10126,	-31164,	6813,	-32052,	3425,	-32588,	0,	-32768,
	10126,	-31164,	5126,	-32365,	0,	-32768,	-5126,	-32365,
	-10126,	-31164,	-14876,	-29197,	-19261,	-26510,	-23170,	-23170,
	21926,	-24351,	20622,	-25466,	19261,	-26510,	17847,	-27482,
	16384,	-28378,	14876,	-29197,	13328,	-29935,	11743,	-30592,
	-3425,	-32588,	-6813,	-32052,	-10126,	-31164,	-13328,	-29935,
	-16384,	-28378,	-19261,	-26510,	-21926,	-24351,	-24351,	-21926,
	-26510,	-19261,	-29197,	-14876,	-31164,	-10126,	-32365,	-5126,
	-32768,	0,	-32365,	5126,	-31164,	10126,	-29197,	14876,
	10126,	-31164,	8481,	-31651,	6813,	-32052,	5126,	-32365,
	3425,	-32588,	1715,	-32723,	0,	0,	0,	0,
	-26510,	-19261,	-28378,	-16384,	-29935,	-13328,	-31164,	-10126,
	-32052,	-6813,	-32588,	-3425,	0,	0,	0,	0,
	-26510,	19261,	-23170,	23170,	-19261,	26510,	-14876,	29197,
	-10126,	31164,	-5126,	32365,	0,	0,	0,	0,

};

// ****************** stage 2 ****************** 

ALIGN(64) const int16_t fft120_tw2[] = 
{
	32767,	0,	32052,	-6813,	32767,	0,	29935,	-13328,
	32767,	0,	26510,	-19261,	32767,	0,	21926,	-24351,
	32767,	0,	16384,	-28378,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	29935,	-13328,	26510,	-19261,	21926,	-24351,	10126,	-31164,
	10126,	-31164,	-10126,	-31164,	-3425,	-32588,	-26510,	-19261,
	-16384,	-28378,	-32768,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	21926,	-24351,	0,	0,	-3425,	-32588,	0,	0,
	-26510,	-19261,	0,	0,	-32052,	6813,	0,	0,
	-16384,	28378,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,

};

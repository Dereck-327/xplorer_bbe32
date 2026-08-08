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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
NatureDSP_Baseband library. FFT
Radix-2  FFT on complex_float, float data
IntegrIT, 2006-2017
*/

#ifndef _FFT_FP_TW_H_
#define _FFT_FP_TW_H_

#include "NatureDSP_types.h"

extern const union ufloat32uint32 cfftf16_twd1[24];
extern const union ufloat32uint32 cfftf32_twd1[48];
extern const union ufloat32uint32 cfftf64_twd1[96];
extern const union ufloat32uint32 cfftf64_twd2[24];
extern const union ufloat32uint32 cfftf128_twd1[192];
extern const union ufloat32uint32 cfftf128_twd2[48];
extern const union ufloat32uint32 cfftf256_twd1[384];
extern const union ufloat32uint32 cfftf256_twd2[96];
extern const union ufloat32uint32 cfftf256_twd3[32];
extern const union ufloat32uint32 cfftf512_twd1[768];
extern const union ufloat32uint32 cfftf512_twd2[192];
extern const union ufloat32uint32 cfftf512_twd3[64];
extern const union ufloat32uint32 cfftf1024_twd1[1536];
extern const union ufloat32uint32 cfftf1024_twd2[384];
extern const union ufloat32uint32 cfftf1024_twd3[128];
extern const union ufloat32uint32 cfftf1024_twd4[24];
extern const union ufloat32uint32 cfftf2048_twd1[3072];
extern const union ufloat32uint32 cfftf2048_twd2[768];
extern const union ufloat32uint32 cfftf2048_twd3[256];
extern const union ufloat32uint32 cfftf2048_twd4[48];
extern const union ufloat32uint32 cfftf4096_twd1[6144];
extern const union ufloat32uint32 cfftf4096_twd2[1536];
extern const union ufloat32uint32 cfftf4096_twd3[512];
extern const union ufloat32uint32 cfftf4096_twd4[96];
extern const union ufloat32uint32 cfftf4096_twd5[24];
extern const union ufloat32uint32 cfftf8192_twd1[12288];
extern const union ufloat32uint32 cfftf8192_twd2[3072];
extern const union ufloat32uint32 cfftf8192_twd3[1024];
extern const union ufloat32uint32 cfftf8192_twd4[192];
extern const union ufloat32uint32 cfftf8192_twd5[48];
extern const union ufloat32uint32 cfftf16384_twd1[24576];
extern const union ufloat32uint32 cfftf16384_twd2[6144];
extern const union ufloat32uint32 cfftf16384_twd3[2048];
extern const union ufloat32uint32 cfftf16384_twd4[384];
extern const union ufloat32uint32 cfftf16384_twd5[96];
extern const union ufloat32uint32 cfftf16384_twd6[24];
extern const union ufloat32uint32 cfftf32768_twd1[49152];
extern const union ufloat32uint32 cfftf32768_twd2[12288];
extern const union ufloat32uint32 cfftf32768_twd3[4096];
extern const union ufloat32uint32 cfftf32768_twd4[768];
extern const union ufloat32uint32 cfftf32768_twd5[192];
extern const union ufloat32uint32 cfftf32768_twd6[48];

extern const union ufloat32uint32 rfftf16_twd1[8];
extern const union ufloat32uint32 rfftf32_twd1[16];
extern const union ufloat32uint32 rfftf64_twd1[32];
extern const union ufloat32uint32 rfftf128_twd1[64];
extern const union ufloat32uint32 rfftf256_twd1[128];
extern const union ufloat32uint32 rfftf512_twd1[256];
extern const union ufloat32uint32 rfftf1024_twd1[512];
extern const union ufloat32uint32 rfftf2048_twd1[1024];
extern const union ufloat32uint32 rfftf4096_twd1[2048];
extern const union ufloat32uint32 rfftf8192_twd1[4096];
extern const union ufloat32uint32 rfftf16384_twd1[8192];
extern const union ufloat32uint32 rfftf32768_twd1[16384];

extern const union ufloat32uint32 tfftf8192_twd1[1032];
extern const union ufloat32uint32 tfftf8192_twd2[3072];
extern const union ufloat32uint32 tfftf8192_twd3[1024];
extern const union ufloat32uint32 tfftf8192_twd4[192];
extern const union ufloat32uint32 tfftf8192_twd5[48];
extern const union ufloat32uint32 tfftf16384_twd1[2056];
extern const union ufloat32uint32 tfftf16384_twd2[6144];
extern const union ufloat32uint32 tfftf16384_twd3[2048];
extern const union ufloat32uint32 tfftf16384_twd4[384];
extern const union ufloat32uint32 tfftf16384_twd5[96];
extern const union ufloat32uint32 tfftf16384_twd6[24];
extern const union ufloat32uint32 tfftf32768_twd1[4104];
extern const union ufloat32uint32 tfftf32768_twd2[12288];
extern const union ufloat32uint32 tfftf32768_twd3[4096];
extern const union ufloat32uint32 tfftf32768_twd4[768];
extern const union ufloat32uint32 tfftf32768_twd5[192];
extern const union ufloat32uint32 tfftf32768_twd6[48];


#endif //#ifndef _FFT_FP_TW_H_




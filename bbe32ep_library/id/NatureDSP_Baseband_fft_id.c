/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
 * NatureDSP_Baseband Library API
 * FFT Routines
 * Annotations
 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband Library FFT routines. */
#include "NatureDSP_Baseband_fft.h"
/* Common utility declarations. */
#include "common.h"

ANNOTATE_FUN(cfftas16,    "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas32,    "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas64,    "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas128,   "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas256,   "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas512,   "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas1024,  "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas2048,  "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas4096,  "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas8192,  "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas16384, "Radix-2 forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas32768, "Radix-2 forward FFT (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(cfftas16_norm,    "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas32_norm,    "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas64_norm,    "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas128_norm,   "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas256_norm,   "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas512_norm,   "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas1024_norm,  "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas2048_norm,  "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas4096_norm,  "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas8192_norm,  "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas16384_norm, "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cfftas32768_norm, "Radix-2 forward FFT for normalized data (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(cifftas16,    "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas32,    "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas64,    "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas128,   "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas256,   "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas512,   "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas1024,  "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas2048,  "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas4096,  "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas8192,  "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas16384, "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas32768, "Radix-2 inverse FFT (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(cifftas16_norm,    "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas32_norm,    "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas64_norm,    "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas128_norm,   "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas256_norm,   "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas512_norm,   "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas1024_norm,  "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas2048_norm,  "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas4096_norm,  "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas8192_norm,  "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas16384_norm, "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cifftas32768_norm, "Radix-2 inverse FFT for normalized data (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(tfft8192,   "Radix-2 forward FFT with reduced twiddle table (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(tfft16384,  "Radix-2 forward FFT with reduced twiddle table (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(tfft32768,  "Radix-2 forward FFT with reduced twiddle table (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(tifft8192,  "Radix-2 inverse FFT with reduced twiddle table (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(tifft16384, "Radix-2 inverse FFT with reduced twiddle table (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(tifft32768, "Radix-2 inverse FFT with reduced twiddle table (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(bcfft16,  "Blockwise radix-2 forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcfft32,  "Blockwise radix-2 forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcfft64,  "Blockwise radix-2 forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcfft128, "Blockwise radix-2 forward FFT (complex 16-bit data, no scaling)");

ANNOTATE_FUN(bcifft16,  "Blockwise radix-2 inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcifft32,  "Blockwise radix-2 inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcifft64,  "Blockwise radix-2 inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcifft128, "Blockwise radix-2 inverse FFT (complex 16-bit data, no scaling)");

ANNOTATE_FUN(cnfft12,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft24,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft36,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft48,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft60,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft72,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft96,   "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft108,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft120,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft144,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft180,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft192,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft216,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft240,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft288,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft300,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft324,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft360,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft384,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft432,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft480,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft540,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft576,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft600,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft648,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft720,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft768,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft864,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft900,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft960,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft972,  "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft1080, "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft1152, "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft1200, "Mixed radix forward FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cnfft1536, "Mixed radix forward FFT (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(cinfft12,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft24,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft36,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft48,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft60,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft72,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft96,   "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft108,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft120,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft144,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft180,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft192,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft216,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft240,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft288,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft300,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft324,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft360,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft384,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft432,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft480,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft540,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft576,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft600,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft648,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft720,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft768,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft864,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft900,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft960,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft972,  "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft1080, "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft1152, "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft1200, "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");
ANNOTATE_FUN(cinfft1536, "Mixed radix inverse FFT (complex 16-bit data, auto scaling)");

ANNOTATE_FUN(bcnfft12,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft24,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft36,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft48,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft60,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft72,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft96,   "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft108,  "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcnfft120,  "Blockwise mixed radix forward FFT (complex 16-bit data, no scaling)");

ANNOTATE_FUN(bcinfft12,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft24,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft36,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft48,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft60,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft72,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft96,  "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft108, "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");
ANNOTATE_FUN(bcinfft120, "Blockwise mixed radix inverse FFT (complex 16-bit data, no scaling)");

ANNOTATE_FUN(rfft16,    "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft32,    "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft64,    "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft128,   "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft256,   "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft512,   "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft1024,  "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft2048,  "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft4096,  "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft8192,  "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft16384, "Radix-2 forward FFT (real 16-bit data, auto scaling)");
ANNOTATE_FUN(rfft32768, "Radix-2 forward FFT (real 16-bit data, auto scaling)");

ANNOTATE_FUN(rifft16,    "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft32,    "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft64,    "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft128,   "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft256,   "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft512,   "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft1024,  "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft2048,  "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft4096,  "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft8192,  "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft16384, "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");
ANNOTATE_FUN(rifft32768, "Radix-2 inverse FFT (forms real 16-bit data, auto scaling)");

ANNOTATE_FUN(cfft16f,    "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft32f,    "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft64f,    "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft128f,   "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft256f,   "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft512f,   "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft1024f,  "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft2048f,  "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft4096f,  "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft8192f,  "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft16384f, "Radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(cfft32768f, "Radix-2 forward FFT (complex floating-point data)");

ANNOTATE_FUN(cifft16f,    "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft32f,    "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft64f,    "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft128f,   "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft256f,   "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft512f,   "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft1024f,  "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft2048f,  "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft4096f,  "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft8192f,  "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft16384f, "Radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(cifft32768f, "Radix-2 inverse FFT (complex floating-point data)");

ANNOTATE_FUN(tfft8192f,   "Radix-2 forward FFT with reduced twiddle table (complex floating-point data)");
ANNOTATE_FUN(tfft16384f,  "Radix-2 forward FFT with reduced twiddle table (complex floating-point data)");
ANNOTATE_FUN(tfft32768f,  "Radix-2 forward FFT with reduced twiddle table (complex floating-point data)");

ANNOTATE_FUN(tifft8192f,  "Radix-2 inverse FFT with reduced twiddle table (complex floating-point data)");
ANNOTATE_FUN(tifft16384f, "Radix-2 inverse FFT with reduced twiddle table (complex floating-point data)");
ANNOTATE_FUN(tifft32768f, "Radix-2 inverse FFT with reduced twiddle table (complex floating-point data)");

ANNOTATE_FUN(bcfft16f,  "Blockwise radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(bcfft32f,  "Blockwise radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(bcfft64f,  "Blockwise radix-2 forward FFT (complex floating-point data)");
ANNOTATE_FUN(bcfft128f, "Blockwise radix-2 forward FFT (complex floating-point data)");

ANNOTATE_FUN(bcifft16f,  "Blockwise radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(bcifft32f,  "Blockwise radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(bcifft64f,  "Blockwise radix-2 inverse FFT (complex floating-point data)");
ANNOTATE_FUN(bcifft128f, "Blockwise radix-2 inverse FFT (complex floating-point data)");

ANNOTATE_FUN(brfft16f,  "Blockwise radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(brfft32f,  "Blockwise radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(brfft64f,  "Blockwise radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(brfft128f, "Blockwise radix-2 forward FFT (real floating-point data)");

ANNOTATE_FUN(brifft16f,  "Blockwise radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(brifft32f,  "Blockwise radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(brifft64f,  "Blockwise radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(brifft128f, "Blockwise radix-2 inverse FFT (forms real floating-point data)");

ANNOTATE_FUN(rfft16f,    "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft32f,    "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft64f,    "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft128f,   "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft256f,   "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft512f,   "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft1024f,  "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft2048f,  "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft4096f,  "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft8192f,  "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft16384f, "Radix-2 forward FFT (real floating-point data)");
ANNOTATE_FUN(rfft32768f, "Radix-2 forward FFT (real floating-point data)");

ANNOTATE_FUN(rifft16f,    "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft32f,    "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft64f,    "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft128f,   "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft256f,   "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft512f,   "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft1024f,  "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft2048f,  "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft4096f,  "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft8192f,  "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft16384f, "Radix-2 inverse FFT (forms real floating-point data)");
ANNOTATE_FUN(rifft32768f, "Radix-2 inverse FFT (forms real floating-point data)");

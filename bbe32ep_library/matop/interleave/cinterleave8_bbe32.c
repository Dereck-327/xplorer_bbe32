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
  NatureDSP_Baseband library. Matrix Operations
    M-to-1 complex/real streams interleave
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"
/*-------------------------------------------------------------------------
M-to-1 complex/real streams interleave

Description: interleave element by element M=2,3,4 or 8 streams into a
single stream. Use rinterleave<M>() functions for real data, and
cinterleave<M>() functions - for complex data.

Representation:
<r|c>interleave<M>   16-bit fixed-point data
<r|c>interleave<M>f  IEEE-754 Std single precision floating-point data

Parameters:
Input:
M                    Number of streams
N                    Number of elements per each input stream
x0[N],...,x<M-1>[N]  Input data streams
x[M]                 M pointers to input data streams
Output:
y[M*N]               Interleaved data streams

Restrictions:
x0,...,x<M-1>,
x[0..M-1],y          Must not overlap and must be aligned on 32-byte boundary
N                    Must be a multiple of:
                       16 for real fixed-point data 
                        8 for complex fixed-point data and real floating-point data
                        4 for complex floating-point data
-------------------------------------------------------------------------*/

/* M=8 */
void cinterleave8 ( complex_fract16 * restrict y, const complex_fract16 * restrict x0,
                                          const complex_fract16 * restrict x1,
                                          const complex_fract16 * restrict x2,
                                          const complex_fract16 * restrict x3,
                                          const complex_fract16 * restrict x4,
                                          const complex_fract16 * restrict x5,
                                          const complex_fract16 * restrict x6,
                                          const complex_fract16 * restrict x7, int N )
{
  const xb_vecNx16* restrict w0 = (const xb_vecNx16*)x0;
  const xb_vecNx16* restrict w1 = (const xb_vecNx16*)x1;
  const xb_vecNx16* restrict w2 = (const xb_vecNx16*)x2;
  const xb_vecNx16* restrict w3 = (const xb_vecNx16*)x3;
  const xb_vecNx16* restrict w4 = (const xb_vecNx16*)x4;
  const xb_vecNx16* restrict w5 = (const xb_vecNx16*)x5;
  const xb_vecNx16* restrict w6 = (const xb_vecNx16*)x6;
  const xb_vecNx16* restrict w7 = (const xb_vecNx16*)x7;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  int n;
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x0, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x1, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x2, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x3, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x4, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x5, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x6, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x7, 2 * BBE_SIMD_WIDTH);
  NASSERT((N&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVNX16_IP(X0, w0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, w1, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, w2, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, w3, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X4, w4, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X5, w5, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X6, w6, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X7, w7, 2 * BBE_SIMD_WIDTH);

    DEINTLV1(Y6, Y2, X6, X2);
    DEINTLV1(Y7, Y3, X7, X3);
    DEINTLV1(Y4, Y0, X4, X0);
    DEINTLV1(Y5, Y1, X5, X1);
    DEINTLV1(X6, X4, Y6, Y4);
    DEINTLV1(X7, X5, Y7, Y5);
    DEINTLV1(Y5, Y4, X5, X4);
    DEINTLV1(Y7, Y6, X7, X6);
    DEINTLV1(X2, X0, Y2, Y0);
    DEINTLV1(X3, X1, Y3, Y1);
    DEINTLV1(Y1, Y0, X1, X0);
    DEINTLV1(Y3, Y2, X3, X2);

    BBE_SVNX16_IP(Y0, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y1, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y2, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y3, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y4, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y5, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y6, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y7, z, 2 * BBE_SIMD_WIDTH);
  }
} /* cinterleave8() */

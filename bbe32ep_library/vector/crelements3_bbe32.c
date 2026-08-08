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
  NatureDSP_Baseband library. Vector Mathematics
    Remove the Elements at Given Indices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Remove the Elements at Given Indices 

Description: <r|c>relements<M> copies real or complex elements from the input
vector x to the output vector y, with each M-th element (x[0],x[M],x[2M],...)
being discarded.

Representation: 16-bit fixed-point data

Parameters:
Input:
x[N*M]      Input data vector
N           Size of input data vector divided by M
Output:
y[N*(M-1)]  Output data vector

Restrictions:
x,y         Aligned on 32-byte boundary
x,y         Must not overlap
N           Multiple of 16 for real data, or a multiple of 8 for complex data
-------------------------------------------------------------------------*/

void crelements3 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N )
{
    int n;
    static const int16_t ALIGN(32) SEL[][BBE_SIMD_WIDTH] =
    {
        { 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 16, 17, 20, 21, 22, 23 },
        { 10, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 28, 29, 30, 31 }
    };
    const xb_vecNx16 * restrict pX0 = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pX1 = pX0 + 3;
    const xb_vecNx16 * restrict pS = (const xb_vecNx16 *)SEL;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    xb_vecNx16 x0, x1, x2, t0, t1;
    vselN      sel0, sel1;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

    x0 = BBE_LVNX16_I(pS, 0 * 2 * BBE_SIMD_WIDTH);
    x1 = BBE_LVNX16_I(pS, 1 * 2 * BBE_SIMD_WIDTH);
    sel0 = BBE_MOVVSELNX16(x0, 0);
    sel1 = BBE_MOVVSELNX16(x1, 0);

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX0, 4 * 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SELNX16(x1, x0, sel0);
        BBE_SVNX16_IP(t0, pY, 2 * BBE_SIMD_WIDTH);
        t1 = BBE_SELNX16(x2, x1, sel1);
        BBE_SVNX16_IP(t1, pY, 2 * BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(x0, pX1, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX1, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX1, 4 * 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SELNX16(x1, x0, sel0);
        BBE_SVNX16_IP(t0, pY, 2 * BBE_SIMD_WIDTH);
        t1 = BBE_SELNX16(x2, x1, sel1);
        BBE_SVNX16_IP(t1, pY, 2 * BBE_SIMD_WIDTH);
    }
    if (N & (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVNX16_IP(x0, pX0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX0, 4 * 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SELNX16(x1, x0, sel0);
        BBE_SVNX16_IP(t0, pY, 2 * BBE_SIMD_WIDTH);
        t1 = BBE_SELNX16(x2, x1, sel1);
        BBE_SVNX16_IP(t1, pY, 2 * BBE_SIMD_WIDTH);
    }
} /* crelements3() */

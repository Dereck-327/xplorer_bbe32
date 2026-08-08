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
    Move the Elements at Given Indices
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
Move the Elements at Given Indices 

Description: <r|c>melements<M> Functions retrieve each M-th real or complex
element of the input vector x (x[0],x[M],x[2*M],...) and stores it to the
output vector y.

Representation: 16-bit fixed-point data

Parameters:
Input:
x[N*M]  Input data vector
N       Size of output data vector
Output:
y[N]    Output data vector

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap
N       Multiple of 16 for real data, or a multiple of 8 for complex data
-------------------------------------------------------------------------*/

void cmelements3(complex_fract16 * restrict y, const complex_fract16 * restrict x, int N)
{
    int n;
    static const int16_t ALIGN(32) SEL[][BBE_SIMD_WIDTH] =
    {
        { 0, 1, 6, 7, 12, 13, 18, 19, 24, 25, 30, 31, 4, 5, 10, 11 },
        { 2, 3, 8, 9, 14, 15, 20, 21, 26, 27, 0, 1, 6, 7, 12, 13 },
        { 4, 5, 10, 11, 16, 17, 22, 23, 28, 29, 2, 3, 8, 9, 14, 15 }
    };
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    xb_vecNx16 x0, x1, x2, x3, x4, x5, t0, t1, t2;
    valign y_zalign;
    vselN sel0, sel1, sel2;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

    x0 = BBE_LVNX16_I((const xb_vecNx16*)SEL, 0 * 2 * BBE_SIMD_WIDTH);
    x1 = BBE_LVNX16_I((const xb_vecNx16*)SEL, 1 * 2 * BBE_SIMD_WIDTH);
    x2 = BBE_LVNX16_I((const xb_vecNx16*)SEL, 2 * 2 * BBE_SIMD_WIDTH);
    sel0 = BBE_MOVVSELNX16(x0, 0);
    sel1 = BBE_MOVVSELNX16(x1, 0);
    sel2 = BBE_MOVVSELNX16(x2, 0);

    y_zalign = BBE_ZALIGN();

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x3, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x4, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x5, pX, 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SELNX16(x1, x0, sel0);
        BBE_SAVNX16_XP(t0, y_zalign, pY, 2 * 12);
        t1 = BBE_SELNX16(x3, x2, sel2);
        BBE_SAVNX16_XP(t1, y_zalign, pY, 2 * 10);
        t2 = BBE_SELNX16(x5, x4, sel1);
        BBE_SAVNX16_XP(t2, y_zalign, pY, 2 * 10);
    }
    if (N & (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX, 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SELNX16(x1, x0, sel0);
        BBE_SAVNX16_XP(t0, y_zalign, pY, 2 * 12);
        t1 = BBE_SELNX16(x2, x2, sel2);
        BBE_SAVNX16_XP(t1, y_zalign, pY, 2 * 4);
    }
    BBE_SAVNX16POS_FP(y_zalign, pY);
} /* cmelements3() */

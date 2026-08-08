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

void rrelements24 ( int16_t* restrict y, const int16_t* restrict x, int N )
{
    int n;
    static const int16_t ALIGN(32) SEL[BBE_SIMD_WIDTH] = { 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17 };
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pS = (const xb_vecNx16 *)SEL;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    xb_vecNx16 x0, x1, x2, y0, y1, y2;
    valign y_align;
    vselN sel0;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N%BBE_SIMD_WIDTH == 0);

    y_align = BBE_ZALIGN();
    x0 = BBE_LVNX16_I(pS, 0 * 2 * BBE_SIMD_WIDTH);
    sel0 = BBE_MOVVSELNX16(x0, 0);

    for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 3)); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, pX, 2 * BBE_SIMD_WIDTH);
        y0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
        BBE_SAVNX16_XP(y0, y_align, pY, 2 * BBE_SIMD_WIDTH);
        y1 = BBE_SELNX16(x2, x1, sel0);
        BBE_SAVNX16_XP(y1, y_align, pY, 2 * BBE_SIMD_WIDTH);
        y2 = BBE_SELNX16I(x2, x2, BBE_SELI_ROTATE_RIGHT_2);
        BBE_SAVNX16_XP(y2, y_align, pY, 2 * 14);
    }
    BBE_SAVNX16POS_FP(y_align, pY);
} /* rrelements24() */

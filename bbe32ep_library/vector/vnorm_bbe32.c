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
    Vector Normalization
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
Vector Normalization

Description: This function scale all data in the input vector x by a power
of 2 factor and returns results in the output vector y. Most often it is used
in conjunction with common block exponent function vbexp() to normalise vector
data, i.e. to shift all elements of a vector to the left such that the minimum
number of redundant sign bits over the output vector is zero. It can, however,
be utilised for shifting data by arbitrary number of bit positions: for positive
shift amount it is a saturating left shift, for negative shift amount it is a
sign-extending right shift.

Representation:
vnorm  16-bit signed fixed-point data

Parameters:
x[N]   Input data
t      Shift amount, [-16..16]
N      Length of input/output vectors
Output:  
y[N]   Output data

Restrictions:
x,y    Aligned on 32-byte boundary
x,y    Must not overlap
N      Multiple of 16
-------------------------------------------------------------------------*/

void vnorm ( int16_t * restrict y,
       const int16_t * restrict x,
       int t,
       int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    int k;
    xb_vecNx16 x0, y0;
    vsaN       shft;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(t >= -16 && t <= 16);
    NASSERT(N>0 && N%BBE_SIMD_WIDTH == 0);

    shft = BBE_MOVVSA32(t);
    for (k = 0; k<(N >> LOG2_BBE_SIMD_WIDTH); k++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        y0 = BBE_SLSNX16(x0, shft);
        BBE_SVNX16_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    }
} /* vnorm() */

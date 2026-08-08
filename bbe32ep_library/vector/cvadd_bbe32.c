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
    Complex Vector Sum
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
Complex Vectors Sum

Description: These routines perform pairwise summation of complex vectors.

Representation:
cvadd   Signed fixed-point format. Input vectors are comprised of 32-bit
        complex elements with 16-bit real and imaginary components. Output
        vector elements are 32-bit complex values with 16-bit saturated
        real and imaginary parts.
cvaddf  IEEE-754 Std. single precision floating-point format for real/imaginary
        components of 64-bit input/output data

Parameters:
Input:
x[2*N]  Input complex data,Q15
y[2*N]  Input complex data,Q15
N       Length of vectors
Output:
z[2*N]  Output complex data,Q15

Restrictions:
z,x,y   Must not overlap
x,y,z   Aligned on 32-byte boundary
N       Multiple of 8 (cvadd) or 4 (cvaddf) 
-------------------------------------------------------------------------*/

void cvadd ( complex_fract16 * restrict z,
        const complex_fract16 * restrict x,
        const complex_fract16 * restrict y,
        int N)
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pY = (const xb_vecNx16 *)y;
    xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
    int n;
    xb_vecNx16 x0, y0, z0;
    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % (BBE_SIMD_WIDTH / 2) == 0);
    for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

        z0 = BBE_ADDSNX16(x0, y0);
        BBE_SVNX16_IP(z0, pZ, 2 * BBE_SIMD_WIDTH);
    }
} /* cvadd() */

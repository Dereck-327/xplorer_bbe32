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
    Fixed to Floating-Point Conversion
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Fixed to Floating-Point Conversion

Description: Function converts 16-bit signed fixed-point data with 
user-defined point position to 16-bit half-precision floating-point
format (IEEE 754-2008 binary16). 

Notes:
1. IEEE 754-2008 binary16 encoding format provides 1 sign bit, 10 trailing
   significand bits and 5 exponent bits.
2. If absolute input value |x|*2^-q exceeds 65504, the result is saturated to
   +/-Inf (0x7C00 or 0xFC00). In particular, if point position is less than
   or equal to -16 (q<=-16), any nonzero fixed-point value is converted to
   +/-Inf.
3. If absolute input value belongs to the closed range of [2^-24,2^-15], then
   the conversion results in a subnormal half-precision floating-point number
   as defined by IEEE 754-2008.
4. If input value x*2^-q belongs to the open range of (-2^-24,2^-24), then it
   is converted to zero of original sign (0x0000 or 0x8000). Specifically, any
   fixed-point value is converted to +/-0 when the point position is greater
   than or equal to 40 (q>=40).
5. If input value cannot be exactly represented in half-precision floating-
   point format, then it is rounded to the next representable value toward zero.

Parameters:
Input:
x[N]    Input fixed-point data, Q(q)
q       Point position for input data
N       Size of input/output arrays
Output:
y[N]    Output half-precision floating-point values

Restrictions:
y,x     Aligned on 32-byte boundary
y,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/

void fix2hf ( int16_t * restrict y, 
        const int16_t * restrict x, 
        int q,
        int N )
{
    const xb_vecNx16 *restrict X = (const xb_vecNx16 *)x;
    xb_vecNx16 *restrict Y = (xb_vecNx16 *)y;
    int i;
    xb_vecNx16 _3FF = BBE_MOVVA16(0x3FF);
    xb_vecNx16 _8000 = BBE_MOVVINX16(BBE_MOVVI_INT16_MININT);
    xb_vecNx16 _30q = BBE_MOVVA16(30 - q);
    xb_vecNx16 _31 = BBE_MOVVINT16(31);
    xb_vecNx16 zero = 0;
    xb_vecNx16 _1 = BBE_MOVVINT16(1);
    xb_vecNx16 t, vx_sgn, vx_mnt, vx_exp;
    vboolN b;
    vsaN nsa;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % BBE_SIMD_WIDTH == 0);
    //NASSERT(q >= -16 && q <= 42);

    __Pragma("loop_count min=1");
    for (i = 0; i<(N >> LOG2_BBE_SIMD_WIDTH); ++i)
    {
        BBE_LVNX16_XP(vx_mnt, X, 2 * BBE_SIMD_WIDTH);
        // take sign, absolute value and normalize
        vx_sgn = BBE_ANDNX16(vx_mnt, _8000);
        vx_mnt = BBE_ABSNX16(vx_mnt);
        nsa = BBE_NSAUNX16(vx_mnt);
        vx_exp = BBE_MOVVVS(nsa);
        vx_mnt = BBE_SLLNX16(vx_mnt, nsa);
        vx_mnt = BBE_SRLINX16(vx_mnt, 5);
        // renormalize using q
        vx_exp = BBE_SUBNX16(_30q, vx_exp);
        b = BBE_EQNX16(vx_mnt, zero);           // if (x_mnt==0) x_exp = 0;
        vx_exp = BBE_MOVNX16T(zero, vx_exp, b); // if (x_mnt==0) x_exp = 0;
        t = BBE_SUBNX16(_1, vx_exp);
        nsa = BBE_MOVVSV(t, 0);
        t = BBE_SRLNX16(vx_mnt, nsa);
        vx_exp = BBE_MINNX16(vx_exp, _31);
        vx_exp = BBE_MAXNX16(vx_exp, zero);
        // overflow 
        b = BBE_EQNX16(vx_exp, _31);
        vx_mnt = BBE_MOVNX16T(zero, vx_mnt, b);
        // subnormal number 
        b = BBE_EQNX16(vx_exp, zero);
        vx_mnt = BBE_MOVNX16T(t, vx_mnt, b);
        // combine matissa, exponenta and the sign 
        vx_mnt = BBE_ANDNX16(vx_mnt, _3FF);
        vx_exp = BBE_SLLINX16(vx_exp, 10);
        vx_mnt = BBE_ADDNX16(vx_mnt, vx_exp);
        vx_mnt = BBE_ADDNX16(vx_mnt, vx_sgn);
        BBE_SVNX16_IP(vx_mnt, Y, 2 * BBE_SIMD_WIDTH);
    }
} /* fix2hf() */

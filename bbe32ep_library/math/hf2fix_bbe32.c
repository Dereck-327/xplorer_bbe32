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
    Floating to Fixed-Point Conversion
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"

#ifdef BBE_MAXVSN
#define MAXVSN(shft,sft,ShiftMinVsa)     \
{                                        \
    shft = BBE_MOVVSV(sft, 0);           \
    shft = BBE_MAXVSN(shft, ShiftMinVsa); \
}
#else
#define MAXVSN(shft,sft,ShiftMinVsa)                  \
{                                                     \
    xb_vecNx16 maxShft;                               \
    maxShft = BBE_MAXNX16(sft, BBE_MOVVVS(ShiftMinVsa)); \
    shft = BBE_MOVVSV(maxShft, 0);                    \
}
#endif

/*-------------------------------------------------------------------------
Floating to Fixed-Point Conversion

Description: Function converts 16-bit half-precision floating-point data
(IEEE 754-2008 binary16 format) to 16-bit signed fixed-point format with
user-defined point position.

Notes:
1. +/-Inf is converted to 32767 or -32768, respectively.
2. NaN values (zero exponent field, non-zero significand field) are treated as
   +/-Inf depending on the sign bit.
3. If a finite input value is too large in its magnitude to be represented
   in the fixed-point format with user-specified point position, then
   it is saturated to -32768 or 32767 depending on the input sign.
4. If the input value cannot be exactly represented in fixed-point format
   with user-specified point position, then it is rounded to the next
   representable value toward -Inf.
5. The sensible range for the point position q is [-15,39]. That is, for q<=-16
   the result of conversion is always zero or -1 (depending on the input sign),
   unless the input value is +/-Inf or NaN. For q>=39, any non-zero floating-
   point value is converted to either -32768 or 32767, depending on the input sign.
6. Given an arbitrary floating-point value, the optimum point position for the
   floating-to-fixed conversion is q = 14-e, where e is the unbiased exponent
   of the floating-point value. For a normal positive floating-point value, the
   optimum point position provides a normalized fixed-point result. For a
   normal negative value on input, the number of redundant sign bits for the
   conversion result is at most 1.

Parameters:
Input:
x[N]    Input half-precision floating-point values 
q       Point position for output data
N       Size of input/output arrays
Output:
y[N]    Output fixed-point data, Q(q)

Restrictions:
y,x     Aligned on 32-byte boundary
y,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/

void hf2fix ( int16_t * restrict y, 
       const int16_t * restrict x, 
       int q,
       int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;

    xb_vecNx16 x0, y0, exp;
    xb_vecNx16 m0;
    xb_vecNx16 sft, ShiftMax;
    vsaN       shft, ShiftMinVsa;
    xb_vecNx16 t0, t1;
    xb_vecNx16 z0, _1f, _0, z1;
    vboolN     bZeroExp, bMaxExp;
    int n;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % BBE_SIMD_WIDTH == 0);
    q = (q < -16) ? (-16) : ((q > 42) ? (42) : (q));

    z0 = BBE_MOVVA16(1023);
    _1f = BBE_MOVVA16(0x001f); // 31
    _0 = BBE_MOVVA16(0);
    z1 = BBE_MOVVA16(1024);

    // Exponent encoded in a floating point number and requested fixed-point
    // position will be applied in a single bit shift. Shift amount at the
    // most significant element of tbl1 (for a special exponent value of 31)
    // guarantees a saturated output value regardless of actual significand
    // value and fixed-point position, thus -/+inf and qNaN, sNaN at input
    // result in -32768/32767 at output.
    t1 = BBE_MOVVA16(q - 10 - 15);
    ShiftMax = BBE_MOVVA16(q - 10 + 31);
    ShiftMinVsa = BBE_MOVVSA32(q - 10 - 14);

    __Pragma("ymemory(pY)");
    __Pragma("loop_count min=1");
    for (n = 0; n<N / BBE_SIMD_WIDTH; n++)
    {
        // Load SIMD_WIDTH floating point numbers.
        BBE_LVNX16_XP(x0, pX, +2 * BBE_SIMD_WIDTH);

        exp = BBE_SRAINX16(x0, 10);
        exp = BBE_ANDNX16(exp, _1f);
        bZeroExp = BBE_EQNX16(exp, _0);
        bMaxExp = BBE_EQNX16(exp, _1f);

        // Select the bit shift amount and leading bit value. The leading bit
        // is always 1 unless a floating point number is subnormal (zero exponent).
        sft = BBE_ADDNX16(exp, t1); // convert exp range "0, ... 31" to "-14, -14, -13, -12, ... 13, 14, 15, 31" and add to result value "q - 10"
        sft = BBE_MOVNX16T(ShiftMax, sft, bMaxExp);
        MAXVSN(shft, sft, ShiftMinVsa);

        // Extract the 10-bit significand field value.
        m0 = BBE_ANDNX16(x0, z0);
        // Add the leading bit value (at position 10). (add lag: 0 or 1024)
        BBE_ADDNX16F(m0, m0, z1, bZeroExp);
        // Account for the sign.
        t0 = BBE_MULSGNNX16(x0, m0);

        // Apply the bit shift. For +/-inf and qNaN, sNaN this results in a value
        // greater than 32767 or less than or equal to -32768. The shift amount
        // varies between -40..63!
        y0 = BBE_SLSNX16(t0, shft);

        // Save SIMD_WIDTH 16-bit fixed point values. Saturation is applied
        // to convert from 20-bit to 16-bit, so that +/-inf eventually turn
        // into 32767/-32768.
        BBE_SVNX16_IP(y0, pY, +2 * BBE_SIMD_WIDTH);
    }
} /* hf2fix() */

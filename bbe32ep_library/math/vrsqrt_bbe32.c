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
    Inverse Square Root
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
Reciprocal Square Root

Description: Evaluate the reciprocal square root of input value x and store
result to y: y = 1/x^0.5.

Representation:
vrsqrt,srsqrt  32-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-
               defined, provided that it is even.
               Fixed-point format Qy for output data is Qy = 31-Qx/2.
               For example, if Qx == 30 then Qy = 31-30/2 = 16
vrsrtf         IEEE-754 Std. single precision floating-point format

Accuracy:
vrsqrt,srsqrt  1.1e-4 (worst case relative error)
vrsrtf,srsqrtf 2 ULP
vfastrsrtf     3 ULP

Notes (for non-fast versions)::
1. Fixed-point functions return MIN_INT32 (0x80000000) for a negative or
   zero input value.
2. Floating-point reciprocal square root conforms to IEEE-754 Std rSqrt 
   operation in respect to signaling error conditions by means of floating-
   point exceptions.
3. Floating-point reciprocal square root limits the range of allowable
   input values, as follows:
   A) If x<0, functions raise the "invalid" floating-point exception,
      assign EDOM to errno and set output value y to NaN.
   B) If x==+/-0, functions set output value y to +/-HUGE_VALF, raise the
      "divide by zero" floating-point exception, and assign ERANGE to errno.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vrsqrt) or 8 (vrsqrtf,vfastrsqrtf)
-------------------------------------------------------------------------*/

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)

DISCARD_FUN(void, vrsqrt, (int32_t * restrict y,
    const int32_t * restrict x,
    int N))

#else

void vrsqrt(int32_t * restrict y,
              const int32_t * restrict x,
              int N )
{
    int n;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    xb_vecNx40   a_vec, x_vec;
    xb_vecNx16   y0, y1;
    xb_vecNx16   _0, _8000, x0, x1, b_vec, cn_vec;
    vsaN nsa;
    vboolN    b0, b1;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N%BBE_SIMD_WIDTH == 0);

    _0 = 0;
    _8000 = BBE_MOVVINX16(BBE_MOVVI_CQ15_MI);
    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(x1, pX, 2 * BBE_SIMD_WIDTH);
        x_vec = BBE_MOVSWV(x1, x0);
        x0 = BBE_SELNX16I(x1, x0, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
        x1 = BBE_SHFLNX16I(x0, BBE_SHFLI_DOUBLE_1_HI);
        x0 = BBE_SHFLNX16I(x0, BBE_SHFLI_DOUBLE_1_LO);
        b0 = BBE_LTNX16(x0, _0);
        b1 = BBE_LTNX16(x1, _0);
        x_vec = BBE_ADDNX40(x_vec, x_vec);
        nsa = BBE_NSAENX40(x_vec);
        x_vec = BBE_SLLNX40(x_vec, nsa);
        BBE_RSQRTLUNX40_0(a_vec, b_vec, cn_vec, x_vec);
        BBE_RSQRTLUNX40_1(a_vec, b_vec, cn_vec, x_vec);
        BBE_MULUUSNX16(a_vec, cn_vec, b_vec);
        nsa = BBE_SUBSR1SAVSN(26, nsa);
        a_vec = BBE_SRANX40(a_vec, nsa);
        y1 = BBE_MOVSVWH(a_vec);
        y0 = BBE_MOVSVWL(a_vec);
        y0 = BBE_MOVNX16T(_8000, y0, b0);
        y1 = BBE_MOVNX16T(_8000, y1, b1);
        // Here BBE_SVNX16_XP provide better schedule than BBE_SVNX16_IP!
        BBE_SVNX16_XP(y0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_XP(y1, pY, 2 * BBE_SIMD_WIDTH);
    }
} /* vrsqrt() */

#endif

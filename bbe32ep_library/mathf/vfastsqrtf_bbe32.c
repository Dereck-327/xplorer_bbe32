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
  NatureDSP_Baseband library. Math functions
    Square Root
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Square Root

Description: These functions compute the positive square root of input
vector elements.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
vsqrtf, ssqrtf  2 ULP 
vfastsqrtf      3 ULP

Notes:
1. Square root functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. For a negative input value, functions raise the "invalid" floating-point
   exception, assign EDOM to errno and set the respective output value z to NaN.
3. Negative zero (-0) is considered as a valid input, the result is also -0.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vfastsqrtf,( float32_t * restrict y,
        const float32_t * restrict x,
        int N ))
#else
void vfastsqrtf ( float32_t * restrict y,
        const float32_t * restrict x,
        int N )
{
  int n;
  xb_vecN_2xf32 x0, t0, t1, one, half;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  half = BBE_CONSTN_2XF32(3);
  one = BBE_CONSTN_2XF32(1);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 r0, r0_err,z0;
    xb_vecN_2xf32 x0_adj, x0_red;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    x0_red = BBE_NEXP01N_2XF32(x0);/* negated x with reduced exponent range */

    /* Initial rsqrt approximation with exponent range reduction */
    r0 = BBE_SQRT0N_2XF32(x0);
    r0 = BBE_NEGN_2XF32(r0);
    /* compute approximation error */
    t0 = BBE_MULN_2XF32(r0, r0);
    t1 = x0_red; BBE_ADDEXPN_2XF32(t1, half);/* -0.5*x */
    r0_err = half; BBE_MULANN_2XF32(r0_err, t0, t1);/* approximation error is (0.5-0.5*x*r*r) */
    BBE_MULANN_2XF32(r0, r0, r0_err);/* Second recip sqrt approximation */

    /* Compute reduced range sqrt approximation */
    z0 = BBE_MULN_2XF32(x0_red, r0);/* z = x*rsqrt(x) */
    /* Make final adjustment and restore range */
    x0_adj = BBE_MKSADJN_2XF32(x0);
    t0 = x0_red;
    BBE_MULANN_2XF32(t0, z0, z0);
    BBE_ADDEXPMN_2XF32(z0, x0_adj);
    t1 = BBE_MULN_2XF32(half, r0);
    BBE_ADDEXPN_2XF32(t1, x0_adj);
    BBE_DIVNN_2XF32(z0, t0, t1);
    BBE_SVN_2XF32_IP(z0, pY, 2 * BBE_SIMD_WIDTH);
  }
} /* vfastsqrtf() */
#endif

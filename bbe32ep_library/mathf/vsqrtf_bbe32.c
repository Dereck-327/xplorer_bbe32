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
#include <errno.h>
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h" 
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
DISCARD_FUN(void,vsqrtf,( float32_t * restrict y,
        const float32_t * restrict x,
        int N ))
#else
void vsqrtf ( float32_t * restrict y,
        const float32_t * restrict x,
        int N )
{
  int n;
  xb_vecN_2xf32 x0, y0, t0, t1, one, half;
  xb_vecN_2x32Uv flags, flags_n;
  __fenv_t fenv;
  vboolN_2 b_edom, b_fe_inv;
  vboolN_2 b_nan;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  half = BBE_CONSTN_2XF32(3);
  one = BBE_CONSTN_2XF32(1);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    {
      xb_vecN_2xf32 r0, r0_err,z0;
      xb_vecN_2xf32 x0_adj, x0_red;
     
      x0_red = BBE_NEXP01N_2XF32(x0);/* negated x with reduced exponent range */

      /* Initial rsqrt approximation with exponent range reduction */
      r0 = BBE_SQRT0N_2XF32(x0);
      r0 = BBE_NEGN_2XF32(r0);
      /* compute approximation error */
      t0 = BBE_MULN_2XF32(r0, r0);
      t1 = x0_red;
      BBE_ADDEXPN_2XF32(t1, half);/* -0.5*x */
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

      y0=z0;
    }
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    /*
    * Perform additional analysis of input data for Error Handling.
    */
    b_nan = BBE_UNN_2XF32(x0, x0);
    BBE_ABSN_2XF32T(x0, x0, b_nan);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);
  }
  {
    xb_vecN_2xf32 v_edom, v_fe_inv;
    int fe_inv, en_edom;
    xb_vecN_2x32Uv maskx;
    fe_inv = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan, v_ltz;
    vboolN_2 b_snan, b_nan, b_ltz;


    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x1));
    v_ltz = BBE_ANDN_2X32U(flags, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_ltz = BBE_NEQN_2X32U(v_ltz, maskx);

    b_edom = BBE_ORBN_2(b_ltz, b_snan);
    b_edom = BBE_ORBN_2(b_edom, b_nan);

    b_fe_inv = BBE_ORBN_2(b_ltz, b_snan);

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_fe_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_fe_inv);

    /* Retrieve EDOM state: some x<0 OR x==sNaN OR x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    /* Merge FE_INVALID state: some x<0 OR x==sNaN (the latter is detected by hardware) */
    fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);

    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };

    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise FE_INVALID (x<0 or x==sNaN) and/or FE_DIVBYZERO (x==0). */
    __feraiseexcept(fe_inv );
  }
} /* vsqrtf() */
#endif

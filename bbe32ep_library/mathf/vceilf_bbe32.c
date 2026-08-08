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
    Floating-Point Ceil
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Floating-Point Ceil

Description: These functions return the smallest integral value that is not
less than an input value.

Data format: IEEE-754 Std. single precision floating-point.

Note:
Ceil functions conform to ANSI C requirements on standard math library
functions in respect to treatment of errno and floating-point exceptions.

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
DISCARD_FUN(void,vceilf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vceilf ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0, x1, y0, y1, v_edom;
  vboolN_2 b_nan;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  xb_vecN_2xf32  * restrict pY = (xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  v_edom = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH )); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x1, pX, 2 * BBE_SIMD_WIDTH);

    y0 = BBE_FICEILN_2XF32(x0);
    y1 = BBE_FICEILN_2XF32(x1);
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(y1, pY, 2 * BBE_SIMD_WIDTH);

    b_nan = BBE_UNN_2XF32(x1, x0);
    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), v_edom, b_nan);

  }
  if (N&(BBE_SIMD_WIDTH/2))
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), v_edom, b_nan);
    y0 = BBE_FICEILN_2XF32(x0);
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
  }
  {
    int en_edom;
    en_edom = 0;

    /* Retrieve EDOM state:  x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);

    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };
  }
} /* vceilf() */
#endif

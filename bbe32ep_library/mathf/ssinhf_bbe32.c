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
    Hyperbolic Sine
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
/* Table for polynomial approximation to sinh(x). */
#include "sinhf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* Constants and polynomial coeffs for exp(x) approximation. */
#include "expf_tbl.h"
#include "pow2f_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Hyperbolic Sine

Description: These functions compute hyperbolic sine of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vsinhf(), ssinhf()
3 ULP for vfastsinhf()

Notes for non-fast versions:
1. Hyperbolic sine functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Due to limited dynamic range of single precision floating-point format,
   hyperbolic sine result for an input value x such that |x|>89.41599 is
   sign(x)*HUGE_VALF.

Input domain for 'fast' version vfastsinhf():
|x|<89.41599
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
#if HAVE_VFPU
float32_t ssinhf ( float32_t x )
{
  float32_t one, quarter;

  float32_t r, z, p0, dy, y1;
  int sx;
  xb_int32v e1, e2;
  /* Powers of input value; floating-point polynomial value */
  float32_t y;
  vbool1 b_le1, b_inf;
  xtfloat big = BBE_MOVF32_FROM32(0x7f7f7f7f); /* 0x7f7f7f7f, 3.4e38f */
  xtfloat inf = BBE_MOVF32_FROM32(0x7F800000);

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  sx = XT_RFR(x);
  x = XT_ABS_S(x);
  if (BBE_MOVAB1(XT_UN_S(x, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, x);
    return x;
  }
  b_le1 = XT_OLT_S(sinhf_maxarg.f, x);
  if (BBE_MOVAB1(b_le1))
  {
    __Pragma("frequency_hint never");
    b_inf = XT_OEQ_S(x, inf);
    if (BBE_MOVAB1(~b_inf))
    {
      errno = ERANGE;
    }
    x = big;
    XT_MOVT_S(x, inf, b_inf);
    XT_MOVLTZ_S(x, -x, sx);
    x = XT_ADD_S(x, x);    /* to generate FE_OVERFLOW when needed */
    return x;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  one = (float32_t)XT_CONST_S(1);
  quarter = (float32_t)XT_CONST_S(3);
  quarter = XT_MUL_S(quarter, quarter);
  if (x>one)
  {
    /* For a large input value, sinh(x) is computed from exp(|x|)/2. */

    /* z <- exp(x)/2 */
    /* log2(n) */
    /* scale input to 1/ln(2) */
    p0 = XT_MUL_S(x, log2_e[0].f);
    p0 = XT_FIROUND_S(p0);
    dy = XT_NEG_S(p0);
    XT_MADD_S(dy, x, log2_e[0].f);
    XT_MADD_S(dy, x, log2_e[1].f);
    /* compute 2^x */
    /* Load and replicate polynomial coefficients */

    y = pow2f_coef[0].f; y1 = pow2f_coef[1].f;
    XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[2].f;
    XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[3].f;
    XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[4].f;
    XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[5].f;
    XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[6].f;
    XT_MADD_S(y1, dy, y); y = y1;
    /* resulted scaling */
    p0 = XT_MAX_S(XT_MIN_S(p0, 129.f), -151.f);

    /* Apply exponential part to the result */
    {
      xb_vecNx16 tmp, v1, v2;
      tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(BBE_MOVN_2XF32_FROMF32(p0), 0));
      tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254 - 30 - 1));
      v1 = BBE_SRLINX16(tmp, 1);
      v2 = BBE_SUBNX16(tmp, v1);
      v1 = BBE_SLLINX16(v1, 7);
      v2 = BBE_SLLINX16(v2, 7);
      v1 = BBE_SELNX16I(v1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      v2 = BBE_SELNX16I(v2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      e1 = BBE_MOV32_FROMN_2X32(BBE_MOVN_2X32_FROMNX16(v2));
      e2 = BBE_MOV32_FROMN_2X32(BBE_MOVN_2X32_FROMNX16(v1));
    }
    /*
    * Convert (y*2^(ex-30))/2 to floating-point p == exp(x)/2
    */
    r = XT_MUL_S(y, 1073741824.f);
    y = XT_MUL_S(r, BBE_MOVF32_FROM32(e2));
    z = XT_MUL_S(y, BBE_MOVF32_FROM32(e1));

    r = XT_RECIP_S(z);   /* Initial approximation for 1/z */

    /* z <- z + 0.25/z == exp(x)/2 - exp(-x)/2 */
    XT_MSUB_S(z, quarter, r);

  }
  else
  {
    float32_t x2, t;
    /*
    * Use polynomial approximation for small input values. This branch is
    * also used for a NaN on input.
    */

    x2 = XT_MUL_S(x, x);
    z = polysinhf_tbl[0].f;
    t = polysinhf_tbl[1].f; XT_MADD_S(t, x2, z); z = t;
    t = polysinhf_tbl[2].f; XT_MADD_S(t, x2, z); z = t;
    z = XT_MUL_S(x2, z);
    XT_MADD_S(x, x, z); z = x;

  }
  /* Restore the input sign. */
  XT_MOVLTZ_S(z, -z, sx);
  BBE_MOVSCFV(SCF);
  return z;
} /* ssinhf() */
#else
DISCARD_FUN(float32_t,ssinhf,( float32_t x ))
#endif

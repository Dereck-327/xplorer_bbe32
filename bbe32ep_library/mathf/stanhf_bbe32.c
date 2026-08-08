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
    Hyperbolic Tangent
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
/* Polynomial coeffs for tanh(x) approximation */
#include "tanhf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* Constants and polynomial coeffs for exp(x) approximation. */
#include "expf_tbl.h"
#include "pow2f_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Hyperbolic Tangent

Description: These functions compute hyperbolic tangent of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vtanhf(), stanhf()
3 ULP for vfasttanhf()

Note for non-fast version
Hyperbolic tangent functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-
point exceptions.

Input domain for 'fast' version vfasttanhf():
|x|<Inf
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
float32_t stanhf ( float32_t x )
{
  float32_t zero, one, two, half, z, r, eps;
  float32_t y;
  float32_t p0, dy, y1;
  xb_int32v ux;
  xb_int32v e1, e2;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  if (BBE_MOVAB1(XT_UN_S(x, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, x);
    return x;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  zero = (float32_t)XT_CONST_S(0);
  one = (float32_t)XT_CONST_S(1);
  two = (float32_t)XT_CONST_S(2);
  half = (float32_t)XT_CONST_S(3);
  ux = BBE_MOV32_FROMF32(x); 
  ux = BBE_OPERATOR_AND32(ux, 0x80000000);
  x = XT_ABS_S(x);
  if (x > halfln3.f)
  {
    /*
    * For a large input value tanh(x) is computed from exp(2*x)/2, using
    * the following identity: tanh(x) == 1 - 2/(exp(2*x)+1)
    */
    r = zero; XT_MADDN_S(r, two, x); x = r;
    x = XT_MIN_S(x, 80.f);

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
    y = XT_MUL_S(y, BBE_MOVF32_FROM32(e1));
    z = XT_ADD_S(y, half);
    /* Initial approximation for 1/y */
    r = XT_RECIP0_S(z);
    /* 2 Newton-Raphson iterations for 1/z  */
    eps = one; XT_MSUB_S(eps, z, r);
    XT_MADD_S(r, r, eps);
    eps = one; XT_MSUB_S(eps, z, r);
    XT_MADD_S(r, r, eps);
    z = XT_SUB_S(one, r);
  }
  else
  {
    /*
    * Use polynomial approximation for small input values. This branch is
    * also used for a NaN on input.
    */

    float32_t x2, x3, tn0, tn1, tn2, tn3;
    x2 = XT_MUL_S(x, x);
    x3 = XT_MUL_S(x, x2);
    tn0 = polytanhf_tbl[0].f;
    tn1 = polytanhf_tbl[1].f;
    tn2 = polytanhf_tbl[2].f;
    tn3 = polytanhf_tbl[3].f;
    XT_MADD_S(tn1, tn0, x2);
    XT_MADD_S(tn2, tn1, x2);
    XT_MADD_S(tn3, tn2, x2);
    z = x;
    XT_MADD_S(z, tn3, x3);
  }
  /* Restore the input sign. */
  z = BBE_MOVF32_FROM32(BBE_OPERATOR_OR32(BBE_MOV32_FROMF32(z), ux));    /* apply sign */

  BBE_MOVSCFV(SCF);

  return (z);
} /* stanhf() */
#else
DISCARD_FUN(float32_t,stanhf,( float32_t x ))
#endif

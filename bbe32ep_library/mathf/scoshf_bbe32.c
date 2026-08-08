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
    Hyperbolic Cosine
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
/* cosh onstants and tables, single precision */
#include "coshf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* Constants and polynomial coeffs for exp(x) approximation. */
#include "expf_tbl.h"
#include "pow2f_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Hyperbolic Cosine

Description: These functions compute hyperbolic cosine of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vcoshf(), scoshf()
3 ULP for vfastcoshf()

Notes for non-fast versions:
1. Hyperbolic cosine functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Due to limited dynamic range of single precision floating-point format,
   hyperbolic cosine result for an input value x such that |x|>89.41599 is
   HUGE_VALF.

Input domain for 'fast' version vfastcoshf():
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
float32_t scoshf ( float32_t x )
{
  float32_t quarter;
  float32_t r;
  float32_t y, p0, dy, y1;
  xb_int32v e1, e2;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  xtfloat big = BBE_MOVF32_FROM32(0x7f7f7f7f); /* 0x7f7f7f7f, 3.4e38f */
  xtfloat inf = BBE_MOVF32_FROM32(0x7F800000);
  vbool1 binf, bbig;

  x = XT_ABS_S(x);
  if (BBE_MOVAB1(XT_UN_S(x, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, x);
    return x;
  }
  bbig = XT_OLT_S(coshf_maxarg.f, x);
  if (BBE_MOVAB1(bbig))
  {
    /* replace input with +Inf if it is too big and set ERANGE (if x is not +Inf) */
    __Pragma("frequency_hint never");
    binf = XT_OEQ_S(x, inf);
    if (vbool1_rtor_xtbool(~binf))
    {
      errno = ERANGE;
    }
    x = big;
    XT_MOVT_S(x, inf, binf);
    x = XT_ADD_S(x, x);
    return x;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  quarter = XT_CONST_S(3);
  quarter = XT_MUL_S(quarter, quarter);

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
    tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254-30 -1));
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

  r = XT_RECIP_S(y);   /* Initial approximation for 1/z */

  /* z <- z + 0.25/z == exp(x)/2 - exp(-x)/2 */
  XT_MADD_S(y, quarter, r);
  BBE_MOVSCFV(SCF);
  return y;
} /* scoshf() */
#else
DISCARD_FUN(float32_t,scoshf,( float32_t x ))
#endif

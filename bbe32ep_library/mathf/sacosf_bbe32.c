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
    Arccosine
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
/* pi, pi/2 values */
#include "pif_tbl.h"
/* asin polynomial coeff table */
#include "asinf_tbl.h"
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Arccosine

Description: These functions compute arccosine of input data. Functions output
is in radians.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy:
2 ULP - vacosf(), sacosf()
3 ULP - vfastacosf()

Notes:
1. Arccosine functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
3. Input values should belong to [-1,1], otherwise functions raise the
   "invalid" floating-point exception, assign EDOM to errno and set the
   respective output value to NaN.

Input domain for 'fast' version vfastacosf():
|x|<=1
The output value is not defined outside of this range

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
float32_t sacosf ( float32_t x )
{
  float32_t zero, one, two, half;
  float32_t g, r, s, t;
  float32_t a, b, c, y;
  int sx;
  vbool1 xnan, xinv;
  vbool1 gzero;
  int  xlesshalf;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  one = (float32_t)XT_CONST_S(1);
  half = (float32_t)XT_CONST_S(3);
  zero = (float32_t)XT_CONST_S(0);
  two = (float32_t)XT_CONST_S(2);
  sx = XT_RFR(x);
  x = XT_ABS_S(x);
  xnan = XT_UN_S(x, x);
  xinv = XT_ULT_S(one, x);
  xlesshalf = BBE_MOVAB1(XT_ULT_S(x, half));

  if (BBE_MOVAB1(xnan))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;

    x = XT_ADD_S(x, x);
    return x;
  }
  if (BBE_MOVAB1(xinv))
  {
    __Pragma("frequency_hint never");
    __feraiseexcept(FE_INVALID);
    errno = EDOM;
    return qNaNf.f;
  }

  XT_MOVT_S(x, two, xinv);

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  if (xlesshalf)
  {
    y = x;
    g = XT_MUL_S(y, y);
    /* compute r */
    r = asinftbl[0].f;
    t = asinftbl[1].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[2].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[3].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[4].f; XT_MADD_S(t, r, g); r = t;
    r = XT_MUL_S(r, g);
    XT_MADD_S(y, y, r);
    XT_MOVLTZ_S(y, -y, sx);
    y = XT_SUB_S(pi2f.f, y);
  }
  else
  {
    x = XT_SUB_S(one, x);
    g = XT_MUL_S(x, half);
    /* compute r */
    r = asinftbl[0].f;
    t = asinftbl[1].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[2].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[3].f; XT_MADD_S(t, r, g); r = t;
    t = asinftbl[4].f; XT_MADD_S(t, r, g); r = t;
    r = XT_MUL_S(r, g);

    /* Initial approximation */
    gzero = XT_OEQ_S(g, zero);
    XT_MOVT_S(g, one, gzero);  // to prevent from DIVZ exceptions if x==+/-1
    t = XT_RSQRT0_S(g);
    XT_MOVT_S(t, zero, gzero);

    /* Newton-Raphson iteration */
    a = XT_MUL_S(g, t);
    b = XT_MUL_S(half, t);
    c = one; XT_MSUB_S(c, a, t);
    XT_MADD_S(t, b, c);

    /* Newton-Raphson iteration */
    a = XT_MUL_S(g, t);
    b = XT_MUL_S(half, t);
    c = one; XT_MSUB_S(c, a, t);

    XT_MADD_S(t, b, c);
    y = XT_MUL_S(x, t);
    s = zero; XT_MOVLTZ_S(s, pif.f, sx);
    XT_MADD_S(y, y, r);
    XT_MOVLTZ_S(y, -y, sx);
    y = XT_ADD_S(s, y);
  }

  BBE_MOVSCFV(SCF);

  return y;
} /* sacosf() */
#else
DISCARD_FUN(float32_t,sacosf,( float32_t x ))
#endif

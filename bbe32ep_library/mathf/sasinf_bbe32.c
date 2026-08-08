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
    Arcsine
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
#include <__fenv.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* pi, pi/2 values */
#include "pif_tbl.h"
/* asin polynomial coeff table */
#include "asinf_tbl.h"
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Arcsine

Description: These functions compute arcsine of input data. Functions output
is in radians.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP - vasinf(), sasinf()
3 ULP - vfastasinf()

Notes for non-fast versions:
1. Arcsine functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Input values should belong to [-1,1], otherwise functions raise the
   "invalid" floating-point exception, assign EDOM to errno and set the
   respective output value to NaN.

Input domain for 'fast' version vfastasinf():
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
float32_t sasinf ( float32_t x )
{
  xtfloat g, f, z, one = 1.f, half = 0.5f, zero = 0.f;
  xtbool xbig;
  float32_t r;
  xb_int32v ux;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  if (BBE_MOVAB1(XT_UN_S(x, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, x);
    return x;
  }
  z = XT_ABS_S(x);
  if (BBE_MOVAB1(XT_OLT_S(one, z)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    __feraiseexcept(FE_INVALID);
    return qNaNf.f;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  /*
  Modyfied Cody Algorithm from Software Manual for Elementary Functions
  Uses polynomial approximation of r=asin(x)/x-1 instead of rational
  approximation (which involves division)
  */
  ux = BBE_MOV32_FROMF32(x);
  ux = BBE_OPERATOR_AND32(ux, 0x80000000);/* take sign */
//   ux >>= 31; ux <<= 31;  /* take sign */
  xbig = XT_OLE_S(half, z);
  if (xbig)
  {
    xtfloat t, a, b, c, y, dy, w, gg;
    xtbool gzero;
    z = XT_SUB_S(one, z);
    g = XT_MUL_S(z, half);
    /* compute r */
    gg = g;
    r = asinftbl[0].f;
    /* Initial approximation */
    gzero = XT_OEQ_S(g, zero);
    XT_MOVT_S(g, one, gzero);  // to prevent from DIVZ exceptions if x==+/-1
    t = XT_RSQRT0_S(g);
    XT_MOVT_S(t, zero, gzero);
    /* Newton-Raphson iteration */
    a = XT_MUL_S(g, t);
    b = XT_MUL_S(half, t);
    c = one; XT_MSUB_S(c, a, t);
    w = asinftbl[1].f; XT_MADD_S(w, r, gg); r = w;
    XT_MADD_S(t, b, c);
    /* Newton-Raphson iteration */
    a = XT_MUL_S(g, t);
    w = asinftbl[2].f; XT_MADD_S(w, r, gg); r = w;
    b = XT_MUL_S(half, t);
    c = one; XT_MSUB_S(c, a, t);
    w = asinftbl[3].f; XT_MADD_S(w, r, gg); r = w;
    XT_MADD_S(t, b, c);
    y = XT_MUL_S(z, t);
    w = asinftbl[4].f; XT_MADD_S(w, r, gg); r = w;
    dy = XT_MUL_S(y, y);
    r = XT_MUL_S(r, gg);
    XT_MSUB_S(z, half, dy); dy = z;
    XT_MADD_S(dy, y, r);
    y = XT_ADD_S(dy, y);
    z = XT_SUB_S(pi2f.f, y);
    z = BBE_MOVF32_FROM32(BBE_OPERATOR_OR32(BBE_MOV32_FROMF32(z), ux));    /* apply sign */
  }
  else
  {
    g = XT_MUL_S(z, z);
    r = asinftbl[0].f;
    f = asinftbl[1].f; XT_MADD_S(f, r, g); r = f;
    f = asinftbl[2].f; XT_MADD_S(f, r, g); r = f;
    f = asinftbl[3].f; XT_MADD_S(f, r, g); r = f;
    f = asinftbl[4].f; XT_MADD_S(f, r, g); r = f;
    r = XT_MUL_S(r, g);
    XT_MADD_S(z, z, r);
    z = BBE_MOVF32_FROM32(BBE_OPERATOR_OR32(BBE_MOV32_FROMF32(z), ux));    /* apply sign */
  }

  BBE_MOVSCFV(SCF);

  return z;

} /* sasinf() */
#else
DISCARD_FUN(float32_t,sasinf,( float32_t x ))
#endif

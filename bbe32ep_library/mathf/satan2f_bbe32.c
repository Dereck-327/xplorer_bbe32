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
    Full Arctangent (Floating-Point)
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
#include "inff_tbl.h"
/* atan polynomial coeff table */
#include "atanf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Full Arctangent (Floating-Point)

Description: These functions calculate full-quadrant arctangent of 
ratio y/x and output results in radians.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
     y    |   x   |  Result   |  Extra Conditions    
  --------+-------+-----------+---------------------
   +/-0   | -0    | +/-pi     |
   +/-0   | +0    | +/-0      |
   +/-0   |  x    | +/-pi     | x<0
   +/-0   |  x    | +/-0      | x>0
   y      | +/-0  | -pi/2     | y<0
   y      | +/-0  |  pi/2     | y>0
   +/-y   | -inf  | +/-pi     | finite y>0
   +/-y   | +inf  | +/-0      | finite y>0
   +/-inf | x     | +/-pi/2   | finite x
   +/-inf | -inf  | +/-3*pi/4 | 
   +/-inf | +inf  | +/-pi/4   |

Notes for non-fast versions:
1. Full arctangent functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Full arctangent functions assign EDOM to errno whenever x==0 and y==0.
   The resulting value depends on signs of x and y, see the Special Cases
   above.

Accuracy: 
2 ULP for vatan2f(),satan2f()
3 ULP for vfastatan2f()

Input domain for 'fast' version:
1.1755e-038 < |x| < Inf
1.1755e-038 < |y| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
y[N]    Numerator values
x[N]    Denominator values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,y,x   Aligned on 32-byte boundary
z,y,x   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t satan2f ( float32_t y, float32_t x )
{
  float32_t z;
  float32_t  zero, one, half;
  const union ufloat32uint32 * p;
  float32_t num, den, eps;
  float32_t cf0, cf1, cf2, cf3;
  float32_t x2;
  int sx, sy;
  vbool1 big, b_subn;
  float32_t a;
  int cond;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  if (BBE_MOVAB1(XT_UN_S(y, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, y);
    return x;
  }

  zero = (float32_t)XT_CONST_S(0);
  one = (float32_t)XT_CONST_S(1);
  half = (float32_t)XT_CONST_S(3);
  sx = XT_RFR(x);
  sy = XT_RFR(y);

  if (XT_ANDB(XT_OEQ_S(x, zero), XT_OEQ_S(y, zero)))
  {
    __Pragma("frequency_hint never");
    z = pif.f;
    XT_MOVLTZ_S(z, -z, sy);
    XT_MOVGEZ_S(z, y, sx);
    errno = EDOM;
    return z;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  x = XT_ABS_S(x);
  y = XT_ABS_S(y);

  big = XT_OLT_S(y, x);
  num = XT_MIN_S(x, y);
  den = XT_MAX_S(x, y);

  /* Scale data up whenever the smaller of x and y appears subnormal. This
  * avoids two potential issues:
  *   - overflow of 1/den
  *   - loss of precision in num/den when num is subnormal. */
  {
    float32_t t;
    b_subn = XT_OLT_S(num, realminf.f);
    t = zero;
    XT_MADDN_S(t, den, 8388608.f);
    XT_MOVT_S(den, t, b_subn);
    t = zero;
    XT_MADDN_S(t, num, 8388608.f);
    XT_MOVT_S(num, t, b_subn);
  }

  /* The condition is inverted to allow for NaN propagation. */
  if (den >= plusInff.f)
  {
    __Pragma("frequency_hint never");
    y = zero; XT_MADDN_S(y, y, num);    // MADDN protects from INV exception here!

    XT_MOVT_S(y, pi4f.f, XT_OLE_S(plusInff.f, num));
    XT_MOVF_S(y, XT_SUB_S(pi2f.f, y), big);
    XT_MOVLTZ_S(y, XT_SUB_S(pif.f, y), sx);
    XT_MOVLTZ_S(y, -y, sy);

    BBE_MOVSCFV(SCF);
    return y;
  }
  /* Initial approximation for 1/den */
  y = XT_RECIP0_S(den);
  /* Newton-Raphson iteration */
  eps = one;
  XT_MSUB_S(eps, den, y);
  XT_MADD_S(y, eps, y);
  /* Quotient approximation: nom/den */
  x = XT_MUL_S(y, num);
  /* Newton-Raphson iteration */
  eps = num;
  XT_MSUB_S(eps, den, x);
  XT_MADD_S(x, eps, y);

  p = (x<half) ? atanftbl1 : atanftbl2;

  /*
  * Approximate atan(x)/x-1. Use a combination of Estrin's and Horner's
  * evaluation schemes.
  */

  cf0 = p[1].f; XT_MADD_S(cf0, p[0].f, x);
  cf1 = p[3].f; XT_MADD_S(cf1, p[2].f, x);
  cf2 = p[5].f; XT_MADD_S(cf2, p[4].f, x);
  cf3 = p[7].f; XT_MADD_S(cf3, p[6].f, x);

  x2 = XT_MUL_S(x, x);
  XT_MADD_S(cf1, cf0, x2);
  XT_MADD_S(cf2, cf1, x2);
  XT_MADD_S(cf3, cf2, x2);
  y = cf3;

  /* convert result to true atan(x) */
  XT_MADD_S(x, x, y); y = x;
  a = 0.f;
  XT_MOVF_S(a, pi2f.f, big);
  XT_MOVLTZ_S(a, XT_SUB_S(pif.f, a), sx);
  XT_MOVLTZ_S(a, -a, sy);
  cond = (BBE_MOVAB1(~big) << 31) ^ sx ^ sy;
  XT_MOVLTZ_S(y, -y, cond);
  y = XT_ADD_S(y, a);

  BBE_MOVSCFV(SCF);
  return y;

} /* satan2f() */
#else
DISCARD_FUN(float32_t,satan2f,( float32_t y, float32_t x ))
#endif

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
  NatureDSP_Baseband library. Complex Math functions
    Phase Angle Of Complex Number
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(float32_t, sargf, ( complex_float x ))
#else
#include "pif_tbl.h"
#include "atanf_tbl.h"

/*-------------------------------------------------------------------------
Phase Angle Of Complex Number

Description: These functions compute phase angles of input complex numbers
and return the results in radians: [-pi,pi].

Representation: IEEE-754 Std. single precision floating-point format

Accuracy: 
2 ULP - vargf, sargf
3 ULP - vfastargf

Input domain for 'fast' version vfastargf:
1.1755e-038<|real(x)|<Inf
1.1755e-038<|imag(x)|<Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
x[N]  Complex numbers
N     Length of vectors
Output:
z[N]  Phase angles

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/

float32_t sargf ( complex_float x )
{
  const float32_t *       pTbl;

  float32_t one, zero, half, t;
  float32_t yin, xin;
  float32_t num, den, g, r, p, p2, eps;
  float32_t t0, t1;
  float32_t cf0, cf1, cf2, cf3;
  float32_t cf4, cf5, cf6, cf7;
  vbool1 bsx, bsy, bxley, beq0, bge05;
  vbool1 binf_num, binf_den, bsubn, b0, b1;

  one = XT_CONST_S(1);
  zero = XT_CONST_S(0);
  half = XT_CONST_S(3);

  /* Evaluate the argument as atan2f(x.im, x.re) */
  union
  {
    complex_float ALIGN(8) z;
    struct { float32_t re, im; }s;
  }z;

  z.z = x;

  xin = z.s.re;
  yin = z.s.im;
  /*
  * Determine the target quadrant from input signs and reduce the problem to the
  * range of (0,pi/4].
  */
  /* Check y == 0 && x == 0 condition. */
  b0 = XT_OEQ_S(yin, zero);
  b1 = XT_OEQ_S(xin, zero);
  beq0 = BBE_OPERATOR_ANDB1(b0, b1);

  {
    uint32_t ux, uy;
    ux = XT_RFR(xin);
    uy = XT_RFR(yin);
    bsx = BBE_MOVBA1(ux>>31);
    bsy = BBE_MOVBA1(uy>>31);
  }

  yin = XT_ABS_S(yin);
  xin = XT_ABS_S(xin);

  /* Check for subnomals in BOTH x and y. */
  t = (xtfloat)BBE_MOVF32_FROM32(0x007FFFFF);
  bsubn = XT_OLE_S(xin, t) | XT_OLE_S(yin, t);

  /* Substitute y=0,x=1 for all-zeros case. The actual result depends on
  * input signs. */
  XT_MOVT_S(xin, one, beq0);

  /* Determine the dominating axis. */
  bxley = XT_OLE_S(xin, yin);
  /* Select the numerator and denominator. */
  num = XT_MIN_S(xin, yin);
  den = XT_MAX_S(yin, xin);

  /* Detect 0/Inf, Inf/Inf conditions. */
  t = (xtfloat)BBE_MOVF32_FROM32(0x7F800000);
  binf_den = XT_OLE_S(t, den);
  binf_num = XT_OLE_S(t, num);

  /* Conditionally scale the numerator and denominator by 2^23 to avoid
  * overflow when inverting a subnormal number. */
  t = XT_WFR(XT_SLLI(XT_MOVI(127 + 23), 23));
  r = XT_MUL_S(num, t); XT_MOVT_S(num, r, bsubn);
  r = XT_MUL_S(den, t); XT_MOVT_S(den, r, bsubn);

  /* Detect 0/Inf, Inf/Inf conditions. */
  t = (xtfloat)BBE_MOVF32_FROM32(0x7F800000);
  binf_den = XT_OLE_S(t, den);
  binf_num = XT_OLE_S(t, num);

  /*
  * Calculate the quotient num/den.
  */
#if 1
  /* Initial appromimation for 1/den. */
  r = XT_RECIP0_S(den);
  /* Newton-Raphson iteration for 1/den. */
  eps = one;
  XT_MSUB_S(eps, den, r);
  XT_MADD_S(r, eps, r);
  /* Approximation for the quotient. */
  p = XT_MUL_S(num, r);
  /* Refine the quotient by a modified Newton-Raphson iteration. */
  eps = num;
  XT_MSUB_S(eps, den, p);
  XT_MADD_S(p, eps, r);

  r = XT_MUL_S(zero, num);
  XT_MOVT_S(r, one, binf_num);
  XT_MOVT_S(p, r, binf_den);
#else
  p = XT_DIV_S(num, den);
  r = XT_MUL_S(zero, num);
  XT_MOVT_S(r, one, binf_num);
  XT_MOVT_S(p, r, binf_den);
#endif

  /*
  * Load and select two alternative coeff sets for p<0.5 and p>=0.5
  */
  bge05 = XT_OLE_S(half, p);
  {
    int tbl;
    tbl = (int)&atanftbl2[0].f;
    XT_MOVEQZ(tbl, (int)&atanftbl1[0].f, (int)bge05);
    pTbl = (float32_t*)tbl;
  }

  cf0 = pTbl[0];
  cf1 = pTbl[1];
  cf2 = pTbl[2];
  cf3 = pTbl[3];
  cf4 = pTbl[4];
  cf5 = pTbl[5];
  cf6 = pTbl[6];
  cf7 = pTbl[7];
  /*
  * Calculate the polynomial: g <- atan(p)/p. Use a combination of Estrin's
  * and Horner's evaluation schemes.
  */
  p2 = XT_MUL_S(p, p);
#if 1
  XT_MADD_S(cf1, cf0, p);
  XT_MADD_S(cf3, cf2, p);
  XT_MADD_S(cf5, cf4, p);
  XT_MADD_S(cf7, cf6, p);

  t0 = cf1;
  t1 = cf3; XT_MADD_S(t1, t0, p2);
  t0 = cf5; XT_MADD_S(t0, t1, p2);
  t1 = cf7; XT_MADD_S(t1, t0, p2);
  g = t1;
#else
  XT_MADD_S(cf1, cf0, p);
  XT_MADD_S(cf2, cf1, p);
  XT_MADD_S(cf3, cf2, p);
  XT_MADD_S(cf4, cf3, p);
  XT_MADD_S(cf5, cf4, p);
  XT_MADD_S(cf6, cf5, p);
  XT_MADD_S(cf7, cf6, p);

  g = cf7;
#endif

  /*
  * Deduce atan(p): p <- g*p + p. Just in case, substitute pi/4 or zero for
  * Inf/Inf and 0/Inf cases, respectively.
  */
  /*  Restore the quadrant.     */

  /* |x|<=|y|: p <- pi/2 - p */
  {
    vbool1 SQ;
    float32_t A, Z;
    SQ = bxley^bsx^bsy;
    A = pi2f.f;
    XT_MOVT_S(A, zero, BBE_NOTB1(bxley));
    Z = XT_SUB_S(pif.f, A);
    XT_MOVT_S(A, Z, bsx);
    Z = XT_NEG_S(A);
    XT_MOVT_S(A, Z, bsy);

    r = XT_NEG_S(p);
    XT_MOVT_S(p, r, SQ);
    XT_MADD_S(p, p, g);
    p = XT_ADD_S(p, A);
  }
  return (float32_t)p;
} /* sargf() */

#endif/* !HAVE_VFPU */

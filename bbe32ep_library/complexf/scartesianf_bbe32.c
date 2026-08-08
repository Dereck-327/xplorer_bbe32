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
    Cartesian To Polar Conversion
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
DISCARD_FUN(float32_t, scartesianf, ( complex_float x, float32_t * phase ))
#else
#include "pif_tbl.h"
#include "atanf_tbl.h"

/*-------------------------------------------------------------------------
Cartesian To Polar Conversion

Description: These functions convert Cartesian coordinates of points on
the complex plane to the polar system (magnitude and phase).

Representation: IEEE-754 Std. single precision floating-point format

Input domain for 'fast' version vfastcartesianf():
|real(x)|>1.1755e-038
|imag(x)|>1.1755e-038
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy: 
2 ULP for vcartesianf(),scartesianf()
3 ULP for vfastcartesianf()

Parameters:
Input:
x[N]       Cartesian coordinates
N          Length of vectors
Output:
z[N]       Magnitude data
phase[N]   Phase data

Restrictions:
z,phase,x  Aligned on 32-byte boundary
z,phase,x  Must not overlap
N          Multiple of 8
-------------------------------------------------------------------------*/

float32_t scartesianf ( complex_float x, float32_t * phase )
{
  const float32_t *       pTbl;

  float32_t one, zero, half, t;
  float32_t yin, xin, X_re, X_im;
  float32_t num, den, g, r, p, p2, eps;
  float32_t t0, t1;
  float32_t cf0, cf1, cf2, cf3;
  float32_t cf4, cf5, cf6, cf7;
  vbool1 bsx, bsy, bxley, beq0, bge05;
  vbool1 binf_num, binf_den, bsubn, b0, b1;
  union
  {
    complex_float ALIGN(8) z;
    struct { float32_t re, im; }s;
  } z;

  NASSERT(phase);
  
  one = XT_CONST_S(1);
  zero = XT_CONST_S(0);
  half = XT_CONST_S(3);

  /* Evaluate the argument as atan2f(x.im, x.re) */
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

  X_im = yin = XT_ABS_S(yin);
  X_re = xin = XT_ABS_S(xin);

  /* Check for subnomals in BOTH x and y. */
  t = XT_WFR(0x007FFFFF);
  bsubn = XT_OLE_S(xin, t) | XT_OLE_S(yin, t);

  /* Substitute y=0,x=1 for all-zeros case. The actual result depends on
  * input signs. */
  XT_MOVT_S(xin, one, beq0);

  /* Determine the dominating axis. */
  bxley = XT_OLE_S(xin, yin);
  /* Select the numerator and denominator. */
  num = XT_MIN_S(xin, yin);
  den = XT_MAX_S(yin, xin);

  /* Conditionally scale the numerator and denominator by 2^23 to avoid
  * overflow when inverting a subnormal number. */
  t = XT_WFR(XT_SLLI(XT_MOVI(127 + 23), 23));
  r = XT_MUL_S(num, t); XT_MOVT_S(num, r, bsubn);
  r = XT_MUL_S(den, t); XT_MOVT_S(den, r, bsubn);

  /* Detect 0/Inf, Inf/Inf conditions. */
  t = XT_WFR(0x7F800000);
  binf_den = XT_OLE_S(t, den);
  binf_num = XT_OLE_S(t, num);

  /*
  * Calculate the quotient num/den.
  */

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
  XT_MADD_S(cf1, cf0, p);
  XT_MADD_S(cf3, cf2, p);
  XT_MADD_S(cf5, cf4, p);
  XT_MADD_S(cf7, cf6, p);

  t0 = cf1;
  t1 = cf3; XT_MADD_S(t1, t0, p2);
  t0 = cf5; XT_MADD_S(t0, t1, p2);
  t1 = cf7; XT_MADD_S(t1, t0, p2);
  g = t1;

  /*
  * Deduce atan(p): p <- g*p + p. Just in case, substitute pi/4 or zero for
  * Inf/Inf and 0/Inf cases, respectively.
  */
  /*  Restore the quadrant.     */

  /* |x|<=|y|: p <- pi/2 - p */
  {
    vbool1 SQ;
    xtfloat A, Z;
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

  *phase = (float32_t)p;

  /* cabsf() part */
  {
    xtfloat X0, Sp, Sm;
    int Exp_re, Exp_im, Exp_abs;
    Exp_re = XT_RFR(X_re)>>23;
    Exp_im = XT_RFR(X_im)>>23;
    Exp_re = Exp_re & 0xFF;
    Exp_im = Exp_im & 0xFF;
    Exp_abs= XT_MAX(Exp_re,Exp_im);
    Exp_abs= XT_MIN(Exp_abs,256-3);
    Exp_abs= XT_MAX(Exp_abs,    1);
    Sp=XT_WFR((256-2-Exp_abs)<<23);
    Sm=XT_WFR((      Exp_abs)<<23);
    /* compute normalized result */
    X_re=XT_MUL_S(X_re,Sp);
    X_im=XT_MUL_S(X_im,Sp);
    X0=XT_MUL_S(X_re,X_re);
    XT_MADD_S(X0,X_im,X_im);
    /* square root and denormalization */
    {
      float32_t z, x_red, x_adj, r, r_err, t0, t1;
      float32_t zero, half;

      zero = XT_CONST_S(0);
      half = XT_CONST_S(3);

      x_red = XT_NEXP01_S(X0);/* x with reduced exponent range */

      /* Initial rsqrt approximation with exponent range reduction */
      r = XT_SQRT0_S(X0);
      /* compute approximation error */
      t0 = XT_MUL_S(r, r);
      t1 = x_red; XT_ADDEXP_S(t1, half);/* -0.5*x */
      r_err = half; XT_MADDN_S(r_err, t0, t1);/* approximation error is (0.5-0.5*x*r*r) */
      XT_MADDN_S(r, r, r_err);/* Second recip sqrt approximation */

      /* Compute reduced range sqrt approximation */
      z = zero;
      XT_MSUBN_S(z, x_red, r);/* z = x*rsqrt(x) */

      /* Make final adjustment and restore range */
      x_adj = XT_MKSADJ_S(X0);
      XT_MADDN_S(x_red, z, z);
      XT_ADDEXPM_S(z, x_adj);
      t0 = zero;
      XT_MSUBN_S(t0, half, r);
      XT_ADDEXP_S(t0, x_adj);
      XT_DIVN_S(z, x_red, t0);
      X0 = z;
    } 
    X0 = XT_MUL_S(X0, Sm);
    return (float32_t)X0;
  } /* cabs() */
} /* scartesianf() */

#endif/* !HAVE_VFPU */

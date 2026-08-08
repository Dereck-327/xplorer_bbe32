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
    Raise To a Power
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
/* Tables */
#include "pow2f_tbl.h"
#include "sqrt2f_tbl.h"
#include "nanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Raise To a Power

Description: These functions compute the value of the first argument x 
raised to the power of the second argument y.

Data format: IEEE-754 Std. single precision floating-point.

Special cases:
      x   |   y    | Result |  Extra Conditions    
  --------+--------+--------+---------------------
    +/-0  | y      | +/-inf | odd y<0
    +/-0  | y      | +inf   | even y<0
    +/-0  | y      | +/-0   | odd y>0
    +/-0  | y      | 0      | even y>0
    +/-1  | +/-inf | 1      | 
    1     | y      | 1      | any y including NaN
    x     | +/-0   | 1      | any x including NaN
    x     | y      | NaN    | finite x<0 and finite non-integer y (see note 2)
    x     | -inf   | +inf   | |x|<1
    x     | -inf   | 0      | |x|>1
    x     | +inf   | 0      | |x|<1
    x     | +inf   | +inf   | |x|>1
    -inf  | y      | -0     | y an odd integer <0
    -inf  | y      | 0      | y<0 and not an odd integer
    -inf  | y      | -inf   | y an odd integer >0
    -inf  | y      | +inf   | y>0 and not an odd integer
    +inf  | y      | 0      | y<0
    +inf  | y      | +inf   | y>0

Accuracy: 2 ULP under condition that |y|<=100

Notes:
1. Raise to a power functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-point
   exceptions.
2. If x<0 is finite, y is finite and not an integer value, then the respective
   result z is set to NaN, errno is assigned the value EDOM, and the "invalid"
   floating-point exception is raised.

Parameters:
Input:
x[N]    Input data
y[N]    Power values
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
/* function computes log(x) in form of y0+y1 */
inline_ void approx_log2(float32_t* py0, float32_t* py1, float32_t x)
{
  vbool1 b_denorm, b_small;
  int e;
  float32_t y, y0, y1, w0, w1;

  /* take mantissa and exponent of input */
  {
    int32_t ux;
    e = -126;
    /* process denormalized inputs */
    b_denorm = XT_OLT_S(x, realminf.f);

    XT_MOVT_S(x, XT_MUL_S(x, 16777216.f), b_denorm);
    XT_MOVT(e, XT_SUB(e, 24), b_denorm);

    ux = XT_RFR(x);
    e += (ux >> 23);
    ux &= ~0xff800000;
    ux |= 0x3f000000;
    x = XT_WFR(ux);
  }

  b_small = XT_OLT_S(x, sqrt0_5f.f);
  XT_MOVT_S(x, XT_MUL_S(x, 2.f), b_small);
  XT_MOVT(e, XT_SUB(e, 1), b_small);

  x = XT_SUB_S(1.f, x);
  /* compute polynomial, last 2 coefficients are computed in extended precision */
    {   /*  alternative decomposition - 1 multiply more, but less dependency */
      
      float32_t p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
      p0 = log2f_coef[0].f;
      p1 = log2f_coef[1].f;
      p2 = log2f_coef[2].f;
      p3 = log2f_coef[3].f;
      p4 = log2f_coef[4].f;
      p5 = log2f_coef[5].f;
      p6 = log2f_coef[6].f;
      p7 = log2f_coef[7].f;
      p8 = log2f_coef[8].f;
      p9 = log2f_coef[9].f;

      #if 1
      XT_MADD_S(p1, x, p0);
      XT_MADD_S(p2, x, p1);
      XT_MADD_S(p3, x, p2);
      XT_MADD_S(p4, x, p3);
      XT_MADD_S(p5, x, p4);
      XT_MADD_S(p6, x, p5);
      XT_MADD_S(p7, x, p6);
      XT_MADD_S(p8, x, p7);
      XT_MADD_S(p9, x, p8);
      y=p9;
      #else
      {
        float32_t x2;
        x2 = XT_MUL_S(x, x);
        XT_MADD_S(p2, x2, p0);
        XT_MADD_S(p4, x2, p2);
        XT_MADD_S(p6, x2, p4);
        XT_MADD_S(p8, x2, p6);
        XT_MADD_S(p3, x2, p1);
        XT_MADD_S(p5, x2, p3);
        XT_MADD_S(p7, x2, p5);
        XT_MADD_S(p9, x2, p7);
        XT_MADD_S(p9, x, p8);
        y = p9;
      }
      #endif
    }
    y0 = XT_MUL_S(x, y);
    y1 = y0; XT_MSUB_S(y1, x, y);
    w0 = XT_ADD_S(y0, 0.5f);
    w1 = XT_SUB_S(w0, 0.5f);
    w1 = XT_SUB_S(y0, w1);
    w1 = XT_SUB_S(w1, y1);
    y0 = w0; y1 = w1;
    w0 = XT_MUL_S(x, y0);
    w1 = w0; XT_MSUB_S(w1, x, y0);
    y0 = w0;
    XT_MSUB_S(w1, y1, x); y1 = w1;
    w0 = XT_ADD_S(y0, 1.0f);
    w1 = XT_SUB_S(w0, 1.0f);
    w1 = XT_SUB_S(y0, w1);
    w1 = XT_SUB_S(w1, y1);
    y0 = w0; y1 = w1;
    x = XT_NEG_S(x);
    w0 = XT_MUL_S(y0, x);
    w1 = w0; XT_MSUB_S(w1, y0, x);
    y0 = w0;
    XT_MSUB_S(w1, y1, x); y1 = w1;
    /* multiply by log2(e) */
    w0 = XT_MUL_S(y0, log2f_coef[12].f);
    w1 = w0;
    XT_MSUB_S(w1, y0, log2f_coef[12].f);
    XT_MADD_S(w1, y1, log2f_coef[12].f);
    XT_MSUB_S(w1, y0, log2f_coef[13].f);
    y0 = w0; y1 = w1;

    /* add e */
    w0 = XT_ADD_S(y0, (float32_t)e);
    w1 = XT_SUB_S(w0, (float32_t)e);
    w1 = XT_SUB_S(y0, w1);
    y1 = XT_SUB_S(w1, y1);
    y0 = w0;
    /* return results */
    py0[0] = y0;
    py1[0] = y1;
}
/* compute 2^x, x=-0.5...0.5 */
inline_ float32_t approx_pow2(float32_t x)
{
  float32_t p0, p1, p2, p3, p4, p5, p6;
  p0 = pow2f_coef[0].f;
  p1 = pow2f_coef[1].f;
  p2 = pow2f_coef[2].f;
  p3 = pow2f_coef[3].f;
  p4 = pow2f_coef[4].f;
  p5 = pow2f_coef[5].f;
  p6 = pow2f_coef[6].f;
  /* NOTE: do not change the order of computations and way of polynomial decomposition ! */
  XT_MADD_S(p1, x, p0);
  XT_MADD_S(p2, x, p1);
  XT_MADD_S(p3, x, p2);
  XT_MADD_S(p4, x, p3);
  XT_MADD_S(p5, x, p4);
  XT_MADD_S(p6, x, p5);
  return p6;
}

float32_t spowf ( float32_t x, float32_t y )
{
  vbool1 b_zero, b_one, b_Inf, b_NaN1, b_NaN2;
  vbool1 b_sx, b_xnan, b_xinf, b_xeqz, b_xeq1;
  vbool1 b_sy, b_ynan, b_yinf, b_yeqz, b_yint;
  vbool1 b_snan, b_xyfin, b_zinf;
  vbool1 b_generic;
  xb_int32v x_i, y_i, zero_i;
  float32_t z;
  int32_t iy;
  int sx, yint, yodd;

  vbool1 b_fe_inv, b_fe_divz, b_fe_ovfl;
  vbool1 b_edom, b_erange;
  int fe_state;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  float32_t zero = XT_CONST_S(0);
  float32_t one = XT_CONST_S(1);
  float32_t half = XT_CONST_S(3);
  zero_i = BBE_MOV32_FROMF32(zero);
  x_i = BBE_MOV32_FROMF32(x);
  y_i = BBE_MOV32_FROMF32(y);

  fe_state = __fetestexcept(FE_INVALID | FE_OVERFLOW);
  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  /* Take sign of x */
  sx = XT_RFR(x);

  yint = ((float32_t)XT_FITRUNC_S(y) == y);
  iy = (int)XT_TRUNC_S(y, 0);
  XT_MOVEQZ(iy, 0, iy - MAX_INT32);
  yodd = yint & (iy & 1); /* check if y is odd */

  sx = sx & (yodd << 31);

  /* process special values first */
  {
    b_sx = XT_OLT_S(x, zero);
    b_sy = XT_OLT_S(y, zero);

    /* Check for integer y */
    b_yint = XT_OEQ_S(XT_FITRUNC_S(y), y);
    {
      xb_int32Uv hx, hy, hxy;


      hx = XT_CLSFY_S(x);
      hy = XT_CLSFY_S(y);

      hxy = BBE_OPERATOR_OR32U(hx, hy);
      hxy = BBE_OPERATOR_AND32U(hxy, 0x40);
      b_snan = BBE_OPERATOR_NEQ32U(hxy, 0);
    }
    x = XT_ABS_S(x); /* reset sign */

    /* process special numbers */
    b_yeqz = XT_OEQ_S(zero, y);                              /*  y ==0     */
    b_yinf = XT_OEQ_S(XT_ABS_S(y), plusInff.f);              /* |y|==INF   */
    b_xeqz = XT_OEQ_S(x, zero);                              /*  x==0      */
    b_xinf = XT_OEQ_S(x, plusInff.f);                        /* |x|==INF   */
    b_xeq1 = XT_OEQ_S(x, one);                               /* |x|==1     */

    b_NaN1 = BBE_ANDNOTB1(b_sx, b_yint);                     /* x<0 && y is not an integer --> z=NaN                */
    b_NaN2 = XT_UN_S(x, y);                                  /* isnan(x)||isnan(y) --> z=NaN                        */
    b_one = BBE_OPERATOR_ANDB1(b_xeq1, BBE_ORNOTB1(b_yinf, b_sx));/* |x|==1 && ( |y|==Inf || x>0 )                       */
    b_one = BBE_OPERATOR_ORB1(b_one, b_yeqz);                         /* ( |x|==1 && ( |y|==Inf || x>0 ) ) || y==0 --> z=1.0 */
    b_Inf = BBE_ANDNOTB1(b_xinf, b_sy);                     /* x==INF && y>0 --> z=INF                             */
    b_Inf = BBE_OPERATOR_ORB1(b_Inf, BBE_OPERATOR_ANDB1(b_xeqz, b_sy));      /* x==0   && y<0 --> z=INF                             */
    b_zero = BBE_ANDNOTB1(b_xeqz, b_sy);                     /* x==0   && y>0 --> z=0.0                             */
    b_zero = BBE_OPERATOR_ORB1(b_zero, BBE_OPERATOR_ANDB1(b_xinf, b_sy));      /* x==INF && y<0 --> z=0.0                             */

    /* Save special number */
    z = half;
    XT_MOVT_S(z, qNaNf.f, b_NaN1);
    XT_MOVT_S(z, zero, b_zero);
    XT_MOVT_S(z, plusInff.f, b_Inf);
    XT_MOVT_S(z, qNaNf.f, b_NaN2);
    XT_MOVT_S(z, one, b_one);
    /* Look if we have a special or generic case. */
    b_generic = XT_OEQ_S(z, half);

    {
      /* EDOM conditions:
      * A) x<0 && x is finite && y is not an integer && y is finite
      * B) x is NaN and y!=0
      * C) y is NaN and x!=1
      * We use that (A or B or C) if z==NaN. */
      b_edom = XT_UN_S(z, z);

      /* Check if x, y are finite */
      b_xyfin = BBE_NOTB1(BBE_OPERATOR_ORB1(b_xinf, b_yinf));
      b_xyfin = BBE_ANDNOTB1(b_xyfin, b_NaN2);

      /* x<0 && x is finite && y is not an integer && y is finite --> raise "invalid" exception */
      b_fe_inv = BBE_ANDNOTB1(b_NaN1, b_xeqz);
      b_fe_inv = BBE_OPERATOR_ANDB1(b_fe_inv, b_xyfin);
      /* "invalid" exception should be also raised if either input is a signalling NaN. */
      /* Check x and y for signalling NaN. */
      b_xnan = XT_UN_S(x, x);
      b_ynan = XT_UN_S(y, y);

      b_fe_inv = BBE_OPERATOR_ORB1(b_fe_inv, b_snan);
      /* x==0 && y<0 --> raise "divide-by-zero" exception */
      b_fe_divz = BBE_OPERATOR_ANDB1(b_xeqz, b_sy);
    }
  }

  /* Computation of x^y: */
  /* x^y = 2^(y*log2(x)) */
  if ((xtbool)b_generic)
  {
    float32_t x0, x1, xy0, dxy0;
    float32_t c0, c1;
    int32_t e0, e1;
    /* compute log2(x) */
    approx_log2(&x0, &x1, x);
    /* compute y*log2(x) and separate into integer and fractional parts */
    xy0 = XT_FIROUND_S(XT_MUL_S(y, x0));
    dxy0 = XT_NEG_S(xy0);
    XT_MADD_S(dxy0, y, x0);
    XT_MADD_S(dxy0, y, x1);
    dxy0 = XT_MIN_S(dxy0, 1.0f);
    dxy0 = XT_MAX_S(dxy0, -1.0f);
    /* compute 2^fract */
    z = approx_pow2(dxy0);
    /* apply integer part */
    xy0 = XT_MAX_S(xy0, -252.f);
    xy0 = XT_MIN_S(xy0, 254.f);
    e0 = (int)xy0;
    e1 = e0 >> 1;
    e0 = (e0- e1);
    c1 = XT_WFR((e1 + 127) << 23);
    c0 = XT_WFR((e0 + 127) << 23);
    z = XT_MUL_S(z, c1);
    z = XT_MUL_S(z, c0);
  }

  /* Update errno and exceptions state */
    {
      int errno_ = 0;
      int fe_inv = 0, fe_divz = 0, fe_ovfl = 0;

      BBE_MOVSCFV(SCF);

      /* (x==0)&&(y<0) || x is finite && y is finite && z is infinite --> set ERANGE */
      b_zinf = XT_OEQ_S(z, plusInff.f);
      b_erange = BBE_OPERATOR_ANDB1(b_xeqz, b_sy);
      b_erange = BBE_OPERATOR_ORB1(b_erange, BBE_OPERATOR_ANDB1(b_xyfin, b_zinf));

      errno_ = errno;
      XT_MOVT(errno_, ERANGE, b_erange);
      XT_MOVT(errno_, EDOM, b_edom);
      errno = errno_;

      /* x!=0 && x is finite && y is finite && z is infinite --> raise "overflow" exception */
      b_fe_ovfl = BBE_ANDNOTB1(b_xyfin, b_xeqz);
      b_fe_ovfl = BBE_OPERATOR_ANDB1(b_fe_ovfl, b_zinf);

      XT_MOVT(fe_inv, FE_INVALID, b_fe_inv);
      XT_MOVT(fe_divz, FE_DIVBYZERO, b_fe_divz);
      XT_MOVT(fe_ovfl, FE_OVERFLOW, b_fe_ovfl);

      __feclearexcept(FE_INVALID | FE_OVERFLOW);
      __feraiseexcept(fe_state | fe_inv | fe_divz | fe_ovfl);
    }
    /* restore sign */
    {
      xb_int32v sgn, z_i;
      sgn = BBE_MOV32_FROMN_2X32(BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(sx)));
      z_i = BBE_MOV32_FROMF32(z);
      z_i = BBE_OPERATOR_XOR32(z_i, sgn);
      z = BBE_MOVF32_FROM32(z_i);
    }

    return z;

} /* spowf() */
#else
DISCARD_FUN(float32_t,spowf,( float32_t x, float32_t y ))
#endif

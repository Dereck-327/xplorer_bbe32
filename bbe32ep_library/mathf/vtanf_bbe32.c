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
    Tangent
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
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* tan/cotan approximation polynomial coeffs. */
#include "tanf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Angular argument reduction constants. */
#include "sinf_tbl.h"
/*-------------------------------------------------------------------------
Tangent 

Description: These functions compute tangent of input data

Representation:
vtan,stan    Signed fixed-point format
             Input data are 16-bit Q15 angular values normalized by pi,
             i.e. fixed-point functions actually compute tan(pi*x).
             Output data are 32-bit Q16.15 values.
vtanf,stanf  IEEE-754 Std. single precision floating-point format for
             input/output data. Input data are treated as angular values
             specified in radians. Floating-point functions limit the
             rangw of allowable input values, see note 3.

Accuracy:
For the fixed-point functions, accuracy depends on the input value x,
as shown in the table below:
   Range of |x|    | Absolute error  | Relative error
-------------------+-----------------+----------------
 [-pi/4; pi/4]     |     1 (Q15)     |
 [pi/4; 7pi/16]    |    15 (Q15)     |    4.6e-4
 [7pi/16; 31pi/64] |   242 (Q15)     |    1.5e-3
-------------------+-----------------+----------------
2 ULP for vtanf(),stanf()
3 ULP for vfasttanf()

Notes for non-fast versions:
1. Fixed-point function result is not defined if input value x is
   +/-pi/2 (+/-8192 in Q15 normalized by pi).
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating-point functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfasttanf():
|x|<804.2477
The output value is not defined outside of this range or accuracy is 
degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vtan) or 8 (vtanf,vfasttanf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vtanf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vtanf ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0, xabs, x2, x3, x4, g, ztan, zout, jf;
  xb_vecN_2x32v ji;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2, cnt;
  xb_vecN_2xf32 cf0, cf1, cf2, tmp;
  xb_vecN_2xf32 cf3, cf4, cf5, cf6;
  xb_vecN_2x32Uv flags, flags_n, sgn, maskx;

  vboolN_2 bcot, bsx, b_nan, b_outl;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pT = (const xb_vecN_2xf32  *)polytanf_tbl;
  const xtfloat  * restrict pT1 = (const xtfloat  *)pi2fc;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)y;
  xb_vecN_2xf32  * restrict pW = (xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Table of different constants used in computations */
  static const ALIGN(32) int32_t c_tbl[] =
  {
    (int32_t)0x3f22f983,/* 2/pi */
    (int32_t)0,
  };
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    BBE_ABSN_2XF32T(x0, x0, b_nan);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);

    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x1));
    sgn = BBE_ANDN_2X32U(flags_n, maskx);
    bsx = BBE_NEQN_2X32U(sgn, BBE_ZERON_2X32U());
    /*
    * Take absolute value of the argument and remember the actual sign. Also
    * check for a domain violation.
    */

    xabs = BBE_ABSN_2XF32(x0);
    /* Check for outliers. */
    cnt = BBE_MOVN_2XF32_FROMF32(tanf_maxval); cnt = BBE_REPN_2XF32(cnt, 0);
    b_outl = BBE_OLTN_2XF32(cnt, xabs);

    /*
    * Argument reduction.
    */

    jf = BBE_MULN_2XF32(xabs, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);

    pi2fc0 = BBE_LSN_2XF32_I(pT1, 0); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(pT1, 4); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);
    pi2fc2 = BBE_LSN_2XF32_I(pT1, 8); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);

    BBE_MULSN_2XF32(xabs, jf, pi2fc0);
    BBE_MULSN_2XF32(xabs, jf, pi2fc1);
    BBE_MULSN_2XF32(xabs, jf, pi2fc2);

    xabs = BBE_MOVN_2XF32T(qNaNf.f, xabs, b_outl);
    BBE_SVN_2XF32_IP(xabs, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");

  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(xabs, pR, 2 * BBE_SIMD_WIDTH);
   /*
    * Compute the polynomial approximation g(p^2) = tan(p)/p-1. We use a combination
    * of Horner's rule and Estrin's method to evaluate the polynomial.
    */
    cf0 = BBE_LSN_2XF32_I((const xtfloat *)pT, 0 * 4); cf0 = BBE_REPN_2XF32(cf0, 0);
    cf1 = BBE_LSN_2XF32_I((const xtfloat *)pT, 1 * 4); cf1 = BBE_REPN_2XF32(cf1, 0);
    cf2 = BBE_LSN_2XF32_I((const xtfloat *)pT, 2 * 4); cf2 = BBE_REPN_2XF32(cf2, 0);
    cf3 = BBE_LSN_2XF32_I((const xtfloat *)pT, 3 * 4); cf3 = BBE_REPN_2XF32(cf3, 0);
    cf4 = BBE_LSN_2XF32_I((const xtfloat *)pT, 4 * 4); cf4 = BBE_REPN_2XF32(cf4, 0);
    cf5 = BBE_LSN_2XF32_I((const xtfloat *)pT, 5 * 4); cf5 = BBE_REPN_2XF32(cf5, 0);
    cf6 = BBE_LSN_2XF32_I((const xtfloat *)pT, 6 * 4); cf6 = BBE_REPN_2XF32(cf6, 0);
    /* Compute polynomials by Horner's method. */
    x2 = BBE_MULN_2XF32(xabs, xabs);
    x3 = BBE_MULN_2XF32(x2, xabs);
    x4 = BBE_MULN_2XF32(x2, x2);

    BBE_MULAN_2XF32(cf2, cf1, x2); cf1 = cf2;
    BBE_MULAN_2XF32(cf4, cf3, x2); cf2 = cf4;
    BBE_MULAN_2XF32(cf6, cf5, x2); cf3 = cf6;

    g = cf0;
    BBE_MULAN_2XF32(cf1, g, x4); g = cf1;
    BBE_MULAN_2XF32(cf2, g, x4); g = cf2;
    BBE_MULAN_2XF32(cf3, g, x4); g = cf3;

    tmp = BBE_LSN_2XF32_I((const xtfloat *)c_tbl, 4);
    ztan = BBE_SUBN_2XF32(xabs, tmp);
    /* Free term of the polynomial in p^2 is zero, thus we obtain the 3rd power. */
    BBE_MULAN_2XF32(ztan, g, x3);
    BBE_SVN_2XF32_IP(ztan, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(ztan, pR, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    xabs = BBE_ABSN_2XF32(x0);
    /*
    * Argument reduction.
    */
    jf = BBE_MULN_2XF32(xabs, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);
    /*
    * Conditionally reciprocate the cotan approximation.
    */
    /* Look if the input value range is odd-numbered. */

    {
      xb_vecN_2x32v tmp, _1;
      _1 = BBE_MOVN_2X32U_FROM32U(0x1); _1 = BBE_REPN_2X32(_1, 0);
      tmp = BBE_ANDN_2X32(ji, _1);
      bcot = BBE_EQN_2X32(tmp, _1);
    }

    /* For an odd-numbered range, we have actually approximated the cotan(x) = 1/tan(x). */
    zout = ztan;BBE_RECIPN_2XF32T(zout, ztan, bcot);
    /* Adjust the output sign. */
    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    BBE_NEGN_2XF32T(zout, zout, bcot);
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(zout)));
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }
  {
    xb_vecN_2xf32 v_edom, v_fe_inv;
    int fe_inv, en_edom;
    fe_inv = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan;
    vboolN_2 b_snan, b_nan, b_edom, b_fe_inv;


    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);

    b_edom = b_nan;
    b_fe_inv = b_snan;

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_fe_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_fe_inv);
    /* Retrieve EDOM state: some x<0 OR x==sNaN OR x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    /* Merge FE_INVALID state: some x<0 OR x==sNaN (the latter is detected by hardware) */
    fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);


    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };

    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise FE_INVALID (x<0 or x==sNaN) and/or FE_DIVBYZERO (x==0). */
    __feraiseexcept(fe_inv);
  }

} /* vtanf() */
#endif

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
    Cotangent
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
/* Value of pi, pi/2, etc. */
#include "pif_tbl.h" 
#include "inff_tbl.h" 
/*-------------------------------------------------------------------------
Cotangent 

Description: These functions compute cotangent of input data.

Data format: IEEE-754 Std. single precision floating-point.
Input data are treated as angular values specified in radians.

Accuracy: 
2 ULP - vcotf(), scotf()
3 ULP - vfastcotf()

Notes for non-fast versions:
1. Cotangent functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Cotangent functions functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfastcotf():
|x|<804.2477, x!=0
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
N     Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vcotf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vcotf ( float32_t * restrict y, const float32_t * restrict x, int N )
{

  int n;
  xb_vecN_2xf32 x0, xabs, x2, x3, x4, g, ztan, zout, jf;
  xb_vecN_2x32v ji;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2, cnt;
  xb_vecN_2xf32 cf0, cf1, cf2, one_f, eps;
  xb_vecN_2xf32 cf3, cf4, cf5, cf6;
  xb_vecN_2x32Uv flags, flags_n, sgn, maskx;
  vboolN_2 bcot, bsx, b_nan, b_outl, b_refine, b_ovfl, b_inf, b_fe_ovfl;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32  * restrict pT = (const xb_vecN_2xf32  *)polytanf_tbl;
  xb_vecN_2xf32  * restrict pW = (xb_vecN_2xf32  *)y;
  static const union ufloat32uint32 _2m128 = { 0x00200000 }; /* 2^-128 */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  b_fe_ovfl = BBE_XORBN_2(b_fe_ovfl, b_fe_ovfl);
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
    /* is_outlier = ( x is NaN OR x is out-of-range ) */
    b_outl = BBE_ORBN_2(b_nan, b_outl);
    /* Replace outliers with a safe value to avoid spurious exceptions. */
   // BBE_CONSTN_2XF32T(xabs, 1, b_outl);
    /*
    * Argument reduction.
    */

    jf = BBE_MULN_2XF32(xabs, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);

    pi2fc0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVVA16C(pi2fc[0].u));
    pi2fc1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVVA16C(pi2fc[1].u));
    pi2fc2 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVVA16C(pi2fc[2].u));

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
  for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
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

    ztan=xabs;
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
      bcot = BBE_NEQN_2X32(tmp, _1);
    }

    /* Newton-Raphson refinement would result in NaN whenver p is zero or
    * small enough for 1/tan(p) to overflow. Lock the refinement procedure! */
    b_refine = BBE_OLEN_2XF32T(_2m128.f, (xabs), bcot);
    b_inf = BBE_ANDNOTBN_2(bcot, b_refine);
    b_ovfl = BBE_UNEQN_2XF32T(ztan, BBE_ZERON_2XF32(), b_inf);
    b_fe_ovfl = BBE_ORBN_2(b_fe_ovfl, b_ovfl);
    /* For an odd-numbered range, we have actually approximated the cotan(x) = 1/tan(x). */
    zout = ztan; BBE_RECIP0N_2XF32T(zout, ztan, bcot);
    /* Conditionally perform two Newton-Raphson iterations for 1/tan. Use PDX_MULSN to
    * suppress undesired exceptions. */
    /* Loop schedule deteriorates if we use BBE_CONSTN_2XF32(1) instead of a literal! */
    one_f = 1.f;
    eps = one_f;
    BBE_MULSN_2XF32(eps, zout, ztan);
    BBE_MULAN_2XF32T(zout, eps, zout, b_refine);
    eps = one_f;
    BBE_MULSN_2XF32(eps, zout, ztan);
    BBE_MULAN_2XF32T(zout, eps, zout, b_refine);

    /* Adjust the sign: for odd-numbered range it must be inverted. */
    /* Adjust the output sign. */

    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    BBE_NEGN_2XF32T(zout, zout, BBE_NOTBN_2(bcot));
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(zout)));
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }

  {
    xb_vecN_2xf32 v_edom, v_fe_inv;
    int fe_inv, en_edom, fe_divz, fe_ovfl, er_erng;
    fe_inv = 0; en_edom = 0; fe_ovfl = 0; fe_divz = 0;
    xb_vecN_2x32Uv v_snan, v_nan, v_fe_ovfl;
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
    v_fe_ovfl = BBE_MOVN_2XF32T(FE_OVERFLOW, BBE_ZERON_2XF32(), b_fe_ovfl);
    /* Retrieve EDOM state: some x<0 OR x==sNaN OR x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    /* Merge FE_INVALID state: some x<0 OR x==sNaN (the latter is detected by hardware) */
    fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);
    /* BBE_RECIP0N_2XF32T may have set the FE_DIVBYZERO status flag. */
    fe_divz = __fetestexcept(FE_DIVBYZERO);
    fe_ovfl = BBE_RMAXNUMN_2XF32(v_fe_ovfl);
    /* Assert ERANGE iff we raise FE_OVERFLOW or FE_DIVBYZERO. */
    er_erng = ((0 != fe_ovfl) || (0 != fe_divz));

    if (0 != er_erng) { __Pragma("frequency_hint never"); errno = ERANGE; };
    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };

    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise the FE_INVALID and/or FE_OVERFLOW and/or FE_DIVBYZERO exception. */
    __feraiseexcept(fe_inv | fe_ovfl | fe_divz);
  }
} /* vcotf() */
#endif

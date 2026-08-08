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

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
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
DISCARD_FUN(void,vfastcotf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vfastcotf ( float32_t * restrict y, const float32_t * restrict x, int N )
{

  int n;
  xb_vecN_2xf32 x0, xabs, x2, x3, x4, g, ztan, zout, jf, tmp;
  xb_vecN_2x32v ji;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2;
  xb_vecN_2xf32 cf0, cf1, cf2;
  xb_vecN_2xf32 cf3, cf4, cf5;
  xb_vecN_2x32Uv sgn;
  vboolN_2 bcot;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xtfloat  * restrict pT1 = (const xtfloat  *)pi2fc;
  const xb_vecN_2xf32  * restrict pT = (const xb_vecN_2xf32  *)fstpolytanf_tbl;
  xb_vecN_2xf32  * restrict pY = (xb_vecN_2xf32  *)y;
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
  };
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    /*
    * Take absolute value of the argument and remember the actual sign. Also
    * check for a domain violation.
    */
    xabs = BBE_ABSN_2XF32(x0);
    /*
    * Argument reduction.
    */
    tmp = BBE_LVN_2XF32_I((const xb_vecN_2xf32 *)c_tbl, 0); tmp = BBE_REPN_2XF32(tmp, 0);
    jf = BBE_MULN_2XF32(xabs, tmp);
    jf = BBE_FIROUNDN_2XF32(jf);

    pi2fc0 = BBE_LSN_2XF32_I(pT1, 0); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(pT1, 4); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);

    BBE_MULSN_2XF32(xabs, jf, pi2fc0);
    BBE_MULSN_2XF32(xabs, jf, pi2fc1);
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
    /* Compute polynomials by Horner's method. */
    x2 = BBE_MULN_2XF32(xabs, xabs);
    x3 = BBE_MULN_2XF32(x2, xabs);
    x4 = BBE_MULN_2XF32(x2, x2);

    BBE_MULAN_2XF32(cf1, cf0, x2);
    BBE_MULAN_2XF32(cf3, cf2, x2);
    BBE_MULAN_2XF32(cf5, cf4, x2);

    g = cf1;
    BBE_MULAN_2XF32(cf3, g, x4); g = cf3;
    BBE_MULAN_2XF32(cf5, g, x4); g = cf5;
    /*x+dx*/
    ztan = (xabs);
    /* y*x+dx+x */
    BBE_MULAN_2XF32(ztan, g, x3);
    BBE_SVN_2XF32_IP(ztan, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(ztan, pR, 2 * BBE_SIMD_WIDTH);
    pi2fc2 = BBE_LSN_2XF32_I(pT1, 8); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);
    /*
    * Take absolute value of the argument and remember the actual sign.
    */
    jf = BBE_ABSN_2XF32(x0);

    /*
    * Argument reduction.
    */
    tmp = BBE_LVN_2XF32_I((const xb_vecN_2xf32 *)c_tbl, 0); tmp = BBE_REPN_2XF32(tmp, 0);
    jf = BBE_MULN_2XF32(tmp, jf);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);
    jf = BBE_MULN_2XF32(jf, pi2fc2);
    ztan = BBE_SUBN_2XF32(ztan, jf);
    /*
    * Conditionally reciprocate the cotan approximation.
    */
    /* Look if the input value range is odd-numbered. */


    /*
    * Conditionally reciprocate the cotan approximation.
    */
    {
      xb_vecN_2x32v tmp, _1;
      _1 = BBE_MOVN_2X32U_FROM32U(0x1); _1 = BBE_REPN_2X32(_1, 0);
      tmp = BBE_ANDN_2X32(ji, _1);
      bcot = BBE_NEQN_2X32(tmp, _1);
    }
    /* For an odd-numbered range, we have actually approximated the cotan(x) = 1/tan(x). */
    zout = ztan; BBE_RECIPN_2XF32T(zout, ztan, bcot);

    /* Adjust the output sign. */
    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(zout)));
    BBE_NEGN_2XF32T(zout, zout, BBE_NOTBN_2(bcot));
    BBE_SVN_2XF32_IP(zout, pY, 2 * BBE_SIMD_WIDTH);
  }

} /* vfastcotf() */
#endif

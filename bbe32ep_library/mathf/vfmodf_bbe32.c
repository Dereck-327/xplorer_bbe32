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
    Modulus
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Modulus

Description: these functions compute the floating-point remainder that
results from dividing the first argument by the second argument. The result
is less than the second argument and has the same sign as the first argument.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
    x    |   y    | Result | Extra Conditions
  -------+--------+--------+-------------------
   +/-0  |   y    |  +/0   | y!=0
    x    | +/-Inf |   x    | x != +/-Inf
  +/-Inf |   y    |  NaN   | for any y
    x    | +/-0   |  NaN   | for any x
    x    |   y    |  NaN   | |x/y|>=2^24


Accuracy: 2 ULP

Notes:
1. Modulus functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. Modulus functions limit the range of allowable input values, as follows:
   A) If |x/y|>=2^24, then the respective result z is set to NaN
   B) If x==+/-Inf and/or y==+/-0, functions set output value z to NaN, raise
      the "invalid" floating-point exception, and assign EDOM to errno.

Parameters:
Input:
x[N]    Input values
y[N]    Modulus values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vfmodf, ( float32_t * restrict z, 
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N ))
#else
void vfmodf ( float32_t * restrict z, 
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N )
{
  int n;
  xb_vecN_2x32Uv flags, flagsx, flagsy, flagsxx, flagsyx,maskx, sgn;
  xb_vecN_2xf32 x0, y0, z0, xabs, yabs;
  xb_vecN_2xf32 q0;
  xb_vecN_2x32v x0_i, y0_i, q31_max;
  static const union ufloat32uint32 _2minus23 = { 0x33800000 };    /* 2^-23 */
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
  xb_vecN_2xf32  * restrict pZ = (      xb_vecN_2xf32  *)z;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)z;
  xb_vecN_2xf32  * restrict pW = (xb_vecN_2xf32  *)z;
  __fenv_t fenv;
  vboolN_2 b_qbig, b_sx, b_sy, b_xinf, b_yinf, b_yeqz, b_ndom, b_dom;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  flagsxx = BBE_ZERON_2X32U();
  flagsyx = BBE_ZERON_2X32U();
  q31_max = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x7FFFFFFF));
  b_dom = BBE_XORBN_2(b_dom, b_dom);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

    /* Take absolute values */
    x0_i = BBE_ANDN_2X32(BBE_MOVN_2X32_FROMN_2XF32(x0), q31_max);
    y0_i = BBE_ANDN_2X32(BBE_MOVN_2X32_FROMN_2XF32(y0), q31_max);
    xabs = BBE_MOVN_2XF32_FROMN_2X32(x0_i);
    yabs = BBE_MOVN_2XF32_FROMN_2X32(y0_i);
    q0 = BBE_DIVN_2XF32(xabs, yabs);
    BBE_SVN_2XF32_IP(q0, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (const xb_vecN_2xf32  *)y;
  pR = (const xb_vecN_2xf32  *)z;
  pW = (xb_vecN_2xf32  *)z;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(q0, pR, 2 * BBE_SIMD_WIDTH);
    b_dom = BBE_UNN_2XF32(x0, x0);
    y0 = BBE_MOVN_2XF32T(x0, y0, b_dom);
    b_dom = BBE_UNN_2XF32(y0, y0);
    x0 = BBE_MOVN_2XF32T(y0, x0, b_dom);

    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x1));
    flagsx = BBE_CLSFYN_2XF32(x0);
    sgn = BBE_ANDN_2X32U(flagsx, maskx); 
    b_sx = BBE_NEQN_2X32U(sgn, BBE_ZERON_2X32U());

    flagsy = BBE_CLSFYN_2XF32(y0);
    sgn = BBE_ANDN_2X32U(flagsy, maskx); 
    b_sy = BBE_NEQN_2X32U(sgn, BBE_ZERON_2X32U());


    flagsxx = BBE_ORN_2X32U(flagsxx, flagsx);
    flagsyx = BBE_ORN_2X32U(flagsyx, flagsy);

    /* Take absolute values */
    x0_i = BBE_ANDN_2X32(BBE_MOVN_2X32_FROMN_2XF32(x0), q31_max);
    y0_i = BBE_ANDN_2X32(BBE_MOVN_2X32_FROMN_2XF32(y0), q31_max);
    xabs = BBE_MOVN_2XF32_FROMN_2X32(x0_i);
    yabs = BBE_MOVN_2XF32_FROMN_2X32(y0_i);

    b_qbig = BBE_OLEN_2XF32(yabs, BBE_MULN_2XF32(xabs, _2minus23.f));
    q0 = BBE_FIFLOORN_2XF32(q0);
    z0 = xabs; BBE_MULSN_2XF32(z0, q0, yabs);
    /* modify q0 by 1 if roundoff error in the division caused wrong result */
    BBE_SUBN_2XF32T(q0, q0, BBE_CONSTN_2XF32(1), BBE_OLTN_2XF32(z0, BBE_CONSTN_2XF32(0)));
    z0 = xabs; BBE_MULSN_2XF32(z0, q0, yabs);


    /* Check for infinities/zeros. */
    b_xinf = BBE_UEQN_2XF32(xabs, plusInff.f);
    b_yinf = BBE_OEQN_2XF32(yabs, plusInff.f);
    b_yeqz = BBE_UEQN_2XF32(yabs, BBE_CONSTN_2XF32(0));

    /* Restore the sign. */
    BBE_NEGN_2XF32T(z0, z0, b_sx);
    /* Pass input values to output whenever x is finite and y is infinite. */
    z0 = BBE_MOVN_2XF32T(x0, z0, BBE_ANDNOTBN_2(b_yinf, b_xinf));
    /* Out-of-domain: large quotient AND finite x AND non-zero y. */
    b_ndom = BBE_ANDNOTBN_2(b_qbig, BBE_ORBN_2(b_yeqz, b_xinf));
    /* Set qNaN output for out-of-domain inputs. */
    z0 = BBE_MOVN_2XF32T(qNaNf.f, z0, b_ndom);

    BBE_SVN_2XF32_IP(z0, pZ, 2 * BBE_SIMD_WIDTH);
  }

  {
    xb_vecN_2xf32 v_edom, v_fe_inv;
    int fe_inv, en_edom;
    fe_inv = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan, v_yeqz, v_xinf;
    vboolN_2 b_snan, b_nan, b_edom, b_fe_inv, b_yeqz, b_xinf;

    flags = BBE_ORN_2X32U(flagsxx,flagsyx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x2));
    v_yeqz = BBE_ANDN_2X32U(flagsyx, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x10));
    v_xinf = BBE_ANDN_2X32U(flagsxx, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_yeqz = BBE_NEQN_2X32U(v_yeqz, maskx);
    b_xinf = BBE_NEQN_2X32U(v_xinf, maskx);

    b_edom = BBE_ORBN_2(b_nan,b_yeqz);
    b_edom = BBE_ORBN_2(b_edom, b_xinf);
    b_fe_inv = BBE_ORBN_2(b_xinf, b_yeqz);
   // b_fe_inv = BBE_ANDNOTBN_2(b_fe_inv,b_nan);
    b_fe_inv = BBE_ORBN_2(b_snan, b_fe_inv);

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
} /* vfmodf() */
#endif

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
    Logarithms
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

#include <errno.h>
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
#include <fenv.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* tables */
#include "log2f_tbl.h"
#include "sqrt2f_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
#include <math.h>
/*-------------------------------------------------------------------------
Logarithms

Description: These function compute base-2, base-10 or natural logarithm of
input data.

Representation:
vlog2,vlogn,vlog10,     Signed fixed-point format
slog2,slogn,slog10      Input data are 32-bit Q16.15, results are 16-bit Q4.11.
                        Here are a few examples for the base-2 logarithm:

                          Function | Input Data Q16.15 <real> | Result Q4.11 <real>
                        -----------+--------------------------+--------------------
                           slog2   | 65536 <2.0>              | 2048 <1.0>
                           slog2   | 2147483647 <65535.99997> | 32767 <15.9995>
                           slog2   | 1 <3.052e-5>             | -30720 <-15.0>
                        -----------+--------------------------+--------------------
vlog2f,vlognf,vlog10f,  IEEE-754 Std. single precision floating-point format
slog2f,slognf,slog10f

Accuracy:
1 LSB for the fixed-point functions
2 ULP for the floating-point functions

Notes:
1. Fixed-point Functions return -32768 for a negative or zero input.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating point functions limit the range of allowable input values:
   A) If x<0, the result is set to NaN, errno is assigned the value EDOM, and
      "invalid" floating-point exception is raised
   B) If x==0, the result is set to minus infinity, errno is assigned the value 
      ERANGE, and "divide-by-zero" floating-point exception is raised

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vlog2,vlogn,vlog10) or 8 (vlog2f,vlognf,vlog10f)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vlog2f, ( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vlog2f ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0, y0, xabs, x2;
  xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
  xb_vecN_2xf32 ef_reg;
  xb_vecN_2x32v e_reg, cnt;
  vboolN_2 bsub, b_lt;
  vboolN_2 b_edom, b_fe_inv, b_fe_divz;
  vboolN_2 b_nan, b_eqz;
  xb_vecN_2x32Uv flags, flags_n;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)y;
  xb_vecN_2xf32  * restrict pW = (xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32      * restrict ptbl = (const xb_vecN_2xf32      *)log2f_tbl;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  static const ALIGN(32) int32_t tbl[] =
  {
    -126,
    (int32_t)0x4b800000, /* 2f^24 */
    24,
    (int32_t)0x007FFFFF,
    (int32_t)0x3F000000,
    (int32_t)0x3f3504f3, 0x00800000, 0 /* sqrt(0.5) */
  };
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 fx, ft;
    xb_vecN_2x32v v, t;

    cnt = BBE_LVN_2X32_I((const xb_vecN_2x32v*)tbl, 0);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    /* frexp */
    fx = x0;
    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 6));
    bsub = BBE_OLTN_2XF32(x0, ft);
    /* Multiply subnormals by 2^23 */
    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 1));
    BBE_MULN_2XF32T(fx, fx, ft, bsub);
    v = BBE_MOVN_2X32_FROMN_2XF32(fx);
    t = (BBE_REPN_2X32(cnt, 3));
    v = BBE_ANDN_2X32(v, t);
    t = (BBE_REPN_2X32(cnt, 4));
    v = BBE_ORN_2X32(v, t);
    fx = BBE_MOVN_2XF32_FROMN_2X32(v);

    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 5));
    b_lt = BBE_OLTN_2XF32(fx, ft);
    BBE_ADDEXPN_2XF32T(fx, BBE_CONSTN_2XF32(2), b_lt); ft = BBE_CONSTN_2XF32(1);
    xabs = BBE_SUBN_2XF32( BBE_CONSTN_2XF32(1),fx);
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
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    BBE_ABSN_2XF32T(x0, x0, b_nan);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);
    BBE_LVN_2XF32_IP(xabs, pR, 2 * BBE_SIMD_WIDTH);
    /* Load and replicate polynomial coefficients */
    BBE_LVN_2XF32_XP(t7, ptbl, 0);
    t0 = BBE_REPN_2XF32(t7, 0);
    t1 = BBE_REPN_2XF32(t7, 1);
    t2 = BBE_REPN_2XF32(t7, 2);
    t3 = BBE_REPN_2XF32(t7, 3);
    t4 = BBE_REPN_2XF32(t7, 4);
    t5 = BBE_REPN_2XF32(t7, 5);
    t6 = BBE_REPN_2XF32(t7, 6);
    t7 = BBE_REPN_2XF32(t7, 7);
    t8 = BBE_LSN_2XF32_I((const xtfloat *)ptbl, 4 * 8); t8 = BBE_REPN_2XF32(t8, 0);
    t9 = BBE_LSN_2XF32_I((const xtfloat *)ptbl, 4 * 9); t9 = BBE_REPN_2XF32(t9, 0);
    //
    // /* Calculate the polynomial approximation */
    x2 = BBE_MULN_2XF32(xabs, xabs);
    //
    BBE_MULAN_2XF32(t1, t0, xabs);
    BBE_MULAN_2XF32(t3, t2, xabs);
    BBE_MULAN_2XF32(t5, t4, xabs);
    BBE_MULAN_2XF32(t7, t6, xabs);
    BBE_MULAN_2XF32(t9, t8, xabs);

    y0 = t1;
    BBE_MULAN_2XF32(t3, y0, x2); y0 = t3;
    BBE_MULAN_2XF32(t5, y0, x2); y0 = t5;
    BBE_MULAN_2XF32(t7, y0, x2); y0 = t7;
    BBE_MULAN_2XF32(t9, y0, x2); y0 = t9;
    y0 = BBE_MULN_2XF32(y0, xabs);
    BBE_SVN_2XF32_IP(y0, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 fx, ft;
    xb_vecN_2x32v v, t;
    xb_vecNx16 h, hx;

    cnt = BBE_LVN_2X32_I((const xb_vecN_2x32v*)tbl, 0);
    BBE_LVN_2XF32_IP(x0, pX,  2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pR, 2 * BBE_SIMD_WIDTH);
    /* frexp */
    fx = x0;

    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 6));
    bsub = BBE_OLTN_2XF32(fx, ft);
    /* Multiply subnormals by 2^23 */
    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 1));
    BBE_MULN_2XF32T(fx, fx, ft, bsub);

    h = BBE_MOVAV16(-126);
    BBE_SUBNX16T(h, h, 24, BBE_JOINBN_2(bsub, bsub));

    v = BBE_MOVN_2X32_FROMN_2XF32(fx);
    hx = BBE_MOVNX16_FROMN_2XF32(fx);
    hx = BBE_SRLINX16(hx, 7);/* 23 - 16 */
    hx = BBE_SELNX16I(hx, hx, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
    h = BBE_ADDNX16(h, hx);

    t = (BBE_REPN_2X32(cnt, 3));
    v = BBE_ANDN_2X32(v, t);
    t = (BBE_REPN_2X32(cnt, 4));
    v = BBE_ORN_2X32(v, t);
    fx = BBE_MOVN_2XF32_FROMN_2X32(v);

    ft = BBE_MOVN_2XF32_FROMN_2X32(BBE_REPN_2X32(cnt, 5));
    b_lt = BBE_OLTN_2XF32(fx, ft);
    BBE_SUBNX16T(h, h, 1, BBE_JOINBN_2(b_lt, b_lt));

    {

      xb_vecNx16 tmp;
      tmp = BBE_SRAINX16(h, 15);
      h = BBE_SELNX16I(tmp, h, BBE_SELI_INTERLEAVE_1_LO);
      e_reg = BBE_MOVN_2X32_FROMNX16(h);
    }
    ef_reg = BBE_FLOATN_2X32(e_reg, 0);
   // BBE_MULAN_2XF32(ef_reg, y0, xabs);
    y0 = BBE_ADDN_2XF32(y0, ef_reg);
    BBE_SVN_2XF32_IP(y0, pW, 2 * BBE_SIMD_WIDTH);
  }

  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(y0, pR, 2 * BBE_SIMD_WIDTH);
    /*
    * check for outbound values (x<=0, Inf and NaN)
    */
    {
      vboolN_2 b_ltz, b_inf;
      BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);

      b_ltz = BBE_ULTN_2XF32(x0, BBE_CONSTN_2XF32(0));
      b_eqz = BBE_OEQN_2XF32(x0, BBE_CONSTN_2XF32(0));
      b_inf = BBE_OEQN_2XF32(plusInff.f, x0);
      y0 = BBE_MOVN_2XF32T(qNaNf.f, y0, b_ltz);
      y0 = BBE_MOVN_2XF32T(minusInff.f, y0, b_eqz);
      y0 = BBE_MOVN_2XF32T(plusInff.f, y0, b_inf);
    }
    BBE_SVN_2XF32_IP(y0, pW, 2 * BBE_SIMD_WIDTH);

  }

  {
    xb_vecN_2xf32 v_edom, v_fe_inv, v_fe_divz;
    int fe_inv, fe_divz, en_edom;
    xb_vecN_2x32Uv maskx;
    fe_inv = 0; fe_divz = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan, v_eqz, v_ltz;
    vboolN_2 b_snan, b_nan, b_eqz, b_ltz;

    maskx = BBE_MOVN_2X32U_FROM32U(0x40); maskx = BBE_REPN_2X32U(maskx, 0);
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROM32U(0x2); maskx = BBE_REPN_2X32U(maskx, 0);
    v_eqz = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROM32U(0x20); maskx = BBE_REPN_2X32U(maskx, 0);
    v_nan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROM32U(0x1); maskx = BBE_REPN_2X32U(maskx, 0);
    v_ltz = BBE_ANDN_2X32U(flags, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_eqz = BBE_NEQN_2X32U(v_eqz, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_ltz = BBE_NEQN_2X32U(v_ltz, maskx);

    b_edom = BBE_ORBN_2(b_ltz, b_snan);
    b_edom = BBE_ORBN_2(b_edom, b_nan);

    b_fe_inv = BBE_ORBN_2(b_ltz, b_snan);
    b_fe_divz = b_eqz;

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_fe_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_fe_inv);
    v_fe_divz = BBE_MOVN_2XF32T(FE_DIVBYZERO, BBE_ZERON_2XF32(), b_fe_divz);

    /* Retrieve EDOM state: some x<0 OR x==sNaN OR x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    /* Merge FE_INVALID state: some x<0 OR x==sNaN (the latter is detected by hardware) */
    fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);
    /* Retrieve FE_DIVBYZERO state: some x==0 */
    fe_divz = BBE_RMAXNUMN_2XF32(v_fe_divz);

    /* For logf(), the ERANGE state is equivalent to FE_DIVBYZERO floating-point exception. */
    if (0 != fe_divz) { __Pragma("frequency_hint never"); errno = ERANGE; };
    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };
    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise FE_INVALID (x<0 or x==sNaN) and/or FE_DIVBYZERO (x==0). */
    __feraiseexcept(fe_inv | fe_divz);
  }
} /* vlog2f() */
#endif

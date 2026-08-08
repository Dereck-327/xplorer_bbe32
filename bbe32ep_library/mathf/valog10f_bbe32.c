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
    Antilogarithm and Exponential
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
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Tables */
#include "pow2f_tbl.h"
#include "alog10f_tbl.h"
/*-------------------------------------------------------------------------
Antilogarithm and Exponential

Description: These functions compute base-10 or natural antilogarithm of
input data.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 2 ULP.

Note:
Antilogarithm functions conform to ANSI C requirements on standard math 
library functions in respect to treatment of errno and floating-point
exceptions.

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
DISCARD_FUN(void,valog10f,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void valog10f ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0, y0;
  xb_vecN_2xf32 c0, c1;
  xb_vecN_2xf32 dy0, zout;
  xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6,tx;
  xb_vecN_2x32v temp0, temp1;
  vboolN_2 b_edom, b_erange;
  vboolN_2 b_nan, b_max, b_snan, b_inf;
  vboolN_2 b_outb, b_inv;
  xb_vecN_2x32Uv flags, flags_n;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pT = (const xb_vecN_2xf32  *)pow2f_coef;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pW = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  b_outb = BBE_XORBN_2(b_outb, b_outb);
  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);
    /* Check input for values that are out of domain/range */
    b_max = BBE_OLEN_2XF32(alog10fminmax[1].f, x0);/* x>=38.5318 or x==NaN */
    b_inf = BBE_UEQN_2XF32(x0, plusInff.f);
    b_max = BBE_ANDNOTBN_2(b_max, b_inf);
    b_outb = BBE_ORBN_2(b_outb, b_max);

    /*
    * Multiply by 1/log10(2)
    */
    c0 = log2_10[0].f;
    c1 = log2_10[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    dy0 = BBE_NEGN_2XF32(y0);
    BBE_MULAN_2XF32(dy0, x0, c0);
    BBE_MULAN_2XF32(dy0, x0, c1);
    dy0 = BBE_MOVN_2XF32T(BBE_ZERON_2XF32(), dy0, b_inf);

    BBE_SVN_2XF32_IP(dy0, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");

  pR = (const xb_vecN_2xf32  *)y;
  pW = (      xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(dy0, pR, 2 * BBE_SIMD_WIDTH);
    /* compute 2^x */
    /* Load and replicate polynomial coefficients */
    BBE_LVN_2XF32_XP(tx, pT, 0);
    t0 = BBE_REPN_2XF32(tx, 0);
    t1 = BBE_REPN_2XF32(tx, 1);
    t2 = BBE_REPN_2XF32(tx, 2);
    t3 = BBE_REPN_2XF32(tx, 3);
    t4 = BBE_REPN_2XF32(tx, 4);
    t5 = BBE_REPN_2XF32(tx, 5);
    t6 = BBE_REPN_2XF32(tx, 6);

    {
      xb_vecN_2xf32 d2;
      d2 = BBE_MULN_2XF32(dy0, dy0);
      BBE_MULAN_2XF32(t1, dy0, t0);
      BBE_MULAN_2XF32(t3, dy0, t2);
      BBE_MULAN_2XF32(t5, dy0, t4);

      BBE_MULAN_2XF32(t3, d2, t1);
      BBE_MULAN_2XF32(t5, d2, t3);
      BBE_MULAN_2XF32(t6, dy0, t5);
      zout = t6;
    }
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (      xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pR, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    b_max = BBE_OLEN_2XF32(x0, alog10fminmax[0].f);
    /*
    * Multiply by 1/log10(2)
    */
    c0 = log2_10[0].f;
    c1 = log2_10[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    /* resulted scaling */
    y0 = BBE_MINN_2XF32(y0, 129.f);
    y0 = BBE_MAXN_2XF32(y0, -151.f);

    {
      xb_vecNx16 tmp, e1, e2;
      tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(y0,0));
      tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254));
      e1 = BBE_SRLINX16(tmp,1);
      e2 = BBE_SUBNX16(tmp, e1);
      e1 = BBE_SLLINX16(e1, 7);
      e2 = BBE_SLLINX16(e2, 7);
      e1 = BBE_SELNX16I(e1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      e2 = BBE_SELNX16I(e2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      temp0 = BBE_MOVN_2X32_FROMNX16(e2);
      temp1 = BBE_MOVN_2X32_FROMNX16(e1);
    }

    c1 = BBE_MOVN_2XF32_FROMN_2X32(temp1);
    c0 = BBE_MOVN_2XF32_FROMN_2X32(temp0);
    zout = BBE_MULN_2XF32(zout, c1);
    zout = BBE_MULN_2XF32(zout, c0);
    zout = BBE_MOVN_2XF32T(qNaNf.f, zout, b_nan);
    zout = BBE_MOVN_2XF32T(BBE_ZERON_2XF32(), zout,b_max);
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }
  /* Process errors */
  {
    xb_vecN_2xf32 v_edom, v_erange;
    int fe_ovfl, fe_inv,en_edom;
    fe_ovfl = 0; en_edom = 0; fe_inv = 0;
    xb_vecN_2x32Uv maskx, v_nan, v_snan, v_inv;

    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx); /* sNaN */
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx); /*NaN*/

    maskx = BBE_ZERON_2X32U();
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);

    b_edom = b_nan;
    b_erange = b_outb;
    b_inv = b_snan;

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_erange = BBE_MOVN_2XF32T(FE_OVERFLOW, BBE_ZERON_2XF32(), b_erange);
    v_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_inv);
    /* Retrieve EDOM state: x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    fe_ovfl = BBE_RMAXNUMN_2XF32(v_erange);
    fe_inv = BBE_RMAXNUMN_2XF32(v_inv);

    /* set errno to ERANGE if x>=38.5318 but not +INF (equivalent to FE_OVERFLOW condition) */
    if (0 != fe_ovfl) { __Pragma("frequency_hint never"); errno = ERANGE; };
    /* set errno to EDOM if x==NaN (EDOM takes precedence over ERANGE!) */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };


    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise the FE_INVALID and/or FE_OVERFLOW exception. */
    __feraiseexcept(fe_inv | fe_ovfl);
  }
} /* valog10f() */
#endif

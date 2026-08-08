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
    Hyperbolic Tangent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Table for polynomial approximation to sinh(x). */
#include "tanhf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Tables */
#include "pow2f_tbl.h"
#include "expf_tbl.h"
/*-------------------------------------------------------------------------
Hyperbolic Tangent

Description: These functions compute hyperbolic tangent of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vtanhf(), stanhf()
3 ULP for vfasttanhf()

Note for non-fast version
Hyperbolic tangent functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-
point exceptions.

Input domain for 'fast' version vfasttanhf():
|x|<Inf
The output value is not defined outside of this range.

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
DISCARD_FUN(void,vfasttanhf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vfasttanhf ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0, y0, dy0, zout, xmag, x1, x2, x3;
  xb_vecN_2xf32 one, quarter;
  xb_vecN_2x32Uv temp0, temp1, sgn;
  vboolN_2 b_inf, b_gt;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pR = (const xb_vecN_2xf32  *)y;
  xb_vecN_2xf32  * restrict pW = (xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32  * restrict pT0;
  const xb_vecN_2xf32  * restrict pT1;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  pT0 = (const xb_vecN_2xf32  *)pow2f_coef;
  pT1 = (const xb_vecN_2xf32  *)polytanhf_tbl;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 c0, c1;
    BBE_LVN_2XF32_IP(x0, pX,  2 * BBE_SIMD_WIDTH);

    x0 = BBE_ABSN_2XF32(x0);
    x0 = BBE_MULN_2XF32(x0, 2.f);
    /*
    * Multiply by 1/ln(2)
    */
    c0 = log2_e[0].f;
    c1 = log2_e[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    dy0 = BBE_NEGN_2XF32(y0);
    BBE_MULAN_2XF32(dy0, x0, c0);
    BBE_MULAN_2XF32(dy0, x0, c1);
    BBE_SVN_2XF32_IP(dy0, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, tx;
    BBE_LVN_2XF32_IP(dy0, pR, 2 * BBE_SIMD_WIDTH);
    /* compute 2^x */
    /* Load and replicate polynomial coefficients */
    BBE_LVN_2XF32_XP(tx, pT0, 0);
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
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 c0, c1, g, r;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pR, 2 * BBE_SIMD_WIDTH);
    /* resulted scaling */
    x0 = BBE_ABSN_2XF32(x0);
    x0 = BBE_MULN_2XF32(x0, 2.f);
    /*
    * Multiply by 1/ln(2)
    */
    c0 = log2_e[0].f;
    c1 = log2_e[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    y0 = BBE_MINN_2XF32(y0, 129.f);
    y0 = BBE_MAXN_2XF32(y0, -151.f);

    {
      xb_vecNx16 tmp, e1, e2;
      tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(y0, 0));
      tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254 - 30 - 1));
      e1 = BBE_SRLINX16(tmp, 1);
      e2 = BBE_SUBNX16(tmp, e1);
      e1 = BBE_SLLINX16(e1, 7);
      e2 = BBE_SLLINX16(e2, 7);
      e1 = BBE_SELNX16I(e1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      e2 = BBE_SELNX16I(e2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
      temp0 = BBE_MOVN_2X32_FROMNX16(e2);
      temp1 = BBE_MOVN_2X32_FROMNX16(e1);
    }
    g = zout;
    zout = BBE_MULN_2XF32(g, 1073741824.f);
    /*
    * Convert (y*2^(ex-30))/2 to floating-point p == exp(x)/2
    */

    c1 = BBE_MOVN_2XF32_FROMN_2X32(temp1);
    c0 = BBE_MOVN_2XF32_FROMN_2X32(temp0);
    zout = BBE_MULN_2XF32(zout, c1);
    r = BBE_CONSTN_2XF32(3);
    zout = BBE_MULN_2XF32(zout, c0);
    zout = BBE_ADDN_2XF32(r, zout);
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pR = (const xb_vecN_2xf32  *)y;
  pW = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pR)");
  __Pragma("ymemory(pW)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 t0, r;
    xb_vecN_2xf32 tn0, tn1, tn2, tn3;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pR, 2 * BBE_SIMD_WIDTH);
    r = BBE_RECIPN_2XF32(zout);

    //BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_inf = BBE_UEQN_2XF32(zout, plusInff.f);
    BBE_CONSTN_2XF32T(r, 0, b_inf);
    xmag = BBE_ABSN_2XF32(x0);
    b_gt = BBE_OLTN_2XF32(halfln3.f, xmag);
    /*
    * Compute the result: z <- 1 - 1/y == 1 - 2/(1+exp(2*x))
    */
    one = BBE_CONSTN_2XF32(1);
    y0 = BBE_SUBN_2XF32(one, r);
    t0 = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(t0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    /*
    * Compute the result: z <- p + 0.25*r == ( exp(x) - exp(-x) )/2
    */
    quarter = BBE_CONSTN_2XF32(3);
    quarter = BBE_MULN_2XF32(quarter, quarter);

    /*
    * Compute the polynomial approximation to tanh(x). This branch also provides
    * the NaN propagation.
    */

    tn0 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 0 * 4); tn0 = BBE_REPN_2XF32(tn0, 0);
    tn1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 1 * 4); tn1 = BBE_REPN_2XF32(tn1, 0);
    tn2 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 2 * 4); tn2 = BBE_REPN_2XF32(tn2, 0);
    tn3 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 3 * 4); tn3 = BBE_REPN_2XF32(tn3, 0);


    /* Preset the output vector with indirect approximation results. */
    x1 = BBE_MOVN_2XF32T(BBE_ZERON_2XF32(), xmag, b_gt);
    zout = BBE_MOVN_2XF32T(y0, x1, b_gt);

    x2 = BBE_MULN_2XF32(x1, x1);
    x3 = BBE_MULN_2XF32(x1, x2);

    BBE_MULAN_2XF32(tn1, tn0, x2);
    BBE_MULAN_2XF32(tn2, tn1, x2);
    BBE_MULAN_2XF32(tn3, tn2, x2);
    BBE_MULAN_2XF32(zout, tn3, x3);

    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(BBE_MOVN_2X32_FROMN_2XF32(zout), sgn));
    BBE_SVN_2XF32_IP(zout, pW, 2 * BBE_SIMD_WIDTH);
  }
} /* vfasttanhf() */
#endif

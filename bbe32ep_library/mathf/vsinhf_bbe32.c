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
    Hyperbolic Sine
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
/* Constants and polynomial coeffs for exp(x) approximation. */
#include "expf_tbl.h"
/* Table for polynomial approximation to sinh(x). */
#include "sinhf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Tables */
#include "pow2f_tbl.h"
#include "expf_tbl.h"
/*-------------------------------------------------------------------------
Hyperbolic Sine

Description: These functions compute hyperbolic sine of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vsinhf(), ssinhf()
3 ULP for vfastsinhf()

Notes for non-fast versions:
1. Hyperbolic sine functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Due to limited dynamic range of single precision floating-point format,
   hyperbolic sine result for an input value x such that |x|>89.41599 is
   sign(x)*HUGE_VALF.

Input domain for 'fast' version vfastsinhf():
|x|<89.41599
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
DISCARD_FUN(void,vsinhf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vsinhf ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  /*
  * Reference C code for a single input value:
  *
  * int sx;
  * float32_t x2, z, y, dy;
  * 
  * sx = takesignf(x);
  * x = sx ? -x : x;
  * / polynomial for small x /
  * x2 = x*x;
  * z = polysinhf_tbl[0].f;
  * z = x2*z + polysinhf_tbl[1].f;
  * z = x2*z + polysinhf_tbl[2].f;
  * z = x2*z;
  * z = x*z + x;
  * / direct formula for bigger x /
  * y = halfexpf(&dy, x);
  * y += dy;
  * y = y - (ldexpf)(1.f / y, -2);
  * 
  * 
  * y = (x>1.f) ? y : z;  / if not > (or nan) use z /
  * y = sx ? -y : y;  / apply sign /
  * return y;  
  *
  * static  float32_t halfexpf(float32_t* dy, float32_t x)
  * {
  *   float32_t x0, x1, xy0, dxy0, z;
  *   union ufloat32uint32 c0, c1;
  *   int32_t e0, e1, d, y, e, t;
  *   union ufloat32uint32 xx;
  *   / log2(n) /
  *   c0.u = 0x3fb8aa3b;
  *   c1.u = 0x32a57060;
  * 
  *   x0 = c0.f;
  *   x1 = c1.f;
  * 
  *   / scale input to 1/ln(2) /
  * 
  *   xy0 = roundf(x*x0);
  *   dxy0 = fmaf(x, x0, -xy0);
  *   dxy0 = fmaf(x, x1, dxy0);
  *   if (isinf(xy0)) dxy0 = 0.f;
  *   dxy0 = MIN(0.50f, dxy0);
  *   dxy0 = MAX(-0.50f, dxy0);
  *   z = approx_pow2(dxy0);
  *   / resulted scaling /
  *   xy0 = MAX(xy0, -151.f);
  *   xy0 = MIN(xy0, 129.f);
  *   e0 = (int)xy0 - 1;
  *   e1 = e0 >> 1;
  *   e0 -= e1;
  *   c1.u = (e1 + 127 - 15) << 23;
  *   c0.u = (e0 + 127 - 15) << 23;
  *   y = (int32_t)(ldexpf)(z, 30);
  * 
  *   / decompose result by 2 floating point values /
  *   d = y & 255;
  *   y -= d;
  *   z = (float32_t)y;
  *   z = z*c1.f;
  *   z = z*c0.f;
  *   x = z;
  *   z = (float32_t)d;
  *   z = z*c1.f;
  *   z = z*c0.f;
  *   dy[0] = z;
  *   return x;  */
  int n;
  xb_vecN_2x32Uv flags, flags_n, sgn;
  xb_vecN_2xf32 x0, y0, dy0, zout, xmag, x1, x2, x3;
  xb_vecN_2xf32 one, quarter;
  xb_vecN_2x32Uv temp0, temp1;
  vboolN_2 b_nan, b_snan, b_inf, b_le1, b_outb, b_olt;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  xb_vecN_2xf32  * restrict pZ = (xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32  * restrict pT0;
  const xb_vecN_2xf32  * restrict pT1;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  b_outb = BBE_XORBN_2(b_outb, b_outb);
  pT0 = (const xb_vecN_2xf32  *)pow2f_coef;
  pT1 = (const xb_vecN_2xf32  *)polysinhf_tbl;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 c0, c1;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);
    x0 = BBE_ABSN_2XF32(x0);
    /*
    * Multiply by 1/ln(2)
    */
    c0 = log2_e[0].f;
    c1 = log2_e[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    b_inf = BBE_OEQN_2XF32(y0, plusInff.f);
    dy0 = BBE_NEGN_2XF32(y0);
    BBE_MULAN_2XF32(dy0, x0, c0);
    BBE_MULAN_2XF32(dy0, x0, c1);
    dy0 = BBE_MOVN_2XF32T(BBE_ZERON_2XF32(),dy0,b_inf);
    BBE_SVN_2XF32_IP(dy0, pY, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  pZ = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, tx;
    BBE_LVN_2XF32_IP(dy0, pY, 2 * BBE_SIMD_WIDTH);
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
    BBE_SVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  pZ = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 c0, c1,g;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);

    x0 = BBE_ABSN_2XF32(x0);
    /*
    * Multiply by 1/ln(2)
    */
    c0 = log2_e[0].f;
    c1 = log2_e[1].f;
    y0 = BBE_MULN_2XF32(x0, c0);
    y0 = BBE_FIROUNDN_2XF32(y0);
    /* resulted scaling */
    y0 = BBE_MINN_2XF32(y0, 129.f);
    y0 = BBE_MAXN_2XF32(y0, -151.f);

    {
      xb_vecNx16 tmp, e1, e2;
      tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(y0, 0));
      tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254-30-1));
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
    zout = BBE_MULN_2XF32(zout, c0);

    BBE_SVN_2XF32_IP(zout, pY, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  pZ = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 t0, r;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pY, 2 * BBE_SIMD_WIDTH);
    xmag = BBE_ABSN_2XF32(x0);
    t0 = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(t0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    r = BBE_RECIPN_2XF32(zout);

    b_inf = BBE_UEQN_2XF32(zout, plusInff.f);
    BBE_CONSTN_2XF32T(r, 0, b_inf);
    /*
    * Compute the result: z <- p + 0.25*r == ( exp(x) - exp(-x) )/2
    */
    quarter = BBE_CONSTN_2XF32(3);
    quarter = BBE_MULN_2XF32(quarter, quarter);
    BBE_MULSN_2XF32(zout, r, quarter);

    /*
    * Perform additional analysis of input data for Error Handling, and
    * saturate resulting values.
    */
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(BBE_MOVN_2X32_FROMN_2XF32(zout), sgn));
    BBE_SVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  pZ = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 sn0, sn1, sn2, t0;
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);

    t0 = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(t0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));

    xmag = BBE_ABSN_2XF32(x0);
    one = BBE_CONSTN_2XF32(1);
    b_le1 = BBE_ULEQN_2XF32(xmag, one);
    /*
    * Compute the polynomial approximation to sinh(x). There is no
    * need for explicit sign processing. Direct approximation also provides the
    * NaN propagation.
    */

    sn0 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 0 * 4); sn0 = BBE_REPN_2XF32(sn0, 0);
    sn1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 1 * 4); sn1 = BBE_REPN_2XF32(sn1, 0);
    sn2 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 2 * 4); sn2 = BBE_REPN_2XF32(sn2, 0);

    x1 = x0;
    /* Preset the output vector with indirect approximation results. */
    zout = BBE_MOVN_2XF32T(x0, zout, b_le1);

    x2 = BBE_MULN_2XF32(x1, x1);
    x3 = BBE_MULN_2XF32(x1, x2);

    BBE_MULAN_2XF32(sn1, sn0, x2);
    BBE_MULAN_2XF32(sn2, sn1, x2);
    BBE_MULAN_2XF32T(zout, sn2, x3, b_le1);

    b_inf = BBE_OEQN_2XF32(xmag, plusInff.f);
    b_olt = BBE_OLTN_2XF32T(sinhf_maxarg.f, xmag, BBE_NOTBN_2(b_inf));
    b_outb = BBE_ORBN_2(b_outb, b_olt); /* Exact result overflows */

    t0 = plusInff.f;
    t0 = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(BBE_MOVN_2X32_FROMN_2XF32(t0), sgn));

    zout = BBE_MOVN_2XF32T(t0, zout, BBE_ORBN_2(b_inf, b_olt));
    BBE_SVN_2XF32_IP(zout, pY, 2 * BBE_SIMD_WIDTH);
  }
  /* Process errors */
  {
    xb_vecN_2xf32 v_edom, v_erange;
    int fe_ovfl, fe_inv, en_edom;
    fe_ovfl = 0; fe_inv = 0; en_edom = 0;
    xb_vecN_2x32Uv maskx, v_nan, v_snan, v_inv;
    vboolN_2 b_edom, b_erange, b_inv;
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx); /* sNaN */
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx); /*NaN*/


    maskx = BBE_ZERON_2X32U();
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);


    b_edom = BBE_ORBN_2( b_nan, b_snan);
    b_erange = b_outb;
    b_inv = b_snan;

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_erange = BBE_MOVN_2XF32T(FE_OVERFLOW, BBE_ZERON_2XF32(), b_erange);
    v_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_inv);
    /* Retrieve EDOM state: x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    fe_ovfl = BBE_RMAXNUMN_2XF32(v_erange);
    fe_inv = BBE_RMAXNUMN_2XF32(v_inv);

    if (0 != fe_ovfl) { __Pragma("frequency_hint never"); errno = ERANGE; };
    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };


    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise the FE_INVALID and/or FE_OVERFLOW exception. */
    __feraiseexcept(fe_inv | fe_ovfl);
  }
} /* vsinhf() */
#endif

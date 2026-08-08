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
    Arcsine
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* pi, pi/2 values */
#include "pif_tbl.h"
/* asin polynomial coeff table */
#include "asinf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Arcsine

Description: These functions compute arcsine of input data. Functions output
is in radians.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP - vasinf(), sasinf()
3 ULP - vfastasinf()

Notes for non-fast versions:
1. Arcsine functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Input values should belong to [-1,1], otherwise functions raise the
   "invalid" floating-point exception, assign EDOM to errno and set the
   respective output value to NaN.

Input domain for 'fast' version vfastasinf():
|x|<=1
The output value is not defined outside of this range

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
DISCARD_FUN(void,vfastasinf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vfastasinf ( float32_t * restrict y, const float32_t * restrict x, int N )
{
  /*
  * Modyfied Cody Algorithm from Software Manual for Elementary Functions
  * Uses polynomial approximation of r=asin(x)/x-1 instead of rational
  * approximation (which involves division)
  *
  * Reference C code:
  *
  *   float32_t g,r,s,y,dy;
  *   int sx;
  *   sx = takesignf( x );
  *   x=sx ? -x:x;
  *   if ( x>1.f )
  *   {
  *     return (0);
  *   }
  *
  *   if(x<0.5f)
  *   {
  *     y=x;
  *     dy=0;
  *     g=y*y;
  *     s=0.f;
  *   }
  *   else
  *   {
  *     g=(1.f-x)*0.5f;
  *     y=-2.f*sqrtf(g);
  *     dy = 0.5f*y*y - (1.f-x);
  *     s=pi2f.f;
  *   }
  *
  *    compute r
  *   r=    asinftbl[0].f;
  *   r=r*g+asinftbl[1].f;
  *   r=r*g+asinftbl[2].f;
  *   r=r*g+asinftbl[3].f;
  *   r=r*g+asinftbl[4].f;
  *   r=r*g;
  *   y+=(dy+y*r);
  *   y+=s;
  *   y=sx ? -y:y;
  *   return y;
  */
  int n;
  xb_vecN_2x32Uv sgn;
  xb_vecN_2xf32 _2g, g, r, y0, dy, zout;
  xb_vecN_2xf32 t0, t1;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pT = (xb_vecN_2xf32  *)asinftbl;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pZ = (xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  xb_vecN_2xf32 x0, one, half, zero;
  vboolN_2 b_ge05, b_neq1;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4;
  if (N <= 0) return;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    /*
    * Classify input values: a) 0<=x<0.5 or b) 0.5<=x<1.0 or c) x==1.0 or d) x>1.0.
    */
    x0 = BBE_ABSN_2XF32(x0);

    one = BBE_CONSTN_2XF32(1);
    half = BBE_CONSTN_2XF32(3);
    b_ge05 = BBE_OLEN_2XF32(half, x0);
    b_neq1 = BBE_UNEQN_2XF32(one, x0);

    /*
    * Range reduction and computation of polynomial argument. We use the following
    * identity to reduce the input range to [0,0.5): asin(x) == pi/2 - 2*asin(((1-x)/2)^0.5)
    *   a) 0.0<=x<0.5: y <- x; g <= x^2
    *   b) 0.5<=x<1.0: y <- -2*((1-x)/2)^0.5; g <- (1-x)/2
    */
    _2g = BBE_SUBN_2XF32(one, x0);
    /* ADDEXP would kill NaN. */
    g = BBE_MULN_2XF32(_2g, half);
    /* 1st: we rely on PDX_RSQRT behavior for negative input to
    *      obtain FE_INVALID when |x|>1
    * 2nd: Use a conditional to avoid 1/sqrt(0) when |x|==1 */
    {
      xb_vecN_2xf32 one, half, r0, r0_err,z0;
      xb_vecN_2xf32 x0_adj, x0_red;
      half = BBE_CONSTN_2XF32(3);
      one = BBE_CONSTN_2XF32(1);
      x0_red = BBE_NEXP01N_2XF32(g);/* negated x with reduced exponent range */

      /* Initial rsqrt approximation with exponent range reduction */
      r0 = BBE_SQRT0N_2XF32(g);
      r0 = BBE_NEGN_2XF32(r0);
      /* compute approximation error */
      t0 = BBE_MULN_2XF32(r0, r0);
      t1 = x0_red; BBE_ADDEXPN_2XF32(t1, half);/* -0.5*x */
      r0_err = half; BBE_MULANN_2XF32(r0_err, t0, t1);/* approximation error is (0.5-0.5*x*r*r) */
      BBE_MULANN_2XF32(r0, r0, r0_err);/* Second recip sqrt approximation */

      /* Compute reduced range sqrt approximation */
      z0 = BBE_MULN_2XF32(x0_red, r0);/* z = x*rsqrt(x) */

      /* Make final adjustment and restore range */
      x0_adj = BBE_MKSADJN_2XF32(g);
      t0 = x0_red;
      BBE_MULANN_2XF32(t0, z0, z0);
      BBE_ADDEXPMN_2XF32(z0, x0_adj);
      t1 = BBE_MULN_2XF32(half, r0);
      BBE_ADDEXPN_2XF32(t1, x0_adj);
      BBE_DIVNN_2XF32(z0, t0, t1);
      r=z0;
    }

    BBE_SVN_2XF32_IP(r, pY, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(r, pY, 2 * BBE_SIMD_WIDTH);
    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));

    x0 = BBE_ABSN_2XF32(x0);
    one = BBE_CONSTN_2XF32(1);
    half = BBE_CONSTN_2XF32(3);
    b_ge05 = BBE_OLEN_2XF32(half, x0);
    b_neq1 = BBE_UNEQN_2XF32(one, x0);

    _2g = BBE_SUBN_2XF32(one, x0);/*1-x*/
    g = BBE_MULN_2XF32(_2g, half); /* (1-x)*0.5 */
        /* y <- ( x>=0.5f ? -2*sqrt((1-x)/2) : x ) */
    t1 = BBE_ADDN_2XF32(r, r);
    y0 = x0; BBE_NEGN_2XF32T(y0, t1, b_ge05);
    /* dy <- ( x>=0.5f ? (1.f-x) - (0.5f*y*y) : 0.f ) */
    BBE_MULSN_2XF32(_2g, t1, r);
    zero = BBE_CONSTN_2XF32(0);
    dy = BBE_MOVN_2XF32T(_2g, zero, b_ge05);

    t0 = BBE_MULN_2XF32(x0, x0);
    g = BBE_MOVN_2XF32T(g, t0, b_ge05);
    /*
    * Compute the polynomial approximation for r = asin(y)/y - 1.
    * Odd-numbered coeffs and the free term are equal to zero.
    */
    cf0 = BBE_LSN_2XF32_I((const xtfloat *)pT, 0 * 4); cf0 = BBE_REPN_2XF32(cf0, 0);
    cf1 = BBE_LSN_2XF32_I((const xtfloat *)pT, 1 * 4); cf1 = BBE_REPN_2XF32(cf1, 0);
    cf2 = BBE_LSN_2XF32_I((const xtfloat *)pT, 2 * 4); cf2 = BBE_REPN_2XF32(cf2, 0);
    cf3 = BBE_LSN_2XF32_I((const xtfloat *)pT, 3 * 4); cf3 = BBE_REPN_2XF32(cf3, 0);
    cf4 = BBE_LSN_2XF32_I((const xtfloat *)pT, 4 * 4); cf4 = BBE_REPN_2XF32(cf4, 0);
                                r = cf0;
    BBE_MULAN_2XF32(cf1, r, g); r = cf1;
    BBE_MULAN_2XF32(cf2, r, g); r = cf2;
    BBE_MULAN_2XF32(cf3, r, g); r = cf3;
    BBE_MULAN_2XF32(cf4, r, g); r = cf4;
    r = BBE_MULN_2XF32(r, g);

    /*
    * z <- y + ( r*y - dy )
    */

    BBE_MULSN_2XF32(dy, r, y0);
    zout = BBE_SUBN_2XF32(y0, dy);
    /*
    * Restore the input range and sign.
    */
    BBE_ADDN_2XF32T(zout, zout, pi2f.f, b_ge05);
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(zout)));

    BBE_SVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* vfastasinf() */
#endif

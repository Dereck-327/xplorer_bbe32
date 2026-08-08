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
    Full Arctangent (Floating-Point)
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
/* atan polynomial coeff table */
#include "atanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Full Arctangent (Floating-Point)

Description: These functions calculate full-quadrant arctangent of 
ratio y/x and output results in radians.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
     y    |   x   |  Result   |  Extra Conditions    
  --------+-------+-----------+---------------------
   +/-0   | -0    | +/-pi     |
   +/-0   | +0    | +/-0      |
   +/-0   |  x    | +/-pi     | x<0
   +/-0   |  x    | +/-0      | x>0
   y      | +/-0  | -pi/2     | y<0
   y      | +/-0  |  pi/2     | y>0
   +/-y   | -inf  | +/-pi     | finite y>0
   +/-y   | +inf  | +/-0      | finite y>0
   +/-inf | x     | +/-pi/2   | finite x
   +/-inf | -inf  | +/-3*pi/4 | 
   +/-inf | +inf  | +/-pi/4   |

Notes for non-fast versions:
1. Full arctangent functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Full arctangent functions assign EDOM to errno whenever x==0 and y==0.
   The resulting value depends on signs of x and y, see the Special Cases
   above.

Accuracy: 
2 ULP for vatan2f(),satan2f()
3 ULP for vfastatan2f()

Input domain for 'fast' version:
1.1755e-038 < |x| < Inf
1.1755e-038 < |y| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
y[N]    Numerator values
x[N]    Denominator values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,y,x   Aligned on 32-byte boundary
z,y,x   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vfastatan2f, ( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N ))
#else
void vfastatan2f ( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N )
{
  int n;
  xb_vecN_2x32Uv sgn;
  xb_vecN_2xf32 g, r, y0, zout;
  xb_vecN_2xf32 den, num, eps, p;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
  const xtfloat        * restrict pTBL;
  xb_vecN_2xf32  * restrict pZ0 = (xb_vecN_2xf32  *)z;
  xb_vecN_2xf32  * restrict pZ1 = (xb_vecN_2xf32  *)z;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  xb_vecN_2xf32 x0, xmag, ymag, one, half, zero;
  vboolN_2 b_sx, b_xley;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7;
  if (N <= 0) return;

  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

    one = BBE_CONSTN_2XF32(1);
    /*
    * Determine the target quadrant from input signs and reduce the problem to the
    * range of (0,pi/4].
    */
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    b_xley = BBE_OLEN_2XF32(xmag, ymag);
    /* Select the numerator and denominator. */
    num = BBE_MINN_2XF32(xmag, ymag);
    den = BBE_MAXN_2XF32(ymag, xmag);
    /*
    * Calculate the quotient num/den.
    */
    /* Initial appromimation for 1/den. */
    r = BBE_RECIP0N_2XF32(den);
    /* Newton-Raphson iteration. */
    eps = one;
    BBE_MULSN_2XF32(eps, den, r);
    BBE_MULAN_2XF32(r, eps, r);
    /* Approximation for the quotient. */
    p = BBE_MULN_2XF32(num, r);
    /* Refine the quotient by a modified Newton-Raphson iteration. */
    eps = num;
    BBE_MULSN_2XF32(eps, den, p);
    BBE_MULAN_2XF32(p, eps, r);
    BBE_SVN_2XF32_IP(p, pZ0, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (const xb_vecN_2xf32  *)y;
  pZ0 = (xb_vecN_2xf32  *)z;
  pTBL = (xtfloat       *)atanftbl_10ord;
  __Pragma("ymemory(pZ0)");
  __Pragma("ymemory(pZ1)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    xb_vecN_2xf32 _80000000;
    _80000000 = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000));
    zero = BBE_CONSTN_2XF32(0);
    one = BBE_CONSTN_2XF32(1);
    half = BBE_CONSTN_2XF32(3);

    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(p, pZ0, 2 * BBE_SIMD_WIDTH);

    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, _80000000);
    b_sx = BBE_EQN_2X32(sgn, _80000000);
    sgn = BBE_MOVN_2X32_FROMN_2XF32(y0);
    sgn = BBE_ANDN_2X32(sgn, _80000000);
    /*
    * Determine the target quadrant from input signs and reduce the problem to the
    * range of (0,pi/4].
    */
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    b_xley = BBE_OLEN_2XF32(xmag, ymag);

    {
      xb_vecN_2xf32 cf8, cf9;
      /*
      * Load coeff set
      */

      cf0 = BBE_LSN_2XF32_I(pTBL, 0 * sizeof(float32_t)); cf0 = BBE_REPN_2XF32(cf0, 0);
      cf1 = BBE_LSN_2XF32_I(pTBL, 1 * sizeof(float32_t)); cf1 = BBE_REPN_2XF32(cf1, 0);
      cf2 = BBE_LSN_2XF32_I(pTBL, 2 * sizeof(float32_t)); cf2 = BBE_REPN_2XF32(cf2, 0);
      cf3 = BBE_LSN_2XF32_I(pTBL, 3 * sizeof(float32_t)); cf3 = BBE_REPN_2XF32(cf3, 0);
      cf4 = BBE_LSN_2XF32_I(pTBL, 4 * sizeof(float32_t)); cf4 = BBE_REPN_2XF32(cf4, 0);
      cf5 = BBE_LSN_2XF32_I(pTBL, 5 * sizeof(float32_t)); cf5 = BBE_REPN_2XF32(cf5, 0);
      cf6 = BBE_LSN_2XF32_I(pTBL, 6 * sizeof(float32_t)); cf6 = BBE_REPN_2XF32(cf6, 0);
      cf7 = BBE_LSN_2XF32_I(pTBL, 7 * sizeof(float32_t)); cf7 = BBE_REPN_2XF32(cf7, 0);
      cf8 = BBE_LSN_2XF32_I(pTBL, 8 * sizeof(float32_t)); cf8 = BBE_REPN_2XF32(cf8, 0);
      cf9 = BBE_LSN_2XF32_I(pTBL, 9 * sizeof(float32_t)); cf9 = BBE_REPN_2XF32(cf9, 0);
      /* Use Horner's evaluation scheme */
      BBE_MULAN_2XF32(cf1, cf0, p);
      BBE_MULAN_2XF32(cf2, cf1, p);
      BBE_MULAN_2XF32(cf3, cf2, p);
      BBE_MULAN_2XF32(cf4, cf3, p);
      BBE_MULAN_2XF32(cf5, cf4, p);
      BBE_MULAN_2XF32(cf6, cf5, p);
      BBE_MULAN_2XF32(cf7, cf6, p);
      BBE_MULAN_2XF32(cf8, cf7, p);
      BBE_MULAN_2XF32(cf9, cf8, p);
      g = cf9;

      /*
      * Deduce atan(p): p <- g*p + p.
      */
    }

    BBE_MULAN_2XF32(p, g, p);
    /*
    * Restore the quadrant.
    */
    /* |x|<=|y|: p <- pi/2 - p */
    BBE_SUBN_2XF32T(p, pi2f.f, p, b_xley);
    /* x<0: z <- pi - p */
    BBE_SUBN_2XF32T(p, pif.f, p, b_sx);
    /* y<0: p <- -p */
    //BBE_NEGN_2XF32T(p, p, b_sy);
    p = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(p)));
    zout = p;
    BBE_SVN_2XF32_IP(zout, pZ1, 2 * BBE_SIMD_WIDTH);

  }

} /* vfastatan2f() */
#endif

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
    Arctangent
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
/* Quiet NaN (single precision). */
#include "nanf_tbl.h"
/* +/-Infinity (single precision). */
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Arctangent 

Description: These functions compute the principal value of arctangent.

Representation:
vatan16,satan16  16-bit signed fixed-point format
                 Input data are Q15. Functions compute atan(x)/(pi/4)
                 and output results in Q15 format.
vatanf,satanf    IEEE-754 Std. single precision floating-point format
vfastatanf       Functions compute atan(x) and output results in radians

Special cases:
    Input | Result 
   -------+--------
    +inf  |  pi/2  (floating-point functions)
    -inf  | -pi/2  

Accuracy:
1 LSB for fixed point fixed point functions
1 ULP for vatanf(), satanf()
2 ULP for vfastatanf()

Notes:
1. These functions are much faster than full-quadrant arctangent atan2,
   so they are preferable when the full phase is not required.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.

Input domain for 'fast' version vfastatanf():
|x|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vatan16) or 8 (vatanf,vfastatanf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vfastatanf, ( float32_t * restrict z, 
         const float32_t * restrict x, 
         int N ))
#else
void vfastatanf ( float32_t * restrict z, 
         const float32_t * restrict x, 
         int N )
{
  /*
  *  float32_t y;
  *  int sx, big;
  *  const union ufloat32uint32* p;
  *  / range reduction /
  *  sx = x<0;
  *  x = sx ? -x : x;
  *  big = x>1.0f;
  *  if (big) x = 1.0f / x;
  *  p = (x<0.5f) ? atanftbl1 : atanftbl2;
  *  / approximate atan(x)/x-1 /
  *  y = p[0].f;
  *  y = x*y + p[1].f;
  *  y = x*y + p[2].f;
  *  y = x*y + p[3].f;
  *  y = x*y + p[4].f;
  *  y = x*y + p[5].f;
  *  y = x*y + p[6].f;
  *  y = x*y + p[7].f;
  *  / convert result to true atan(x) /
  *  y = x*y + x;
  *
  *  if (big) y = pi2f.f - y;
  *  y = sx ? -y : y; / apply original sign /
  *  return y;
  */
  int n;
  xb_vecN_2x32v sgn;
  xb_vecN_2xf32 g, r, y0, y2, zout;
  xb_vecN_2xf32 t0, t1;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pT0 = (xb_vecN_2xf32  *)atanftbl1;
  const xb_vecN_2xf32  * restrict pT1 = (xb_vecN_2xf32  *)atanftbl2;
  xb_vecN_2xf32  * restrict pZ = (xb_vecN_2xf32  *)z;
  xb_vecN_2xf32  * restrict pY = (xb_vecN_2xf32  *)z;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  xb_vecN_2xf32 x0, xmag, one, half, two;
  vboolN_2 b_ge1, b_ge05;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7, cfx0, cfx1;
  if (N <= 0) return;
  one = BBE_CONSTN_2XF32(1);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);

    
    xmag = BBE_ABSN_2XF32(x0);
    b_ge1 = BBE_OLTN_2XF32(one, xmag);

    /*
    * Range reduction:
    *   a) 0.0<=x<1.0: y <- x
    *   b) 1.0<x    : y <- 1/x
    */
    /* Initial approximation for reciprocal value. */
    r = BBE_RECIPN_2XF32(xmag);

    y0 = BBE_MOVN_2XF32T(r, xmag, b_ge1);
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)z;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    one = BBE_CONSTN_2XF32(1);
    two = BBE_CONSTN_2XF32(2);
    half = BBE_CONSTN_2XF32(3);

    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000)));
    xmag = BBE_ABSN_2XF32(x0);
    b_ge1 = BBE_OLTN_2XF32(one, xmag);
    /*
    * Load and select two alternative coeff sets for y<0.5 and y>=0.5
    */ 
    b_ge05 = BBE_OLEN_2XF32(half, y0);
    x0 = y0;
    #if 1
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 0*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 0*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf0 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 1*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 1*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf1 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 2*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 2*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf2 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 3*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 3*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf3 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 4*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 4*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf4 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 5*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 5*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf5 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 6*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 6*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf6 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 7*4); cfx0 = BBE_REPN_2XF32(cfx0, 0); 
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 7*4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf7 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);

    #else
    cfx0 = BBE_LVN_2XF32_I(pT0, 0);
    cfx1 = BBE_LVN_2XF32_I(pT1, 0);

    cf0 = BBE_REPN_2XF32(cfx0, 0); t0 = BBE_REPN_2XF32(cfx1, 0); cf0 = BBE_MOVN_2XF32T(t0, cf0, b_ge05);
    cf1 = BBE_REPN_2XF32(cfx0, 1); t0 = BBE_REPN_2XF32(cfx1, 1); cf1 = BBE_MOVN_2XF32T(t0, cf1, b_ge05);
    cf2 = BBE_REPN_2XF32(cfx0, 2); t0 = BBE_REPN_2XF32(cfx1, 2); cf2 = BBE_MOVN_2XF32T(t0, cf2, b_ge05);
    cf3 = BBE_REPN_2XF32(cfx0, 3); t0 = BBE_REPN_2XF32(cfx1, 3); cf3 = BBE_MOVN_2XF32T(t0, cf3, b_ge05);
    cf4 = BBE_REPN_2XF32(cfx0, 4); t0 = BBE_REPN_2XF32(cfx1, 4); cf4 = BBE_MOVN_2XF32T(t0, cf4, b_ge05);
    cf5 = BBE_REPN_2XF32(cfx0, 5); t0 = BBE_REPN_2XF32(cfx1, 5); cf5 = BBE_MOVN_2XF32T(t0, cf5, b_ge05);
    cf6 = BBE_REPN_2XF32(cfx0, 6); t0 = BBE_REPN_2XF32(cfx1, 6); cf6 = BBE_MOVN_2XF32T(t0, cf6, b_ge05);
    cf7 = BBE_REPN_2XF32(cfx0, 7); t0 = BBE_REPN_2XF32(cfx1, 7); cf7 = BBE_MOVN_2XF32T(t0, cf7, b_ge05);
    #endif

    /*
    * Calcualte the polynomial: g <- atan(y)/y. Use a combination of Horner's and
    * Estrin's evaluation schemes.
    */
    BBE_MULAN_2XF32(cf1, cf0, y0);
    BBE_MULAN_2XF32(cf3, cf2, y0);
    BBE_MULAN_2XF32(cf5, cf4, y0);
    BBE_MULAN_2XF32(cf7, cf6, y0);

    y2 = BBE_MULN_2XF32(y0, y0);

    t0 = cf1;
    t1 = cf3; BBE_MULAN_2XF32(t1, t0, y2);
    t0 = cf5; BBE_MULAN_2XF32(t0, t1, y2);
    t1 = cf7; BBE_MULAN_2XF32(t1, t0, y2);
    g = t1;
    /*
    * Deduce atan(y) and restore the input range.
    */
    /* y <- y + y*g */
    BBE_MULAN_2XF32(x0, x0, g); y0 = x0;
    /* If x>=1: z <- pi/2 - y, otherwise z <- y */
    t0 = BBE_SUBN_2XF32(pi2f.f, y0);
    zout = BBE_MOVN_2XF32T(t0, y0, b_ge1);
    /* Restore the sign. */
    zout = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(zout)));
    BBE_SVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* vfastatanf() */
#endif

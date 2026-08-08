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
  NatureDSP_Baseband library. Complex Math functions
    Phase Angle Of Complex Number
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, vfastargf, ( float32_t * restrict z, const complex_float * restrict x, int N ))
#else

/* pi, pi/2 values */
#include "pif_tbl.h"
/* atan polynomial coeff table */
#include "atanf_tbl.h"

/*-------------------------------------------------------------------------
Phase Angle Of Complex Number

Description: These functions compute phase angles of input complex numbers
and return the results in radians: [-pi,pi].

Representation: IEEE-754 Std. single precision floating-point format

Accuracy: 
2 ULP - vargf, sargf
3 ULP - vfastargf

Input domain for 'fast' version vfastargf:
1.1755e-038<|real(x)|<Inf
1.1755e-038<|imag(x)|<Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
x[N]  Complex numbers
N     Length of vectors
Output:
z[N]  Phase angles

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/

void vfastargf ( float32_t * restrict z, const complex_float * restrict x, int N )
{
  /*
  * Reference C code:
  * const union ufloat32uint32* p;
  * float32_t X, Y, DEN;
  * int sx, sy, big;
  * int n;
  * 
  * for (n = 0; n<N; n++)
  * {
  *   X = crealf(x[n]);
  *   Y = cimagf(x[n]);
  *   
  *   sx = takesignf(X);
  *   sy = takesignf(Y);
  *   X = fabsf(X);
  *   Y = fabsf(Y);
  *   if (X == 0.f && Y == 0.f)
  *   {
  *     X = 1.f;
  *     Y = 0.f;
  *   }
  * 
  *   big = X>Y;
  *   if (big)
  *   {
  *     DEN = 1.0f/X;
  *     X = Y*DEN;
  *   }
  *   else
  *   {
  *     DEN = 1.0f/Y;
  *     X = X*DEN;
  *   }
  *   if (X == plusInff.f || DEN == plusInff.f)
  *   {
  *     X = qNaNf.f;
  *   }
  * 
  *   p = atanftbl_10ord;
  *   Y = p[0].f;
  *   Y = X*Y + p[1].f;
  *   Y = X*Y + p[2].f;
  *   Y = X*Y + p[3].f;
  *   Y = X*Y + p[4].f;
  *   Y = X*Y + p[5].f;
  *   Y = X*Y + p[6].f;
  *   Y = X*Y + p[7].f;
  *   Y = X*Y + p[8].f;
  *   Y = X*Y + p[9].f;
  *   Y = X*Y + X;
  * 
  *   if (!big) Y = pi2f.f -Y;
  *   if (sx)   Y = pif.f - Y;
  *   if (sy)   Y = -Y;
  *   
  *   z[n] = Y;
  */
  const xb_vecN_2xf32  * restrict pX0;
  const xb_vecN_2xf32  * restrict pX1;
  const xtfloat        * restrict pTBL;
  const xb_vecN_2xf32  * restrict pZld;
        xb_vecN_2xf32  * restrict pZst;
  int n;
  xb_vecN_2xf32 g, r, y0, zout;
  xb_vecN_2xf32 t0, den, num, eps, p;
  xb_vecNx16 x0i, y0i;
  vboolN b_s;

  xb_vecN_2xf32 x0, xmag, ymag, one, half, zero;
  vboolN_2 b_sx, b_sy, b_eqz, b_inf_den, b_xley;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7, cf8, cf9;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  /*
   * First stage:
   * make range reduction, calculate the quotinet p=num/den and save it to the output;
   * process special cases (Inf/Inf, 0/Inf or NaN/Inf).
   */

  pX0  = (const xb_vecN_2xf32 *)x;
  pZst = (      xb_vecN_2xf32 *)z;

  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX0, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pX0, 2 * BBE_SIMD_WIDTH);
    BBE_DSELN_2XF32I(y0, x0, y0, x0, BBE_DSELI_DEINTERLEAVE_2);

    zero = BBE_CONSTN_2XF32(0);
    one  = BBE_CONSTN_2XF32(1);
    half = BBE_CONSTN_2XF32(3);
    /*
     * Reduce the problem to the range of (0,pi/4].
     */
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);

    /* Select the numerator and denominator. */
    num = BBE_MINN_2XF32(xmag, ymag);
    den = BBE_MAXN_2XF32(ymag, xmag);


    /*
     * Calculate the quotient num/den.
     */
    /* Initial appromimation for 1/den. */
    r = BBE_RECIP0N_2XF32(den);
    /* Detect 0/Inf, Inf/Inf conditions. */
    b_inf_den = BBE_OEQN_2XF32(zero, r);
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
    /* Just in case, substitute zero or NaN for
     * Inf/Inf, 0/Inf or NaN/Inf cases, respectively. */
    t0 = BBE_MULN_2XF32(num, zero); /* t0 <- zero or NaN */
    p = BBE_MOVN_2XF32T(t0, p, b_inf_den);

    BBE_SVN_2XF32_IP(p, pZst, 2 * BBE_SIMD_WIDTH);
  }

  /*
   * Second stage:
   * Calculate polynomial g=atan(p)/p-1.0 where p is previously computed quotient;
   * compute atan(p)=g*p+p.
   */

  pTBL = (      xtfloat       *)atanftbl_10ord;
  pZld = (const xb_vecN_2xf32 *)z;
  pZst = (      xb_vecN_2xf32 *)z;

  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(p, pZld, 2 * BBE_SIMD_WIDTH);

    /*
     * Load coeff set
     */

    cf0 = BBE_LSN_2XF32_I(pTBL, 0*sizeof(float32_t)); cf0 = BBE_REPN_2XF32(cf0, 0);
    cf1 = BBE_LSN_2XF32_I(pTBL, 1*sizeof(float32_t)); cf1 = BBE_REPN_2XF32(cf1, 0);
    cf2 = BBE_LSN_2XF32_I(pTBL, 2*sizeof(float32_t)); cf2 = BBE_REPN_2XF32(cf2, 0);
    cf3 = BBE_LSN_2XF32_I(pTBL, 3*sizeof(float32_t)); cf3 = BBE_REPN_2XF32(cf3, 0);
    cf4 = BBE_LSN_2XF32_I(pTBL, 4*sizeof(float32_t)); cf4 = BBE_REPN_2XF32(cf4, 0);
    cf5 = BBE_LSN_2XF32_I(pTBL, 5*sizeof(float32_t)); cf5 = BBE_REPN_2XF32(cf5, 0);
    cf6 = BBE_LSN_2XF32_I(pTBL, 6*sizeof(float32_t)); cf6 = BBE_REPN_2XF32(cf6, 0);
    cf7 = BBE_LSN_2XF32_I(pTBL, 7*sizeof(float32_t)); cf7 = BBE_REPN_2XF32(cf7, 0);
    cf8 = BBE_LSN_2XF32_I(pTBL, 8*sizeof(float32_t)); cf8 = BBE_REPN_2XF32(cf8, 0);
    cf9 = BBE_LSN_2XF32_I(pTBL, 9*sizeof(float32_t)); cf9 = BBE_REPN_2XF32(cf9, 0);
    /*
     * Calculate the polynomial: g <- atan(p)/p-1.
     */
#if 0
    /* Use a combination of Estrin's and Horner's evaluation schemes */
    BBE_MULAN_2XF32(cf1, cf0, p);
    BBE_MULAN_2XF32(cf3, cf2, p);
    BBE_MULAN_2XF32(cf5, cf4, p);
    BBE_MULAN_2XF32(cf7, cf6, p);
    BBE_MULAN_2XF32(cf9, cf8, p);

    p2 = BBE_MULN_2XF32(p, p);

    t0 = cf1;
    t1 = cf3; BBE_MULAN_2XF32(t1, t0, p2);
    t0 = cf5; BBE_MULAN_2XF32(t0, t1, p2);
    t1 = cf7; BBE_MULAN_2XF32(t1, t0, p2);
    t0 = cf9; BBE_MULAN_2XF32(t0, t1, p2);
    g = t0;
#else
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
#endif
    /*
     * Deduce atan(p): p <- g*p + p.
     */
    BBE_MULAN_2XF32(p, g, p);

    BBE_SVN_2XF32_IP(p, pZst, 2 * BBE_SIMD_WIDTH);
  }

  /*
   * Third stage:
   * Restore the range and sign of the result.
   */

  pX1  = (const xb_vecN_2xf32 *)x;
  pZld = (const xb_vecN_2xf32 *)z;
  pZst = (      xb_vecN_2xf32 *)z;

  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    /*
     * Restore the quadrant.
     */
    BBE_LVN_2XF32_IP(p, pZld, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x0, pX1, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pX1, 2 * BBE_SIMD_WIDTH);
    BBE_DSELN_2XF32I(y0, x0, y0, x0, BBE_DSELI_DEINTERLEAVE_2);

    /* Get the sign bit of input values */
    x0i = BBE_MOVNX16_FROMN_2XF32(x0);
    y0i = BBE_MOVNX16_FROMN_2XF32(y0);
    x0i = BBE_SELNX16I(y0i, x0i, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
    b_s = BBE_LTNX16(x0i, BBE_ZERONX16());
    BBE_EXTRACTBN(b_sy, b_sx, b_s);

    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    /* Check y == 0 && x == 0 condition. */
    b_eqz = BBE_OEQN_2XF32(ymag, zero);
    b_eqz = BBE_OEQN_2XF32T(xmag, zero, b_eqz);
    /* Determine the dominating axis. */
    b_xley = BBE_OLEN_2XF32(xmag, ymag);

    /* |x|<=|y|: p <- pi/2 - p */
    BBE_SUBN_2XF32T(p, pi2f.f, p, b_xley);
    /* x==0 || y==0: p <- 0 */
    BBE_CONSTN_2XF32T(p, 0, b_eqz);
    /* x<0: z <- pi - p */
    BBE_SUBN_2XF32T(p, pif.f, p, b_sx);
    /* y<0: p <- -p */
    BBE_NEGN_2XF32T(p, p, b_sy);
    zout=p;
    BBE_SVN_2XF32_IP(zout, pZst, 2 * BBE_SIMD_WIDTH);
  }

} /* vfastargf() */

#endif/* !HAVE_VFPU */

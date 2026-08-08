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
DISCARD_FUN(void, vargf, ( float32_t * restrict z, const complex_float * restrict x, int N ))
#else

/* pi, pi/2 values */
#include "pif_tbl.h"
/* atan polynomial coeff table */
#include "atanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"

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

void vargf ( float32_t * restrict z, const complex_float * restrict x, int N )
{
  /*
  * Reference C code:
  *   const union ufloat32uint32* p;
  *   int sx,sy,big;
  *   sx=takesignf(x);
  *   sy=takesignf(y);
  *   x=fabs(x);
  *   y=fabs(y);
  *   if(x==0.f && y==0.f)
  *   {
  *     if ( sx )
  *       return ( sy ? -pif.f : pif.f );
  *     else
  *       return ( sy ? -0.f : 0.f );
  *   }
  *
  *   big=x>y;
  *   if(big)
  *   {
  *       x=y/x;
  *   }
  *   else
  *   {
  *       compare x==y is necessary to support (+/-Inf, +/-Inf) cases
  *       x=(x==y)? 1.0f : x/y;
  *   }
  *   p=(x<0.5f) ? atanftbl1:atanftbl2;
  *   approximate atan(x)/x-1
  *   y=    p[0].f;
  *   y=x*y+p[1].f;
  *   y=x*y+p[2].f;
  *   y=x*y+p[3].f;
  *   y=x*y+p[4].f;
  *   y=x*y+p[5].f;
  *   y=x*y+p[6].f;
  *   y=x*y+p[7].f;
  *   convert result to true atan(x)
  *   y=x*y+x;
  *
  *   if(!big) y=pi2f.f-y;
  *   if(sx)   y=pif.f -y;
  *   if(sy)   y=-y;
  *   return   y;
  */
  const xb_vecN_2xf32  * restrict pX0;
  const xb_vecN_2xf32  * restrict pX1;
  const xb_vecN_2xf32  * restrict pT0;
  const xb_vecN_2xf32  * restrict pT1;
  const xb_vecN_2xf32  * restrict pZld;
        xb_vecN_2xf32  * restrict pZst;
  int n;
  xb_vecN_2xf32 g, r, y0, zout;
  xb_vecN_2xf32 t0, den, num, eps, p;
  xb_vecNx16 x0i, y0i;
  vboolN b_s;

  xb_vecN_2xf32 x0, xmag, ymag, one, half, zero;
  vboolN_2 b_sx, b_sy, b_eqz, b_subn, b_ge05, b_inf_num, b_inf_den, b_xley;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7, cfx0, cfx1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  /*
   * First stage:
   * make range reduction, calculate the quotinet p=num/den and save it to the output;
   * process special cases (Inf/Inf, 0/Inf or NaN/Inf).
   */

  pX0  = (const xb_vecN_2xf32  *)x;
  pZst = (      xb_vecN_2xf32  *)z;

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
    /* Check for subnomals in any of x and y. */
    b_subn = BBE_ORBN_2(BBE_OLTN_2XF32(xmag, realminf.f),
    BBE_OLTN_2XF32(ymag, realminf.f));

    /* Select the numerator and denominator. */
    num = BBE_MINN_2XF32(xmag, ymag);
    den = BBE_MAXN_2XF32(ymag, xmag);
    /* Conditionally scale data up by 2^23 whenever the smaller of x and y appears
     * subnormal. This avoids two potential issues:
     *   - overflow of 1/den
     *   - loss of precision in num/den when num is subnormal. */
    BBE_MULN_2XF32T(num, num, 8388608.f, b_subn);
    BBE_MULN_2XF32T(den, den, 8388608.f, b_subn);

    /* Detect 0/Inf, Inf/Inf conditions. */
    b_inf_den = BBE_OLEN_2XF32(plusInff.f, den);
    b_inf_num = BBE_OLEN_2XF32(plusInff.f, num);
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
    /* Just in case, substitute pi/4, zero or NaN for
     * Inf/Inf, 0/Inf or NaN/Inf cases, respectively. */
    t0 = BBE_MULN_2XF32(num, zero); /* t0 <- zero or NaN */
    BBE_CONSTN_2XF32T(t0, 1, b_inf_num);/* p=1.0 -> atan=pi/4 */
    p = BBE_MOVN_2XF32T(t0, p, b_inf_den);

    BBE_SVN_2XF32_IP(p, pZst, 2 * BBE_SIMD_WIDTH);
  }

  /*
   * Second stage:
   * Calculate polynomial g=atan(p)/p-1.0 where p is previously computed quotient;
   * compute atan(p)=g*p+p.
   */

  pX1  = (const xb_vecN_2xf32 *)x;
  pZld = (const xb_vecN_2xf32 *)z;
  pZst = (      xb_vecN_2xf32 *)z;
  pT0  = (      xb_vecN_2xf32 *)atanftbl1;
  pT1  = (      xb_vecN_2xf32 *)atanftbl2;

  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(p, pZld, 2 * BBE_SIMD_WIDTH);
    /*
     * Load and select two alternative coeff sets for p<0.5 and p>=0.5
     */
    b_ge05 = BBE_OLEN_2XF32(half, p);

    BBE_LVN_2XF32_XP(cfx0, pT0, 0);
    BBE_LVN_2XF32_XP(cfx1, pT1, 0);

    cf0 = BBE_REPN_2XF32(cfx0, 0); t0 = BBE_REPN_2XF32(cfx1, 0); cf0 = BBE_MOVN_2XF32T(t0, cf0, b_ge05);
    cf1 = BBE_REPN_2XF32(cfx0, 1); t0 = BBE_REPN_2XF32(cfx1, 1); cf1 = BBE_MOVN_2XF32T(t0, cf1, b_ge05);
    cf2 = BBE_REPN_2XF32(cfx0, 2); t0 = BBE_REPN_2XF32(cfx1, 2); cf2 = BBE_MOVN_2XF32T(t0, cf2, b_ge05);
    cf3 = BBE_REPN_2XF32(cfx0, 3); t0 = BBE_REPN_2XF32(cfx1, 3); cf3 = BBE_MOVN_2XF32T(t0, cf3, b_ge05);
    cf4 = BBE_REPN_2XF32(cfx0, 4); t0 = BBE_REPN_2XF32(cfx1, 4); cf4 = BBE_MOVN_2XF32T(t0, cf4, b_ge05);
    cf5 = BBE_REPN_2XF32(cfx0, 5); t0 = BBE_REPN_2XF32(cfx1, 5); cf5 = BBE_MOVN_2XF32T(t0, cf5, b_ge05);
    cf6 = BBE_REPN_2XF32(cfx0, 6); t0 = BBE_REPN_2XF32(cfx1, 6); cf6 = BBE_MOVN_2XF32T(t0, cf6, b_ge05);
    cf7 = BBE_REPN_2XF32(cfx0, 7); t0 = BBE_REPN_2XF32(cfx1, 7); cf7 = BBE_MOVN_2XF32T(t0, cf7, b_ge05);
    /*
     * Calcualte the polynomial: g <- atan(p)/p-1. Use a combination of Estrin's
     * and Horner's evaluation schemes.
     */
#if 0
    BBE_MULAN_2XF32(cf1, cf0, p);
    BBE_MULAN_2XF32(cf3, cf2, p);
    BBE_MULAN_2XF32(cf5, cf4, p);
    BBE_MULAN_2XF32(cf7, cf6, p);

    y2 = BBE_MULN_2XF32(p, p);

    t0 = cf1;
    t1 = cf3; BBE_MULAN_2XF32(t1, t0, y2);
    t0 = cf5; BBE_MULAN_2XF32(t0, t1, y2);
    t1 = cf7; BBE_MULAN_2XF32(t1, t0, y2);
    g = t1;
#else
    BBE_MULAN_2XF32(cf1, cf0, p);
    BBE_MULAN_2XF32(cf2, cf1, p);
    BBE_MULAN_2XF32(cf3, cf2, p);
    BBE_MULAN_2XF32(cf4, cf3, p);
    BBE_MULAN_2XF32(cf5, cf4, p);
    BBE_MULAN_2XF32(cf6, cf5, p);
    BBE_MULAN_2XF32(cf7, cf6, p);
    g = cf7;
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

} /* vargf() */

#endif/* !HAVE_VFPU */

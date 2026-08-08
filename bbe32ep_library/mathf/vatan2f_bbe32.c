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
#include <errno.h>

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
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
DISCARD_FUN(void,vatan2f,( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N ))
#else
void vatan2f ( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N )
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

    int n;
  __fenv_t fenv;
  xb_vecN_2x32Uv flags, flagsx, flagsy, maskx, sgn;
  xb_vecN_2xf32 g, r, y0, y2, zout;
  xb_vecN_2xf32 t0, t1, den, num, eps, p;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32  * restrict pT0 = (     xb_vecN_2xf32  *)atanftbl1;
  const xb_vecN_2xf32  * restrict pT1 = (     xb_vecN_2xf32  *)atanftbl2;
        xb_vecN_2xf32  * restrict pZ0 = (      xb_vecN_2xf32  *)z;
        xb_vecN_2xf32  * restrict pZ1 = (      xb_vecN_2xf32  *)z;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  xb_vecN_2xf32 x0, xmag, ymag, one, half, zero;
  vboolN_2 b_sx, b_eqz_xy, b_eqz, b_subn, b_ge05, b_inf_num, b_inf_den, b_xley;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7, cfx0, cfx1;
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  b_eqz_xy = BBE_XORBN_2(b_eqz_xy, b_eqz_xy);
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

    flagsx = BBE_CLSFYN_2XF32(x0);
    flagsy = BBE_CLSFYN_2XF32(y0);
    flags = BBE_ORN_2X32U(flags, flagsx);
    flags = BBE_ORN_2X32U(flags, flagsy);
    zero = BBE_CONSTN_2XF32(0);
    one  = BBE_CONSTN_2XF32(1);
    half = BBE_CONSTN_2XF32(3);
    /*
    * Determine the target quadrant from input signs and reduce the problem to the
    * range of (0,pi/4].
    */
    /* Check y == 0 && x == 0 condition. */
    b_eqz = BBE_OEQN_2XF32(y0, zero);
    b_eqz = BBE_OEQN_2XF32T(x0, zero, b_eqz);
    b_eqz_xy = BBE_ORBN_2(b_eqz_xy, b_eqz);
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    /* Check for subnomals in any of x and y. */
    b_subn = BBE_ORBN_2(BBE_OLTN_2XF32(xmag, realminf.f),
    BBE_OLTN_2XF32(ymag, realminf.f));

    /* Substitute y=0,x=1 for all-zeros case. The actual result depends on
    * input signs. */
    BBE_CONSTN_2XF32T(xmag, 1, b_eqz);
    /* Determine the dominating axis. */
    b_xley = BBE_OLEN_2XF32(xmag, ymag);
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
    BBE_SVN_2XF32_IP(p, pZ0, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (const xb_vecN_2xf32  *)y;
  pZ0 = (xb_vecN_2xf32  *)z;
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

    /*
    * Determine the target quadrant from input signs and reduce the problem to the
    * range of (0,pi/4].
    */
    /* Check y == 0 && x == 0 condition. */
    b_eqz = BBE_OEQN_2XF32(y0, zero);
    b_eqz = BBE_OEQN_2XF32T(x0, zero, b_eqz);
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    /* Check for subnomals in any of x and y. */
    b_subn = BBE_ORBN_2(BBE_OLTN_2XF32(xmag, realminf.f),
      BBE_OLTN_2XF32(ymag, realminf.f));
    /* Substitute y=0,x=1 for all-zeros case. The actual result depends on
    * input signs. */
    BBE_CONSTN_2XF32T(xmag, 1, b_eqz);
    b_xley = BBE_OLEN_2XF32(xmag, ymag);
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
    * Load and select two alternative coeff sets for p<0.5 and p>=0.5
    */
    b_ge05 = BBE_OLEN_2XF32(half, p);

    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 0 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 0 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf0 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 1 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 1 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf1 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 2 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 2 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf2 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 3 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 3 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf3 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 4 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 4 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf4 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 5 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 5 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf5 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 6 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 6 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf6 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);
    cfx0 = BBE_LSN_2XF32_I((const xtfloat *)pT0, 7 * 4); cfx0 = BBE_REPN_2XF32(cfx0, 0);
    cfx1 = BBE_LSN_2XF32_I((const xtfloat *)pT1, 7 * 4); cfx1 = BBE_REPN_2XF32(cfx1, 0); cf7 = BBE_MOVN_2XF32T(cfx1, cfx0, b_ge05);

    /*
    * Calcualte the polynomial: g <- atan(p)/p. Use a combination of Estrin's
    * and Horner's evaluation schemes.
    */
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

    /*
    * Deduce atan(p): p <- g*p + p. Just in case, substitute pi/4, zero or NaN for
    * Inf/Inf, 0/Inf or NaN/Inf cases, respectively.
    */
    t0 = BBE_MULN_2XF32(num, zero); /* t0 <- zero or NaN */
    t1 = BBE_MOVN_2XF32T(pi4f.f, t0, b_inf_num);
    p = BBE_MOVN_2XF32T(t1, p, b_inf_den);
    BBE_MULAN_2XF32T(p, g, p, BBE_NOTBN_2(b_inf_den));
    BBE_SVN_2XF32_IP(p, pZ1, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (const xb_vecN_2xf32  *)y;
  pZ1 = (xb_vecN_2xf32  *)z;
  pZ0 = (xb_vecN_2xf32  *)z;
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
    BBE_LVN_2XF32_IP(p, pZ1, 2 * BBE_SIMD_WIDTH);

    sgn = BBE_MOVN_2X32_FROMN_2XF32(x0);
    sgn = BBE_ANDN_2X32(sgn, _80000000);
    b_sx = BBE_EQN_2X32(sgn, _80000000);
    sgn = BBE_MOVN_2X32_FROMN_2XF32(y0);
    sgn = BBE_ANDN_2X32(sgn, _80000000);
    /*
    * Determine the target quadrant from input signs and reduce the problem to the
    * range of (0,pi/4].
    */
    /* Check y == 0 && x == 0 condition. */
    b_eqz = BBE_OEQN_2XF32(y0, zero);
    b_eqz = BBE_OEQN_2XF32T(x0, zero, b_eqz);
    ymag = BBE_ABSN_2XF32(y0);
    xmag = BBE_ABSN_2XF32(x0);
    /* Check for subnomals in any of x and y. */
    b_subn = BBE_ORBN_2(BBE_OLTN_2XF32(xmag, realminf.f),
      BBE_OLTN_2XF32(ymag, realminf.f));
    /* Substitute y=0,x=1 for all-zeros case. The actual result depends on
    * input signs. */
    BBE_CONSTN_2XF32T(xmag, 1, b_eqz);
    b_xley = BBE_OLEN_2XF32(xmag, ymag);

    /*
    * Restore the quadrant.
    */
    /* |x|<=|y|: p <- pi/2 - p */
    BBE_SUBN_2XF32T(p, pi2f.f, p, b_xley);
    /* x<0: z <- pi - p */
    BBE_SUBN_2XF32T(p, pif.f, p, b_sx);
    /* y<0: p <- -p */
    p = BBE_MOVN_2XF32_FROMN_2X32(BBE_XORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(p)));
    zout=p;
    BBE_SVN_2XF32_IP(zout, pZ0, 2 * BBE_SIMD_WIDTH);
    /*
    * Perform additional analysis of input data for Error Handling
    */

  }

  {
    xb_vecN_2xf32 v_edom, v_fe_inv;
    int fe_inv, en_edom;
    fe_inv = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan;
    vboolN_2 b_snan, b_nan, b_edom, b_fe_inv;


    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20));
    v_nan = BBE_ANDN_2X32U(flags, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);

    b_edom = BBE_ORBN_2(b_nan, b_eqz_xy);
    b_fe_inv = b_snan;

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
} /* vatan2f() */
#endif

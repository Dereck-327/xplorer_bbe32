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
    Sine/Cosine
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
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Sine/Cosine 

Description: These functions compute sine or cosine of input data

Representation:
vsine,vcos,    16-bit signed fixed-point format Q15 for input/output data
ssine,scos     It is assumed that input angular values are normalized by pi.
               That is, Fixed-point functions actually compute sin(pi*x) or
               cos(pi*x)
vsinef,vcosf,  IEEE-754 Std. single precision floating-point format for
ssinef,scosf   input/output data. Input data are treated as angular values
               specified in radians. Floating-point functions limit the
               rangw of allowable input values, see note 2.

Accuracy:
2 LSB - vsine(),vcos(), ssine(), scos(),
2 ULP - vsinef(), vcosf(), ssinef(), scosf(),
3 ULP - vfastsinef(),vfastcosf()

Notes for non-fast versions:
1. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Floating-point functions require that input value belongs to the 
   closed range [-102940.0,102940.0], otherwise the respective result
   is NaN.

Input domain for 'fast' versions vfastsinef(),vfastcosf()
|x|<804.2477
The output value is not defined outside of this range or accuracy is degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vsine,vcos) or 8 (vsinef,vcosf,vfastsinef,vfastcosf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vcosf,( float32_t * restrict y, const float32_t * restrict x, int N ))
#else
void vcosf ( float32_t * restrict y, const float32_t * restrict x, int N )
{ 

  int n;
  xb_vecN_2xf32 x0, y0, xabs, x2, zsn, zcs, jf;
  xb_vecN_2x32v ji;
  xb_vecNx16 h0, h1;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2;
  xb_vecN_2xf32 sn0, sn1, sn2;
  xb_vecN_2xf32 cs0, cs1, cs2;
  xb_vecN_2x32Uv flags, flags_n, maskx;
  vboolN s;
  vboolN_2 bneg, bsc, b_nan, b_outl;
  __fenv_t fenv;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
        xb_vecN_2xf32  * restrict pZ = (      xb_vecN_2xf32  *)y;
  const xb_vecN_2xf32  * restrict pT = (const xb_vecN_2xf32  *)polysinf_tbl;
  const xtfloat  * restrict pT1 = (const xtfloat  *)pi2fc;
  xb_vecN_2xf32  * restrict pY = (xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  pT = (const xb_vecN_2xf32  *)polysinf_tbl;
  pT1 = (const xtfloat  *)pi2fc;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    pT = (const xb_vecN_2xf32  *)polysinf_tbl;
    pT1 = (const xtfloat  *)pi2fc;

    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    /*
    * Take absolute value of the argument and remember the actual sign. Also
    * check for a domain violation.
    */
    xabs = BBE_ABSN_2XF32(x0);
    x0 = BBE_MOVN_2XF32T(x0, xabs, b_nan);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);

    /* Check for outliers. */

    b_outl = BBE_OLTN_2XF32T(sinf_maxval.f, xabs, BBE_NOTBN_2(b_nan));
    /*
    * Argument reduction.
    */

    jf = BBE_MULN_2XF32(xabs, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);

    pi2fc0 = BBE_LSN_2XF32_I(pT1, 0); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(pT1, 4); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);
    pi2fc2 = BBE_LSN_2XF32_I(pT1, 8); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);

    BBE_MULSN_2XF32(xabs, jf, pi2fc0);
    BBE_MULSN_2XF32(xabs, jf, pi2fc1);
    BBE_MULSN_2XF32(xabs, jf, pi2fc2);

    /* Set bad ekements to NaN to
    *   1) achieve NaN propagation and
    *   2) force NaN result for outliers. */
    xabs = BBE_MOVN_2XF32T(qNaNf.f, xabs, b_outl);
    BBE_SVN_2XF32_IP(xabs, pY, 2 * BBE_SIMD_WIDTH);
  }
  __Pragma("no_reorder");
  pX = (const xb_vecN_2xf32  *)x;
  pY = (xb_vecN_2xf32  *)y;
  __Pragma("ymemory(pY)");
  __Pragma("ymemory(pZ)");
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(xabs, pY, 2 * BBE_SIMD_WIDTH);

    /*
    * Take absolute value of the argument and remember the actual sign.
    */
    jf = BBE_ABSN_2XF32(x0);

    jf = BBE_MULN_2XF32(jf, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);
    /*
    * Compute sine/cosine via minmax polynomial.
    */

    sn0 = BBE_LSN_2XF32_I((const xtfloat *)pT, 0); sn0 = BBE_REPN_2XF32(sn0, 0);
    sn1 = BBE_LSN_2XF32_I((const xtfloat *)pT, 4); sn1 = BBE_REPN_2XF32(sn1, 0);
    sn2 = BBE_LSN_2XF32_I((const xtfloat *)pT, 8); sn2 = BBE_REPN_2XF32(sn2, 0);

    cs0 = BBE_LSN_2XF32_I((const xtfloat *)polycosf_tbl, 0); cs0 = BBE_REPN_2XF32(cs0, 0);
    cs1 = BBE_LSN_2XF32_I((const xtfloat *)polycosf_tbl, 4); cs1 = BBE_REPN_2XF32(cs1, 0);
    cs2 = BBE_LSN_2XF32_I((const xtfloat *)polycosf_tbl, 8); cs2 = BBE_REPN_2XF32(cs2, 0);

    x2 = BBE_MULN_2XF32(xabs, xabs);
    /* Compute polynomials by Horner's method. */
    BBE_MULAN_2XF32(sn1, sn0, x2);
    BBE_MULAN_2XF32(sn2, sn1, x2);
    sn2 = BBE_MULN_2XF32(sn2, x2);
    zsn = xabs; BBE_MULAN_2XF32(zsn, sn2, xabs);

    BBE_MULAN_2XF32(cs1, cs0, x2);
    BBE_MULAN_2XF32(cs2, cs1, x2);
    zcs = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(zcs, cs2, x2);

    /*
    * Select sine/cosine and restore the sign.
    */
    /* Determine the sign and sin/cos selector from the input range. */
    h0 = BBE_MOVNX16_FROMN_2X32(ji);

    h1 = BBE_SLLINX16(h0, 14);
    {
      xb_vecNx16 tt;
      xb_vecN_2x32v tmp, _1;
      _1 = BBE_MOVN_2X32U_FROM32U(0x1); _1 = BBE_REPN_2X32(_1, 0);
      tmp = BBE_ANDN_2X32(ji, _1);
      bsc = BBE_EQN_2X32(tmp, _1);

      //  tt = BBE_LSNX16_I((const int16_t *)c_tbl, 0); tt = BBE_REPNX16(tt, 0);
      tt = BBE_MOVAV16(0x4000);
      h1 = BBE_ADDNX16(h1, tt);
    }

    h0 = BBE_SHFLNX16I(h1, BBE_SHFLI_DUPLICATE_1_EVEN);
    s = BBE_LTNX16(h0, BBE_ZERONX16());
    bneg = BBE_MOVN_2_FROMN(s);

    /* Select sine or cosine. */
    y0 = BBE_MOVN_2XF32T(zsn, zcs, bsc);
    /* Adjust the output sign. */
    BBE_NEGN_2XF32T(y0, y0, bneg);
    BBE_SVN_2XF32_IP(y0, pZ, 2 * BBE_SIMD_WIDTH);
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

    b_edom = b_nan;
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

} /* vcosf() */
#endif

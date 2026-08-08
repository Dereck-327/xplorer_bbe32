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
    Cotangent
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
/* Angular argument reduction constants. */
#include "sinf_tbl.h"
/* tan/cotan approximation polynomial coeffs. */
#include "tanf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Cotangent 

Description: These functions compute cotangent of input data.

Data format: IEEE-754 Std. single precision floating-point.
Input data are treated as angular values specified in radians.

Accuracy: 
2 ULP - vcotf(), scotf()
3 ULP - vfastcotf()

Notes for non-fast versions:
1. Cotangent functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Cotangent functions functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfastcotf():
|x|<804.2477, x!=0
The output value is not defined outside of this range or accuracy is 
degraded

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
#if HAVE_VFPU
float32_t scotf ( float32_t x )
{
  float32_t z, zc, jf, g, eps;
  xb_int32v ji;
  float32_t p, p2, p3, p4, p8;
  float32_t cf0, cf1, cf2, cf3, cf4, cf5, cf6;
  vbool1 b_sx, b_sz, b_tan, b_refine, b_outl;
  static const union ufloat32uint32 _2m128 = { 0x00200000 }; /* 2^-128 */

  xb_int32v SCF; /* Floating-point Status and Control Register values. */


  /* Get the input sign and take the absolute value. */
  /* Compute the sign adjustment term. */
  {
    xb_int32Uv hx;
    hx = XT_CLSFY_S(x);
    hx = BBE_OPERATOR_AND32U(hx, 0x1);
    b_sx = BBE_OPERATOR_NEQ32U(hx, 0);
  }
  p = XT_ABS_S(x);

  /* Check if the input value is a NaN, or does not belong to the allowed range, */
  b_outl = XT_ULT_S(tanf_maxval, p);

  /* Assert EDOM for a NaN on input. In addition, we rely on XT_UN_S to raise
  * the FE_INVALID exception for a signalling NaN. */
  if (BBE_MOVAB1(XT_UN_S(x, x)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    x = XT_ADD_S(x, x);
    return x;
  }

  /* FE_DIVBYZERO for a zero on input. */
  if (BBE_MOVAB1(XT_OEQ_S(x, XT_CONST_S(0))))
  {
    __Pragma("frequency_hint never");
    errno = ERANGE;
    x = XT_RSQRT0_S(x);
    return x;
  }
  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  BBE_MOVSCFV(XTENSA_FE_SET_EXCEPTION_ENABLE(SCF, FE_OVERFLOW));
  /*
  * Argument reduction.
  */

  /* Determine the pi/2-wide range the input value belongs to. */
  jf = XT_FIROUND_S(XT_MUL_S(p, inv2pif.f));
  /* Replace infinities and NaNs with a safe value to prevent XT_TRUNC_S from
  * raising an undesirable FE_INVALID exception. */
  XT_MOVT_S(jf, XT_CONST_S(1), b_outl);
  ji = XT_TRUNC_S(jf, 0);
  /* Force NaN on output for an outlier on input. */
  XT_MOVT_S(p, BBE_MOVF32_FROM32(BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(-1))), b_outl);

  /* Reduce p to [-pi/4,pi/4]. */
  XT_MSUB_S(p, jf, pi2fc[0].f);
  XT_MSUB_S(p, jf, pi2fc[1].f);
  XT_MSUB_S(p, jf, pi2fc[2].f);

  /*
  * Compute the polynomial approximation g(p^2) = tan(p)/p-1. We use
  * Estrin's method to evaluate the polynomial.
  */

  cf0 = polytanf_tbl[0].f;
  cf1 = polytanf_tbl[1].f;
  cf2 = polytanf_tbl[2].f;
  cf3 = polytanf_tbl[3].f;
  cf4 = polytanf_tbl[4].f;
  cf5 = polytanf_tbl[5].f;
  cf6 = polytanf_tbl[6].f;

  p2 = XT_MUL_S(p, p);
  p3 = XT_MUL_S(p2, p);
  p4 = XT_MUL_S(p2, p2);
  p8 = XT_MUL_S(p4, p4);

  XT_MADD_S(cf2, cf1, p2); cf1 = cf2;
  XT_MADD_S(cf4, cf3, p2); cf2 = cf4;
  XT_MADD_S(cf6, cf5, p2); cf3 = cf6;

  XT_MADD_S(cf1, cf0, p4); cf0 = cf1;
  XT_MADD_S(cf3, cf2, p4); cf1 = cf3;

  XT_MADD_S(cf1, cf0, p8); g = cf1;

  /* Free term of the polynomial in p^2 is zero, thus we obtain the 3rd power. */
  z = p; XT_MADD_S(z, g, p3);

  /*
  * For an even-numbered range, we have actually approximated tan(x) = 1/cotan(x).
  * Conditionally reciprocate the tangent approximation.
  */

  /* Look if the input value range is even-numbered. */
  {
    xb_vecNx16 h0, h1;
    vboolN s;
    vboolN btan;
    h0 = BBE_MOVNX16_FROM32(ji);

    h1 = BBE_SLLINX16(h0, 15);
    h0 = BBE_SELNX16I(h1, h1, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    s = BBE_LENX16(BBE_ZERONX16(),h0 );
    BBE_EXTRACTB(btan, btan, s);
    b_tan = BBE_MOV1_FROMN(btan);
  }
//  b_tan = IVP_LE32(IVP_MOVVA32(0), IVP_SLLI32(ji, 31));
  /* Newton-Raphson refinement would result in NaN whenver p is zero or
  * small enough for 1/tan(p) to overflow. Lock the refinement procedure! */
  b_refine = BBE_OPERATOR_ANDB1(b_tan, XT_ULE_S(_2m128.f, XT_ABS_S(p)));

  /* Initial appromimation for 1/tan. We rely on PDX_RECIP0 to raise the
  * FE_DIVBYZERO for a zero on input. */
 
  zc = z;
  if (vbool1_rtor_xtbool(b_tan))
  {
    zc = XT_RECIP0_S(z);
  }
  /* Conditionally perform two Newton-Raphson iterations for 1/tan. Use XT_MSUBN to
  * suppress undesired exceptions. */
  if (BBE_MOVAB1(b_refine))
  {
    eps = XT_CONST_S(1); XT_MSUBN_S(eps, zc, z);
    XT_MADD_S(zc, eps, zc);
    eps = XT_CONST_S(1); XT_MSUBN_S(eps, zc, z);
    XT_MADD_S(zc, eps, zc);
  }
  /* Adjust the sign: for odd-numbered range it must be inverted. */
  b_sz = BBE_OPERATOR_XORB1(b_sx, BBE_NOTB1(b_tan));
  /* Restore the sign. */
  if (vbool1_rtor_xtbool(b_sz))
  {
    zc = XT_NEG_S(zc);
  }

  BBE_MOVSCFV(SCF);
  /*
  * Perform additional analysis for Error Handling
  */
  {
    vbool1 b_eqz, b_ovfl;
    /* Is tan(p) zero? Make use of p==tan(p) for p ~= 0 */
    b_eqz = XT_UEQ_S(XT_CONST_S(0), p);
    /* Check if p is so small that 1/tan(p) overflows. */
    b_ovfl = BBE_OPERATOR_ANDB1(XT_ULE_S(XT_ABS_S(p), _2m128.f), b_tan);
    /* |1/tan(p)| == Inf -> assert ERANGE */
    if (vbool1_rtor_xtbool(b_ovfl))
    {
      __Pragma("frequency_hint never");
      errno = ERANGE;
      /* |1/tan(p)| == Inf AND p!=0 -> raise FE_OVERFLOW */
      if (BBE_MOVAB1(BBE_ANDNOTB1(b_ovfl, b_eqz)))
      {
        __Pragma("frequency_hint never");
        __feraiseexcept(FE_OVERFLOW);
      }
    }
   // BBE_MOVSCFV(SCF);
  }
  return (zc);

} /* scotf() */
#else
DISCARD_FUN(float32_t,scotf,( float32_t x ))
#endif

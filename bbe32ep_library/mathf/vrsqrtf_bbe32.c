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
    Reciprocal Square Root
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h" 
/*-------------------------------------------------------------------------
Reciprocal Square Root

Description: Evaluate the reciprocal square root of input value x and store
result to y: y = 1/x^0.5.

Representation:
vrsqrt,srsqrt  32-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-
               defined, provided that it is even.
               Fixed-point format Qy for output data is Qy = 31-Qx/2.
               For example, if Qx == 30 then Qy = 31-30/2 = 16
vrsrtf         IEEE-754 Std. single precision floating-point format

Accuracy:
vrsqrt,srsqrt  1.1e-4 (worst case relative error)
vrsrtf,srsqrtf 2 ULP
vfastrsrtf     3 ULP

Notes (for non-fast versions)::
1. Fixed-point functions return MIN_INT32 (0x80000000) for a negative or
   zero input value.
2. Floating-point reciprocal square root conforms to IEEE-754 Std rSqrt 
   operation in respect to signaling error conditions by means of floating-
   point exceptions.
3. Floating-point reciprocal square root limits the range of allowable
   input values, as follows:
   A) If x<0, functions raise the "invalid" floating-point exception,
      assign EDOM to errno and set output value y to NaN.
   B) If x==+/-0, functions set output value y to +/-HUGE_VALF, raise the
      "divide by zero" floating-point exception, and assign ERANGE to errno.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vrsqrt) or 8 (vrsqrtf,vfastrsqrtf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vrsqrtf,( float32_t * restrict y,
         const float32_t * restrict x,
         int N ))
#else
void vrsqrtf ( float32_t * restrict y,
         const float32_t * restrict x,
         int N )
{
  int n;
  xb_vecN_2x32v _80000000, x_fr, _7f200000;
  xb_vecN_2xf32 x0, y0, s_inf, one, half, z0_err,t0,_inf;
  vboolN_2 binf, beqz;
  xb_vecN_2x32Uv flags, flags_n;
  __fenv_t fenv;
  vboolN_2 b_edom, b_fe_inv, b_fe_divz;
  vboolN_2 b_nan;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  _7f200000 = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x7f200000));/* 1.25*2^127 */
  _80000000 = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000));
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  flags = BBE_ZERON_2X32U();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 0*2 * BBE_SIMD_WIDTH);

    x_fr = BBE_MOVN_2X32_FROMN_2XF32(x0);
    _inf = plusInff.f;
    binf = BBE_OEQN_2XF32(_inf, x0);
    beqz = BBE_OEQN_2XF32(x0, BBE_ZERON_2XF32());

    /* sign(x)*Inf */
    x_fr = BBE_ANDN_2X32(x_fr, _80000000);
    x_fr = BBE_ORN_2X32(x_fr, _7f200000);
    s_inf = BBE_MOVN_2XF32_FROMN_2X32(x_fr);
    
    half = BBE_CONSTN_2XF32(3);
    one = BBE_CONSTN_2XF32(1);
    y0 = BBE_RSQRT0N_2XF32(x0);
    /* perform additional processing to get *
    * correct result and avoid unnecessary *
    * exceptions (if x==Inf or x==+/-0)    */
    x0 = BBE_MOVN_2XF32T(BBE_ZERON_2XF32(), x0, binf);
    y0 = BBE_MOVN_2XF32T(s_inf, y0, beqz);

    /* compute approximation error */
    z0_err = one;
    t0 = BBE_MULN_2XF32(x0, y0);
    BBE_MULSN_2XF32(z0_err, t0, y0);

    /* 2-nd rsqrt approximation */
    t0 = BBE_MULN_2XF32(half, y0);
    BBE_MULANN_2XF32(y0, t0, z0_err);

    /* compute approximation error */
    z0_err = one;
    t0 = BBE_MULN_2XF32(x0, y0);
    BBE_MULSN_2XF32(z0_err, t0, y0);

    /* 3-rd rsqrt approximaiton */
    t0 = BBE_MULN_2XF32(half, y0);
    BBE_MULANN_2XF32(y0, t0, z0_err);
    /*
    * Perform additional analysis of input data for Error Handling.
    */
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b_nan = BBE_UNN_2XF32(x0, x0);
    BBE_ABSN_2XF32T(x0, x0, b_nan);
    flags_n = BBE_CLSFYN_2XF32(x0);
    flags = BBE_ORN_2X32U(flags, flags_n);
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
  }
  {
    xb_vecN_2xf32 v_edom, v_fe_inv, v_fe_divz;
    int fe_inv, fe_divz, en_edom;
    xb_vecN_2x32Uv maskx;
    fe_inv = 0; fe_divz = 0; en_edom = 0;
    xb_vecN_2x32Uv v_snan, v_nan, v_eqz, v_ltz;
    vboolN_2 b_snan, b_nan, b_eqz, b_ltz;


    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40)); 
    v_snan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x2)); 
    v_eqz = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x20)); 
    v_nan = BBE_ANDN_2X32U(flags, maskx);
    maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x1)); 
    v_ltz = BBE_ANDN_2X32U(flags, maskx);

    maskx = BBE_ZERON_2X32U();
    b_snan = BBE_NEQN_2X32U(v_snan, maskx);
    b_eqz = BBE_NEQN_2X32U(v_eqz, maskx);
    b_nan = BBE_NEQN_2X32U(v_nan, maskx);
    b_ltz = BBE_NEQN_2X32U(v_ltz, maskx);

    b_edom = BBE_ORBN_2(b_ltz, b_snan);
    b_edom = BBE_ORBN_2(b_edom, b_nan);

    b_fe_inv = BBE_ORBN_2(b_ltz, b_snan);
    b_fe_divz = b_eqz;

    v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
    v_fe_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_fe_inv);
    v_fe_divz = BBE_MOVN_2XF32T(FE_DIVBYZERO, BBE_ZERON_2XF32(), b_fe_divz);

    /* Retrieve EDOM state: some x<0 OR x==sNaN OR x==qNaN */
    en_edom = BBE_RMAXNUMN_2XF32(v_edom);
    /* Merge FE_INVALID state: some x<0 OR x==sNaN (the latter is detected by hardware) */
    fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);
    /* Retrieve FE_DIVBYZERO state: some x==0 */
    fe_divz = BBE_RMAXNUMN_2XF32(v_fe_divz);

    /* For logf(), the ERANGE state is equivalent to FE_DIVBYZERO floating-point exception. */
    if (0 != fe_divz) { __Pragma("frequency_hint never"); errno = ERANGE; };
    /* EDOM takes precedence over ERANGE! */
    if (0 != en_edom) { __Pragma("frequency_hint never"); errno = EDOM; };

    /* Restore exception enable flags and status flags, suppress undesired status flags. */
    __fesetenv(&fenv);
    /* Raise FE_INVALID (x<0 or x==sNaN) and/or FE_DIVBYZERO (x==0). */
    __feraiseexcept(fe_inv | fe_divz);
  }
} /* vrsqrtf() */
#endif

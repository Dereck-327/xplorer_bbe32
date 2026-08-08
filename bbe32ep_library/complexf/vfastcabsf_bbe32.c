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
    Magnitude Of Complex Number
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
DISCARD_FUN(void, vfastcabsf, ( float32_t * restrict z, const complex_float * restrict x, int N ))
#else

/*-------------------------------------------------------------------------
Magnitude Of Complex Number

Description: These functions multiply complex input number by its conjugate
and take the square root of the result.

Representation:
vcabs          16-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-defined.
               Fixed-point format for output data is Qz = Qx+sh/2, where sh
               is the shift control argument. If a resulting value is too
               large to be represented in Qz format, then it is saturated
               by 32767.
vcabsf,scabsf  IEEE-754 Std. single precision floating-point format


Input domain for 'fast' version vfastcabsf():
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
1 LSB - vcabs
2 ULP - vcabsf,scabsf
3 ULP - vfastcabsf

Parameters:
Input:
x[N]  Complex numbers
sh    (vcabs only) bi-directional shift control, an even integer from
      the range [-2,30]. Positive value corresponds to a left shift
N     Length of vectors
Output:
z[N]  Magnitudes

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vcabs) or 8 (vcabsf,vfastcabsf)
-------------------------------------------------------------------------*/

void vfastcabsf ( float32_t * restrict z, const complex_float   * restrict x, int N )
{
  const xb_vecN_2xf32 *restrict X;
        xb_vecN_2xf32 *restrict Z;
  xb_vecN_2xf32 Xre, Xim, Zout;
  xb_vecN_2xf32 half, t0, t1, x0, x0_red, x0_adj, r0, r0_err;
  int n;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N%(BBE_SIMD_WIDTH/2) == 0);

  half = BBE_CONSTN_2XF32(3);

  X = (const xb_vecN_2xf32 *)x;
  Z = (      xb_vecN_2xf32 *)z;

  for ( n = 0; n < (N>>(LOG2_BBE_SIMD_WIDTH-1)); n++ )
  {
    BBE_LVN_2XF32_IP( t0, X, 2*BBE_SIMD_WIDTH );
    BBE_LVN_2XF32_IP( t1, X, 2*BBE_SIMD_WIDTH );
    BBE_DSELN_2XF32I(Xim, Xre, t1, t0, BBE_DSELI_DEINTERLEAVE_2);

    /* Find the magnitude of the normalized complex data */
    /* Zout = sqrt(Xre*Xre + Xim*Xim) */
    x0 = BBE_MULN_2XF32(Xre, Xre);
    BBE_MULAN_2XF32(x0, Xim, Xim);
    x0_red = BBE_NEXP01N_2XF32(x0);/* negated x with reduced exponent range */
    /* Initial rsqrt approximation with exponent range reduction */
    r0 = BBE_SQRT0N_2XF32(x0);
    r0 = BBE_NEGN_2XF32(r0);
    /* compute approximation error */
    t0 = BBE_MULN_2XF32(r0, r0);
    t1 = x0_red; BBE_ADDEXPN_2XF32(t1, half);/* -0.5*x */
    r0_err = half; BBE_MULANN_2XF32(r0_err, t0, t1);/* approximation error is (0.5-0.5*x*r*r) */
    BBE_MULANN_2XF32(r0, r0, r0_err);/* Second recip sqrt approximation */
    /* Compute reduced range sqrt approximation */
    Zout = BBE_MULN_2XF32(x0_red, r0);/* z = x*rsqrt(x) */
    /* Make final adjustment and restore range */
    x0_adj = BBE_MKSADJN_2XF32(x0);
    t0 = x0_red;
    BBE_MULANN_2XF32(t0, Zout, Zout);
    BBE_ADDEXPMN_2XF32(Zout, x0_adj);
    t1 = BBE_MULN_2XF32(half, r0);
    BBE_ADDEXPN_2XF32(t1, x0_adj);
    BBE_DIVNN_2XF32(Zout, t0, t1);
    /*--------------------------------*/

    BBE_SVN_2XF32_IP( Zout, Z, 2*BBE_SIMD_WIDTH );
  }

} /* vfastcabsf() */

#endif/* !HAVE_VFPU */

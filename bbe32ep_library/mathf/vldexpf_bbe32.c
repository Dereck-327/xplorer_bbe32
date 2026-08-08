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
    Modify the Exponent of a Floating-Point Number
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
#include <errno.h>
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Modify the Exponent of a Floating-Point Number

Description: These functions multiply input values by 2^n, where n is an
exponent adjustment term. If the result overflows, functions return 
HUGE_VALF with the proper sign. If, on the contrary, the result underflows,
functions return zero with proper sign.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: exact

Note:
Exponent modification functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-point
exceptions.

Parameters:
Input:
x[N]    Input data
n[N]    Exponent adjustment terms
N       Length of input/output data vectors
Output:
y[N]    Results

Restrictions:
y,x,n   Aligned on 32-byte boundary
y,x,n   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vldexpf,( float32_t * restrict y, 
         const float32_t * restrict x, 
         const int32_t   * restrict n, 
         int N ))
#else
void vldexpf ( float32_t * restrict y, 
         const float32_t * restrict x, 
         const int32_t   * restrict n, 
         int N )
{
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32 *)x;
  const xb_vecN_2x32v  * restrict pN = (const xb_vecN_2x32v *)n;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32 *)y;
  xb_vecN_2xf32 x0, z0;
  xb_vecN_2xf32 y0, y1, y2;
  xb_vecN_2x32v n0;
  xb_vecNx16 h0, e0, e1, e2;
  xb_vecNx16  I0, I1;
  int i = 0;  
  vboolN_2 b_edom;
  __fenv_t fenv;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(n, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  I0 = 127;
  I1 = -126;
  
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  b_edom = BBE_XORBN_2(b_edom,b_edom);
  for (i = 0; i<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); i++)
  {
    /* load inputs */
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2X32_IP(n0, pN, 2 * BBE_SIMD_WIDTH);
    /* Assert EDOM whenever a NaN is encounetered. */
    b_edom = BBE_ORBN_2(b_edom, BBE_UNN_2XF32(x0, x0));

    /* Sensible range of exponent adjustment is [-278,277], as shown below:
    *   - max finite -> zero: pow2(single((2-2^-23)*2^127),-278) == 0
    *   - min subnormal -> Infinity: pow2(single(2^-149),277) == Infinity */
    n0 = BBE_MAXN_2XF32(-278, BBE_MINN_2XF32(277, n0));
    h0 = BBE_MOVNX16_FROMN_2X32(n0);
    /* Multiplication by a power of 2 is the most convenient way to modify
    * the exponent. Maximum representable power of two for the single precision
    * is 127, thus we have to split n into 3 components to cover the full range
    * of exponent adjustment: [-278,277]. The value is partitioned in such a way
    * that a subnormal result may appear only once in three products, excepting
    * trivial multiplucations by 1. This guarantees that the subnormal result is
    * rounded only once. */
    e2 = BBE_MAXNX16(I1, BBE_MINNX16(I0, h0));
    e1 = BBE_MAXNX16(I1, BBE_MINNX16(I0, BBE_SUBNX16(h0, e2)));
    e0 = BBE_SUBNX16(BBE_SUBNX16(h0, e2), e1);

    e0 = BBE_ADDNX16(I0, e0);
    e1 = BBE_ADDNX16(I0, e1);
    e2 = BBE_ADDNX16(I0, e2);

    e0 = BBE_SLLINX16(e0, 7);
    e1 = BBE_SLLINX16(e1, 7);
    e2 = BBE_SLLINX16(e2, 7);

    e0 = BBE_SELNX16I(e0, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
    e1 = BBE_SELNX16I(e1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
    e2 = BBE_SELNX16I(e2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);

    y0 = BBE_MOVN_2XF32_FROMNX16(e0);
    y1 = BBE_MOVN_2XF32_FROMNX16(e1);
    y2 = BBE_MOVN_2XF32_FROMNX16(e2); 

    /* Multiply the input x by those factors. We rely on BBE_MULN_2XF32 is raising
    * floating-point exceptions: FE_OVERFLOW (z overflows) and FE_INVALID (x is sNaN). */
    z0 = BBE_MULN_2XF32(x0, y0);
    z0 = BBE_MULN_2XF32(z0, y1);
    z0 = BBE_MULN_2XF32(z0, y2);

    BBE_SVN_2XF32_IP(z0, pY, 2*BBE_SIMD_WIDTH);
  }
    {
      xb_vecN_2xf32 v_edom;
      int fe_flags;
      fe_flags = 0; 


      v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
      /* Check if FE_OVERFLOW and/or FE_INVALID exception flags have been set by the core loop. */
      fe_flags = __fetestexcept(FE_OVERFLOW | FE_INVALID);

      /* For ldexp(), the ERANGE state is equivalent to FE_OVERFLOW floating-point exception. */
      if (0 != (fe_flags & FE_OVERFLOW)) { __Pragma("frequency_hint never"); errno = ERANGE; };
      /* EDOM takes precedence over ERANGE! */
      if (0 != (int)BBE_RMAXNUMN_2XF32(v_edom)) { __Pragma("frequency_hint never"); errno = EDOM; };

      /* Restore exception enable flags and status flags, suppress undesired status flags. */
      __fesetenv(&fenv);
      /* Raise FE_OVERFLOW and/or FE_INVALID. */
      __feraiseexcept(fe_flags);
    }
} /* vldexpf() */
#endif

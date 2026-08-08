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
    Reciprocal
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Infinities for single precision routines */
#include "inff_tbl.h"
#include <float.h>
/*-------------------------------------------------------------------------
Reciprocal

Description: Evaluate the reciprocal of input value x and store result
to y: y = 1/x.

Representation:
vrecip16,srecip16  16-bit signed fixed-point format
vrecipf            IEEE-754 Std. single precision floating-point format

Fixed-point routines compute reciprocals for Q15 input data, and return the
fractional and exponential parts of the result. Since the reciprocal of 
a 16-bit Q15 is at least 1.0 in magnitude, functions return fractional part
frac in Q(15-exp) format, where exp is the exponential part of the respective
result. Full result can be restored in 32-bit Q15 format by sign extending
the fractional part to 32 bits and shifting it to the left by exp bit positions.
Scalar fixed-point function returns packed 32-bit result, where exponential part
resides in 16 MSBs and fractional part is located in 16 LSBs.

Special cases:
      x    |  Result |  Extra Conditions    
  ---------|---------|---------------------
   +/-Inf  | +/-0    | vrecipf,srecipf
   +/-0    | +/-Inf  |  
  ---------|---------|---------------------
     0     |   not   | vrecip16,srecip16
           | defined |

Input domain for vfastrecipf:
|x|>2.94e-39, |x|<Inf
The output value is not defined outside of this range.

Accuracy:
vrecip16,srecip16  1 LSB of the fractional part
vrecipf,srecipf    1 ULP 
vfastrecipf        2 ULP 

Parameters:
Input:
x[N]           Input data vector
N              Length of vectors
Output:
vrecip16,srecip16
frac[N]        Fractional part of reciprocals, Q(15-exp)
exp[N]         Exponent of reciprocals (1...16)
vrecip:
y[N]           Reciprocals

Restrictions:
y,x,fract,exp  Aligned on 32-byte boundary
y,x,fract,exp  Must not overlap
N              Multiple of 16 (vrecip16), 8 (vrecipf,vfastrecipf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vrecipf,( float32_t * restrict y,
          const float32_t * restrict x,
          int N ))
#else
void vrecipf ( float32_t * restrict y,
          const float32_t * restrict x,
          int N )
{
  int n;
  xb_vecN_2x32v sgn, _80000000, zden, _7f800000;
  xb_vecN_2xf32 x0, y0, xabs, x_fr, _00200000;
  vboolN_2 binf, bsub;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Table of different constants used in computations */
  static const ALIGN(32) int32_t c_tbl[] =
  {
    (int32_t)0x80000000,
    (int32_t)0x7f800000,
    (int32_t)0x00200000,
  };
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    x_fr = BBE_MOVN_2X32_FROMN_2XF32(x0);
    _80000000 = BBE_LSN_2X32_I((const xb_int32v *)c_tbl, 0); _80000000 = BBE_REPN_2X32(_80000000, 0);
    sgn = BBE_ANDN_2X32(x_fr, _80000000);
    _7f800000 = BBE_LSN_2X32_I((const xb_int32v *)c_tbl, 4); _7f800000 = BBE_REPN_2X32(_7f800000, 0);

    zden = BBE_ORN_2X32(sgn, _7f800000);
    xabs = BBE_ABSN_2XF32(x0);
    binf = BBE_OEQN_2XF32(xabs, BBE_MOVN_2XF32_FROMN_2X32(_7f800000));
    _00200000 = BBE_LSN_2XF32_I((const xtfloat *)c_tbl, 8);
    _00200000 = BBE_REPN_2XF32(_00200000, 0); bsub = BBE_OLTN_2XF32(xabs, _00200000);
    y0 = BBE_RECIPN_2XF32(x0);

    y0 = BBE_MOVN_2XF32T(BBE_MOVN_2XF32_FROMN_2X32(sgn), y0, binf);
    y0 = BBE_MOVN_2XF32T(BBE_MOVN_2XF32_FROMN_2X32(zden), y0, bsub);
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
  }
} /* vrecipf() */
#endif

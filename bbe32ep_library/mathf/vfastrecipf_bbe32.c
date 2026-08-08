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
DISCARD_FUN(void,vfastrecipf, ( float32_t * restrict y, 
          const float32_t * restrict x, 
          int N ))
#else
void vfastrecipf ( float32_t * restrict y, 
          const float32_t * restrict x, 
          int N )
{
  int n;
  xb_vecN_2xf32 x0, y0;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    y0 = BBE_RECIPN_2XF32(x0); 
    BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
  }
} /* vfastrecipf() */
#endif

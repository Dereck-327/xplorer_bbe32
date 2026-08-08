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
  NatureDSP_Baseband library. Vector Operations
    Square Root of Sum of Squares
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Vector Operations. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Square Root of Sum of Squares 

Description: These routines compute the magnitude of a vector

Representation:
vmag    Signed fixed-point format
        Input vector elements are 16-bit signed data of arbitrary format
        Qx. Sum of squared values is computed in a 40-bit accumulator,
        which is then shifted to the right by rsh bit positions. Square
        root approximation of the shifted sum is formatted as a 16-bit
        non-negative value with Qx-rsh/2 fractional bits.
vmagf   IEEE-754 Std. single precision floating-point format for
        input vector and magnitude estimation result

Parameters:
Input:
x[N]    Input data
rsh     Right shift amount (vmag)
N       Length of vector 
Returned Value:
Square root of the sum of squares over the input vector

Restrictions:
x       Aligned on 32-byte boundary
N       Multiple of 16 (vmag) or 8 (vmagf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(float32_t,vmagf,( const float32_t * restrict x,
                  int N ))
#else
float32_t vmagf ( const float32_t * restrict x,
                  int N )
{
  int n;
  xtfloat f;
  xb_vecN_2xf32 x0, Acc, tmp;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return 0.f;
  Acc = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    tmp = BBE_MULN_2XF32(x0, x0);
    Acc = BBE_ADDN_2XF32(Acc, tmp);
  }
  f = BBE_RADDN_2XF32(Acc);
  f = XT_SQRT_S(f);
  return f;
} /* vmagf() */
#endif

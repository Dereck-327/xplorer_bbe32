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
    Dot Product of Real Vectors
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
Dot Product of Real Vectors

Description: These routines take two real vectors and calculate their dot 
product.

Representation:
rvdot   Signed fixed-point format
        Input vectors x and y are 16-bit signed data of arbitrary formats
        Qx and Qy. Dot product is computed in a 40-bit accumulator, which
        is then rounded, shifted to the right by rsh bit positions and
        saturated to form a 32-bit result with Qx+Qy-rsh fractional bits.
rvdotf  IEEE-754 Std. single precision floating-point format for
        input vectors and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (rvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (rvdot) or 8 (rvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(float32_t,rvdotf,( const float32_t * restrict x,
                   const float32_t * restrict y,
                   int N ))
#else
float32_t rvdotf ( const float32_t * restrict x,
                   const float32_t * restrict y,
                   int N )
{
  int n;
  xb_vecN_2xf32 x0, y0, x1, y1, Acc0, tmp, Acc1, tmp1;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0); 
  if (N <= 0) return 0.f; 
  Acc0 = BBE_ZERON_2XF32();
  Acc1 = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH)); n++) 
  { 
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    tmp = BBE_MULN_2XF32(x0, y0);
    Acc0 = BBE_ADDN_2XF32(Acc0,tmp);

    BBE_LVN_2XF32_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, pY, 2 * BBE_SIMD_WIDTH);
    tmp1 = BBE_MULN_2XF32(x1, y1);
    Acc1 = BBE_ADDN_2XF32(Acc1, tmp1);
  }
  if (N & (BBE_SIMD_WIDTH / 2))
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    tmp = BBE_MULN_2XF32(x0, y0);
    Acc0 = BBE_ADDN_2XF32(Acc0, tmp);
  }
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  return BBE_RADDN_2XF32(Acc0); 
} /* rvdotf() */
#endif

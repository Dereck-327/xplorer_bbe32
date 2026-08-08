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
    Real Vectors Sum
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
Real Vectors Sum

Description: These routines perform pairwise summation of real vectors.

Representation:
rvadd   Signed fixed-point format. 16-bit inputs, 16-bit saturated results
rvaddf  IEEE-754 Std. single precision floating-point format for input/output
        data

Parameters:
Input:
x[N]   Input vector
y[N]   Input vector
N      Length of vectors 
Output:
z[N]   Sum of input vectirs

Restrictions:
z,x,y  Must not overlap
z,x,y  Aligned on 32-byte boundary
N      Multiple of 16 (rvadd) or 8 (rvaddf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,rvaddf,( float32_t * restrict z,
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N ))
#else
void rvaddf ( float32_t * restrict z,
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N )
{
  int n;
  xb_vecN_2xf32 x0, y0, z0;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pZ = (      xb_vecN_2xf32  *)z;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return; 
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    z0 = BBE_ADDN_2XF32(x0, y0);
    BBE_SVN_2XF32_IP(z0, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* rvaddf() */
#endif

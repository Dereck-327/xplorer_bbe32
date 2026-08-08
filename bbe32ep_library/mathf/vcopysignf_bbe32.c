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
    Copy Sign
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
Copy Sign

Description: These functions copy the sign bit of the second argument to the
first argument and return the result.

Data format: IEEE-754 Std. single precision floating-point.

Note:
In terms of IEEE-754 Std, this is a quiet-computational operation which
treats floating-point numbers and NaNs alike, and does not raise any 
floating-point exceptions.

Parameters:
Input:
x[N]    Input data
y[N]    Sign data
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vcopysignf,( float32_t * restrict z, 
            const float32_t * restrict x, 
            const float32_t * restrict y, 
            int N ))
#else
void vcopysignf ( float32_t * restrict z, 
            const float32_t * restrict x, 
            const float32_t * restrict y, 
            int N )
{
  int n;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32 *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32 *)y;
        xb_vecN_2xf32  * restrict pZ = (      xb_vecN_2xf32 *)z;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  xb_vecN_2xf32 x0, y0, z0;
  xb_vecN_2x32Uv ux, uy, signmask, absmask;
  if (N <= 0) return;
  signmask = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x80000000U));
  absmask = BBE_NOTN_2X32U(signmask);
  for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

    /* Take the sign from one and put it to sign bit of second operand */
    ux = BBE_MOVN_2X32U_FROMNX16(BBE_MOVNX16_FROMN_2XF32(x0));
    uy = BBE_MOVN_2X32U_FROMNX16(BBE_MOVNX16_FROMN_2XF32(y0));
    uy = BBE_ANDN_2X32U(uy, signmask);
    ux = BBE_ANDN_2X32U(ux, absmask);
    ux = BBE_ORN_2X32U(ux, uy);
    z0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_2X32U(ux));
    BBE_SVN_2XF32_IP(z0, pZ, 2 * BBE_SIMD_WIDTH); 
  }
} /* vcopysignf() */
#endif

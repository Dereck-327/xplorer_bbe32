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
    Integer To Floating Value Conversion
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
Integer To Floating Value Conversion

Description: These functions convert integer input values to floating 
values and scale them by 2^t.

Data format: Signed 32-bit integer on input, 
             IEEE-754 Std. single precision floating-point on output.

Parameters:
Input:
x[N]  Input integers
N     Length of input/output data vectors
t     Scale factor
Output:
y[N]  Conversion results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
t     Must belong to [-126,126]
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vint2float,( float32_t * restrict y, 
            const int32_t   * restrict x,
            int t, int N ))
#else
void vint2float ( float32_t * restrict y, 
            const int32_t   * restrict x,
            int t, int N )
{
  int n;
  const xb_vecN_2x32v  * restrict pX = (const xb_vecN_2x32v *)x;
        xb_vecN_2xf32  * restrict pY = (      xb_vecN_2xf32 *)y;
  xb_vecN_2x32v x0;
  xb_vecN_2xf32 zout, fx0, sc;
  uint32_t s;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;

  /* Make the scaling coefficient */
  s = ((uint32_t)(t + 127)) << 23;
  sc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVVA16C(s));/* sc=2^t */
  for (n = 0; n < (N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2X32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    /* Convert input values to the floating-point format */
    fx0 = BBE_FLOATN_2X32(x0, 0);
    /* Scale values */ 
    zout = BBE_MULN_2XF32(fx0,sc);
    /* Save result*/
    BBE_SVN_2XF32_IP(zout, pY, 2 * BBE_SIMD_WIDTH);
  }
} /* vint2float() */
#endif

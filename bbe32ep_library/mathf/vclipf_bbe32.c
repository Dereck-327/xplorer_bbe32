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
    Clipping
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
Clipping

Description: These functions limit the absolute value of inputs. If 
magnitude of input value is less than the second argument then it is left
unchanged. Otherwise it is replaced by absolute value of the second argument
or its negation, depending on the sign of the input value.

Data format: IEEE-754 Std. single precision floating-point.

Parameters:
Input:
x[N]  Input data
y     Limiting value
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vclipf,( float32_t * restrict z, 
        const float32_t * restrict x, 
              float32_t            y,
        int N ))
#else
void vclipf ( float32_t * restrict z, 
        const float32_t * restrict x, 
              float32_t            y,
        int N )
{
  int n;
  const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
        xb_vecN_2xf32 * restrict pZ = (      xb_vecN_2xf32 *)z;


  if (N <= 0) { return; }
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N>0 && N % (BBE_SIMD_WIDTH/2) == 0);
  xb_vecN_2xf32 x0, z0, y0;
  xb_vecN_2xf32 yabs, xabs;
  xb_vecN_2x32v sgn, _80000000;
  vboolN_2 b_nan;
  y0 = y;
  yabs = BBE_ABSN_2XF32(y0);
  b_nan = BBE_UNN_2XF32(y0, y0);
  _80000000 = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000));
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    
    xabs = BBE_ABSN_2XF32(x0); 
    z0 = BBE_MINNUMN_2XF32(xabs, yabs);
    sgn = BBE_ANDN_2X32(BBE_MOVN_2X32_FROMN_2XF32(x0), _80000000);
    z0 = BBE_MOVN_2XF32_FROMN_2X32(BBE_ORN_2X32(sgn, BBE_MOVN_2X32_FROMN_2XF32(z0)));
    z0 = BBE_MOVN_2XF32T(x0,z0,b_nan);
    BBE_SVN_2XF32_IP(z0, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* vclipf() */
#endif

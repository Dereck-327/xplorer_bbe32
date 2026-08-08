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
    Common Block Exponent
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
Common Block Exponent

Description: These functions compute base-2 exponent adjustment term needed
to normalize data in the input vector. Exact meaning of normalization depends
on the data format (see below).

Representation: 
vbexp,vbexp_fast        16-bit signed fixed-point format
                        For each input value functions count the number of
                        redundant sign bits (as if the value was loaded in
                        a 32-bit register) and return the minimum result over
                        the input data vector.
sbexp                   32-bit signed fixed-point format
                        Count the number of redundant sign bits and return 
                        the result.
vbexpf,sbexpf           IEEE-754 Std. single precision floating-point format
                        For each finite input value x, functions estimate the
                        integer E(x), such that 0.5 <= |x|*2^E(x) < 1 and
                        -128 <= E(x) <= 148. The minimum value of E(x) over
                        input data vector is the result.

Special cases:
   x    |  Result |    Extra Conditions    
--------+---------|---------------------------
0       |    0    |
+/-Inf  | -129    | floating-point functions
NaN     |    0    |
--------|---------|---------------------------
0       |   31    |
-32768  |   16    | fixed-point functions
32767   |   16    |

Parameters:
Input:
x[N]    Input data
N       Length of data vector
Returned Value: exponent adjustment term, or zero if N<=0
Restrictions:
vbexp(), vbexpf():
  No restrictions
vbexp_fast():
  x     Aligned on 32-byte boundary
  N     Multiple of 16
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(int,vbexpf,( const float32_t * restrict x, int N ))
#else
int vbexpf ( const float32_t * restrict x, int N )
{
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  unsigned e;
  xb_vecN_2xf32 x0, max0, max1, x1;
  int32_t exp0, u0;
  xtfloat f, fmax0, fmax1;
  xtbool zero, bnan;
  valign x_align;
  int n, M;
  NASSERT(x);
  if (N <= 0) return 0;
  fmax0 = fmax1 = XT_MOV_S(0.f);
  x_align = BBE_LAN_2XF32_PP(pX);
  max0 = BBE_ZERON_2XF32();
  max1 = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH)); n++)
  {
    BBE_LAN_2XF32_IP(x0, x_align, pX);
    BBE_LAN_2XF32_IP(x1, x_align, pX);
    x0 = BBE_ABSN_2XF32(x0);
    x1 = BBE_ABSN_2XF32(x1);
    max0 = BBE_MAXNUMN_2XF32(max0, x0);
    max1 = BBE_MAXNUMN_2XF32(max1, x1);
  }
  
  M = (N & ((BBE_SIMD_WIDTH)-1)) *sizeof(*x);
  BBE_LAVN_2XF32_XP(x0, x_align, pX, M);
  BBE_LAVN_2XF32_XP(x1, x_align, pX, M - (BBE_SIMD_WIDTH/2) * sizeof(*x));
  x0 = BBE_ABSN_2XF32(x0);
  x1 = BBE_ABSN_2XF32(x1);
  max0 = BBE_MAXNUMN_2XF32(max0, x0);
  max1 = BBE_MAXNUMN_2XF32(max1, x1);
  max0 = BBE_MAXNUMN_2XF32(max0, max1);
  f = BBE_RMAXNUMN_2XF32(max0);

  zero = XT_UN_S(f, f); 
  bnan = XT_OEQ_S(f, 0.f);
  zero = XT_ORB(zero, bnan);
  exp0 = 0;
  u0 = XT_RFR(f);
  u0 = ((uint32_t)u0) >> 23;
  if (u0 == 0)
  {
    f = XT_MUL_S(f, 16777216.f);
    exp0 -= 24;
  }
  u0 = XT_RFR(f);
  u0 = ((uint32_t)u0) >> 23;
  exp0 += ((int)(u0)-126);
  e = exp0;
  XT_MOVT(e, 0, zero);
  return -(int32_t)e;
} /* vbexpf() */
#endif

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
    Mean of Vector Elements
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
Mean of Vector Elements

Description: Compute the mean value over all elements of given vector (real 
or complex).

Representation:
rvmean,cvmean    16-bit signed fixed-point format
rvmeanf,cvmeanf  IEEE-754 Std. single precision floating-point format

Parameters:
Input:
x[N]    Input vector
N       Length of input vector, in real or complex samples
Output:
m[1]    Mean value, real or complex

Restrictions:
x       Aligned on 32-byte boundary
x,m     Must not overlap
N       Must be a multiple of either:
          4 (cvmeanf), or
          8 (cvmean, rvmeanf), or
          16 (rvmean)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,cvmeanf,( complex_float   * m, const complex_float   * restrict x, int N ))
#else
void cvmeanf ( complex_float   * m, const complex_float   * restrict x, int N )
{
  int n;
  xb_vecN_2xf32 x0;
  xb_vecN_2xf32 Acc0, Acc1, Acc2, Acc3;

  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 4) == 0);
 
  if (N <= 0) { MOV_COMPLEXFLOAT_FROM_N_4XCF32(m[0], BBE_ZERON_4XCF32()); return;}
  Acc0 = BBE_ZERON_2XF32();
  Acc1 = BBE_ZERON_2XF32();
  Acc2 = BBE_ZERON_2XF32();
  Acc3 = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    Acc0 = BBE_ADDN_2XF32(Acc0, x0);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    Acc1 = BBE_ADDN_2XF32(Acc1, x0);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    Acc2 = BBE_ADDN_2XF32(Acc2, x0);
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    Acc3 = BBE_ADDN_2XF32(Acc3, x0);
  }
  for (n = 0; n<(N & ((BBE_SIMD_WIDTH / 4) * 3)); n += (BBE_SIMD_WIDTH / 4))
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    Acc0 = BBE_ADDN_2XF32(Acc0, x0);
  }
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  Acc2 = BBE_ADDN_2XF32(Acc3, Acc2);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
  Acc1 = BBE_SELN_2XF32I(Acc0, Acc0, BBE_SELI_ROTATE_RIGHT_8);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  Acc1 = BBE_SELN_2XF32I(Acc0, Acc0, BBE_SELI_ROTATE_RIGHT_4);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);

  Acc0 = BBE_DIVN_2XF32(Acc0,N);
  {
    xb_vecN_4xcf32 t;
    t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc0));
    MOV_COMPLEXFLOAT_FROM_N_4XCF32(m[0], t);
  }
} /* cvmeanf() */
#endif

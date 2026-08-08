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
    Find Values Above Threshold
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
Find Values Above Threshold 

Description: function collects indices of all input vector elements that hold
a value greater than the designated threshold.

Representation:
vthreshold   16-bit signed fixed-point format
vthresholdf  IEEE-754 Std. single precision floating-point format

Parameters:
Input:
thr       Threshold value
x[N]      Input data
N         Size of input array
Output:
idx[N+1]  Indices of elements with a value greater the threshold
NOTE: 
extra cell should be reserved in the end of array
Returns:  Number of found elements

Restrictions:
x,idx     Aligned on 32-byte boundary
x,idx     Must not overlap
N         Multiple of 16 (vthreshold) or 8 (vthresholdf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(int,vthresholdf,( int16_t   * restrict idx, 
                  float32_t            thr,
            const float32_t * restrict x,
            int N ))
#else
int vthresholdf ( int16_t   * restrict idx, 
                  float32_t            thr,
            const float32_t * restrict x,
            int N )
{
  int n, M, cnt;
  unsigned i16 = 0;
  const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
        xb_vecNx16    * restrict pI = (      xb_vecNx16    *)idx;
  xb_vecN_2xf32 x0, x1, t0;
  xb_vecNx16 y0, i0;
  vboolN_2   b0, b1,b2;
  vboolN     bx;
  vselN      sel;
  valign     idx_align;

  if (N <= 0) { return(0); }
  NASSERT_ALIGN(idx, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N>0 && N % (BBE_SIMD_WIDTH/2) == 0);
  t0 = thr;
  idx_align = BBE_ZALIGN();cnt=0;
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
    b0 = BBE_OLTN_2XF32(t0, x0);
    b1 = BBE_OLTN_2XF32(t0, x1);
    bx = BBE_JOINBN_2(b1, b0);
    BBE_SQZN(sel, cnt, bx);
    BBE_MOVIDXNX16T(i0, i0, i16, bx);
    y0 = BBE_SHFLNX16(i0, sel);
    BBE_SAVNX16_XP(y0, idx_align, pI, cnt);
  }
  b2 = BBE_XORBN_2(b2, b2); i0 = BBE_ZERONX16();
  if (N&(BBE_SIMD_WIDTH/2)) 
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b0 = BBE_OLTN_2XF32(t0, x0);
    bx = BBE_JOINBN_2(b2, b0);
    BBE_SQZN(sel, cnt, bx);
    BBE_MOVIDXNX16T(i0, i0, i16, bx);
    y0 = BBE_SHFLNX16(i0, sel);
    BBE_SAVNX16_XP(y0, idx_align, pI, cnt);
  }
  BBE_SAVNX16POS_FP(idx_align, pI);
  M = ((int16_t *)pI - idx);
  return M;
} /* vthresholdf() */
#endif

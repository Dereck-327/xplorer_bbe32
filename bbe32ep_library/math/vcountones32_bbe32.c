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
    Count One Bits in a Word
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Tables for 16-bit fixed-point countones(x) functions. */
#include "countones_16b_tbl.h" 
/*-------------------------------------------------------------------------
Count One Bits in a Word

Description: Functions count the number of one bits in the number or each number of a vector.

Data format: 16-bit/32-bit fixed-point format

Accuracy: exact

Parameters:
Input:
x[N]   Input data, 16-bit/32-bit
N      Length of input/output data vectors
Output:
z[N]   Results, 16-bit/32-bit

Restrictions:
z,x,y  Aligned on 32-byte boundary
z,x,y  Must not overlap
N      Multiple of 16 (vcountones16) or 8 (vcountones32)
-------------------------------------------------------------------------*/

void vcountones32 (int32_t *z, const int32_t   *x, int N)
{
  int n;
  const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict TBL;
  xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
  xb_vecNx16 x0, x1, tbl0, t0, t1;
  vselN u0, u1, u2, u3;
  vselN u4, u5, u6, u7;
  xb_vecNx16 y0, y1, y2, y3;
  xb_vecNx16 y4, y5, y6, y7;
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH/2) == 0);
  if (N <= 0) return;

  TBL = (const xb_vecNx16*)countones_16b_tbl;
  BBE_LVNX16_XP(tbl0, TBL, 2 * BBE_SIMD_WIDTH);
  for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, pX, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(t0, t1, x0, x1, BBE_DSELI_DEINTERLEAVE_1);
    u0 = BBE_MOVVSELNX16(t0, 0);
    u1 = BBE_MOVVSELNX16(t0, 4);
    u2 = BBE_MOVVSELNX16(t0, 8);
    u3 = BBE_MOVVSELNX16(t0, 12);
    u4 = BBE_MOVVSELNX16(t1, 0);
    u5 = BBE_MOVVSELNX16(t1, 4);
    u6 = BBE_MOVVSELNX16(t1, 8);
    u7 = BBE_MOVVSELNX16(t1, 12);
    y0 = BBE_SHFLNX16(tbl0, u0);
    y1 = BBE_SHFLNX16(tbl0, u1);
    y2 = BBE_SHFLNX16(tbl0, u2);
    y3 = BBE_SHFLNX16(tbl0, u3);
    y4 = BBE_SHFLNX16(tbl0, u4);
    y5 = BBE_SHFLNX16(tbl0, u5);
    y6 = BBE_SHFLNX16(tbl0, u6);
    y7 = BBE_SHFLNX16(tbl0, u7);
    y0 = BBE_ADDNX16(y0, y1);
    y2 = BBE_ADDNX16(y2, y3);
    y0 = BBE_ADDNX16(y0, y2);
    y4 = BBE_ADDNX16(y4, y5);
    y6 = BBE_ADDNX16(y6, y7);
    y4 = BBE_ADDNX16(y4, y6);
    y0 = BBE_ADDNX16(y0, y4);
    BBE_DSELNX16I(t0, t1, BBE_ZERONX16(), y0, BBE_DSELI_INTERLEAVE_1);
    BBE_SVNX16_XP(t0, pZ, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(t1, pZ, 2 * BBE_SIMD_WIDTH);
  }
  if (N&(BBE_SIMD_WIDTH/2))
  {
    BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(t0, t1, x0, BBE_ZERONX16(), BBE_DSELI_DEINTERLEAVE_1);
    u0 = BBE_MOVVSELNX16(t0, 0);
    u1 = BBE_MOVVSELNX16(t0, 4);
    u2 = BBE_MOVVSELNX16(t0, 8);
    u3 = BBE_MOVVSELNX16(t0, 12);
    u4 = BBE_MOVVSELNX16(t1, 0);
    u5 = BBE_MOVVSELNX16(t1, 4);
    u6 = BBE_MOVVSELNX16(t1, 8);
    u7 = BBE_MOVVSELNX16(t1, 12);
    y0 = BBE_SHFLNX16(tbl0, u0);
    y1 = BBE_SHFLNX16(tbl0, u1);
    y2 = BBE_SHFLNX16(tbl0, u2);
    y3 = BBE_SHFLNX16(tbl0, u3);
    y4 = BBE_SHFLNX16(tbl0, u4);
    y5 = BBE_SHFLNX16(tbl0, u5);
    y6 = BBE_SHFLNX16(tbl0, u6);
    y7 = BBE_SHFLNX16(tbl0, u7);
    y0 = BBE_ADDNX16(y0, y1);
    y2 = BBE_ADDNX16(y2, y3);
    y0 = BBE_ADDNX16(y0, y2);
    y4 = BBE_ADDNX16(y4, y5);
    y6 = BBE_ADDNX16(y6, y7);
    y4 = BBE_ADDNX16(y4, y6);
    y0 = BBE_ADDNX16(y0, y4);
    BBE_DSELNX16I(t0, t1, BBE_ZERONX16(), y0, BBE_DSELI_INTERLEAVE_1);
    BBE_SVNX16_XP(t0, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* vcountones32() */

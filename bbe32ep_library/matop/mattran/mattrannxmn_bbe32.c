/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Matrix Operations
    Real Matrix Transpose
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"


/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

/* get allocated space per one real/complex matrix written in the block order */
static int getSpace(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, 4);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = (((S - 1) >> m) + 1) << m;
  return S;
} /* getSpace() */
/*-------------------------------------------------------------------------
Real Matrix Transpose

Description: These functions perform transposition for each matrix from input
sequence and store results to output sequence. Both the block order and
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
The functions mattrannxnn(), mattrannxmn() and mattrannxmnf() (real matrix
transpose for the block order) may distort the input matrices sequence x[L*S].

Parameters:
Input:
x[L*S]  Sequence of input matrices
N,M     Matrix dimensions 
L       Number of matrices
Output:
y[L*S]  Sequence of output matrices

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void mattrannxmn ( int16_t * restrict y, 
                   int16_t * restrict x, 
                   int N, int M, int L )
{
  static const ALIGN(32) int16_t sel_tab[BBE_SIMD_WIDTH] =
  { 0, 4, 16, 20,
  1, 5, 17, 21,
  2, 6, 18, 22,
  3, 7, 19, 23
  };

  int i, j, k;
  int Sx = getSpace(M*N);
  xb_vecNx16 x0, x1, x2, x3, t0, t1;
  xb_vecNx16 * restrict py0;
  xb_vecNx16 * restrict py1;
  xb_vecNx16 * restrict py2;
  xb_vecNx16 * restrict py3;
  xb_vecNx16 * restrict px0;
  xb_vecNx16 * restrict px1;
  xb_vecNx16 * restrict px2;
  xb_vecNx16 * restrict px3;

  vselN  s0, s1, s2, s3;
  valign v0, v1, v2, v3;

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(N % 4 == 0);
  NASSERT(M % 4 == 0);
  if (L<=0 || N<=0 || M<=0) return;
  /*--------------------------------------------------
  at the first stage, convert maximum possible portion
  to the streaming format, make permutation and convert back
  --------------------------------------------------*/
  if (L >= BBE_SIMD_WIDTH)
  {
    int _L = L&(~(BBE_SIMD_WIDTH - 1));
    rbsmxn(y, x, M, N, _L);
    mattrannxms(x, y, N, M, _L);
    rsbmxn(y, x, M, N, _L);
    x += Sx*_L;
    y += Sx*_L;
    L &= (BBE_SIMD_WIDTH - 1);
  }
  /*--------------------------------------------------
  last stage: process remainder - up to
  BBE_SIMD_WIDTH matrices
  --------------------------------------------------*/
  x0 = BBE_LVNX16_I((const xb_vecNx16*)sel_tab, 0 * 2 * BBE_SIMD_WIDTH);
  x1 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_4);
  x2 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_8);
  x3 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_12);

  s0 = BBE_MOVVSELNX16(x0, 0);
  s1 = BBE_MOVVSELNX16(x1, 0);
  s2 = BBE_MOVVSELNX16(x2, 0);
  s3 = BBE_MOVVSELNX16(x3, 0);

  for (k = 0; k < L; k++)
  {
    for (j = 0; j < (N >> 2); j++)
    {
      py0 = (xb_vecNx16 *)y;
      py1 = (xb_vecNx16 *)XT_ADDX2(M, (uintptr_t)py0);
      py2 = (xb_vecNx16 *)XT_ADDX4(M, (uintptr_t)py0);
      py3 = (xb_vecNx16 *)XT_ADDX2(M, (uintptr_t)py2);
      v0 = v1 = v2 = v3 = BBE_ZALIGN();
      px0 = (xb_vecNx16 *)x;
      px1 = (xb_vecNx16 *)XT_ADDX2(N, (uintptr_t)px0);
      px2 = (xb_vecNx16 *)XT_ADDX4(N, (uintptr_t)px0);
      px3 = (xb_vecNx16 *)XT_ADDX2(N, (uintptr_t)px2);

      __Pragma("loop_count min=1");
      for (i = 0; i < (M >> 2); i++)
      {
      
        x0 = BBE_LV4X16_I(px0, 0);
        x1 = BBE_LV4X16_I(px1, 0);
        x2 = BBE_LV4X16_I(px2, 0); 
        x3 = BBE_LV4X16_I(px3, 0);
        px0 = (xb_vecNx16 *)XT_ADDX8(N, (uintptr_t)px0);
        px1 = (xb_vecNx16 *)XT_ADDX8(N, (uintptr_t)px1);
        px2 = (xb_vecNx16 *)XT_ADDX8(N, (uintptr_t)px2);
        px3 = (xb_vecNx16 *)XT_ADDX8(N, (uintptr_t)px3);

        t0 = BBE_SELNX16I(x1, x0, BBE_SELI_INTERLEAVE_4_LO);
        t1 = BBE_SELNX16I(x3, x2, BBE_SELI_INTERLEAVE_4_LO);

        x0 = BBE_SELNX16(t1, t0, s0);
        x1 = BBE_SELNX16(t1, t0, s1);
        x2 = BBE_SELNX16(t1, t0, s2);
        x3 = BBE_SELNX16(t1, t0, s3);

        BBE_SAVNX16_XP(x0, v0, py0, 8);
        BBE_SAVNX16_XP(x1, v1, py1, 8);
        BBE_SAVNX16_XP(x2, v2, py2, 8);
        BBE_SAVNX16_XP(x3, v3, py3, 8);
      }
      BBE_SAPOS_FP(v0, py0);
      BBE_SAPOS_FP(v1, py1);
      BBE_SAPOS_FP(v2, py2);
      BBE_SAPOS_FP(v3, py3);
      y = (int16_t *)XT_ADDX8(M, (uintptr_t)y);
      x = (int16_t *)XT_ADDI_N((uintptr_t)x, 8);
    }
    x = (int16_t *)XT_ADDX2(Sx - N, (uintptr_t)x);
    y = (int16_t *)XT_ADDX2(Sx - M*N, (uintptr_t)y);
  }
} /* mattrannxmn() */

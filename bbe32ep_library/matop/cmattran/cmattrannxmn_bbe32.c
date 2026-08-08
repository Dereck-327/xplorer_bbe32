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
    Complex Matrix Conjugate Transpose
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"
/*-------------------------------------------------------------------------
Complex Matrix Conjugate Transpose

Description: These functions perform transposition and then take the complex
conjugate for each matrix from input sequence. Results are stored to output
sequence. Both the block order and streaming order are allowed for input/output
matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
1. Complex conjugation of fixed-point data may involve 16-bit saturation of
   imaginary components
2. The functions cmattrannxnn(), cmattrannxmn() and cmattrannxmnf() (conjugate
   transpose for the block order) may distort the input matrices sequence x[L*S].

Parameters:
Input:
x[L*S]  Sequence of input matrices.
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
void cmattrannxmn ( complex_fract16 * restrict y, 
                    complex_fract16 * restrict x, 
                    int N, int M, int L )
{
  static const int16_t ALIGN(32) sel_tab[2 * BBE_SIMD_WIDTH] =
  {
    0 * 2, 0 * 2 + 1, 4 * 2, 4 * 2 + 1,
    8 * 2, 8 * 2 + 1, 12 * 2, 12 * 2 + 1,
    1 * 2, 1 * 2 + 1, 5 * 2, 5 * 2 + 1,
    9 * 2, 9 * 2 + 1, 13 * 2, 13 * 2 + 1,

    2 * 2, 2 * 2 + 1, 6 * 2, 6 * 2 + 1,
    10 * 2, 10 * 2 + 1, 14 * 2, 14 * 2 + 1,
    3 * 2, 3 * 2 + 1, 7 * 2, 7 * 2 + 1,
    11 * 2, 11 * 2 + 1, 15 * 2, 15 * 2 + 1
  };

  int i, j, k;

  xb_vecNx16 x0, x1, x2, x3, t0, t1;
  xb_vecNx16 * restrict py0;
  xb_vecNx16 * restrict py1;
  xb_vecNx16 * restrict py2;
  xb_vecNx16 * restrict py3;
  xb_vecNx16 * restrict px0;
  xb_vecNx16 * restrict px1;
  xb_vecNx16 * restrict px2;
  xb_vecNx16 * restrict px3;

  vselN s0, s1, s2, s3;
  valign vy0, vy1, vy2, vy3;

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(N % 4 == 0);
  NASSERT(M % 4 == 0);
  if (L <= 0 || N <= 0 || M <= 0) return;
  /*--------------------------------------------------
  at the first stage, convert maximum possible portion
  to the streaming format, make permutation and convert back
  --------------------------------------------------*/
  if (L >= BBE_SIMD_WIDTH / 2)
  {
    int _L = L&(~(BBE_SIMD_WIDTH / 2 - 1)); 
    cbsmxn_even(y, x, M*N, _L);
    cmattrannxms(x, y, N, M, _L);
    csbmxn_even(y, x, M*N, _L);
    x += M*N*_L;
    y += M*N*_L;
    L &= (BBE_SIMD_WIDTH / 2 - 1);
  }
  /*--------------------------------------------------
  last stage: process remainder - up to
  BBE_SIMD_WIDTH/2 matrices
  --------------------------------------------------*/
  x0 = BBE_LVNX16_I((const xb_vecNx16*)sel_tab, 0 * 2 * BBE_SIMD_WIDTH);
  x1 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_8);
  x2 = BBE_LVNX16_I((const xb_vecNx16*)sel_tab, 1 * 2 * BBE_SIMD_WIDTH);
  x3 = BBE_SELNX16I(x2, x2, BBE_SELI_ROTATE_RIGHT_8);

  s0 = BBE_MOVVSELNX16(x0, 0);
  s1 = BBE_MOVVSELNX16(x1, 0);
  s2 = BBE_MOVVSELNX16(x2, 0);
  s3 = BBE_MOVVSELNX16(x3, 0);
  if (L <= 0) return;

  for (k = 0; k<L; k++)
  {
    __Pragma("loop_count min=1");
    for (j = 0; j < (N >> 2); j++)
    {
      py0 = (xb_vecNx16 *)y;
      py1 = (xb_vecNx16 *)XT_ADDX2(2 * M, (uintptr_t)py0);
      py2 = (xb_vecNx16 *)XT_ADDX4(2 * M, (uintptr_t)py0);
      py3 = (xb_vecNx16 *)XT_ADDX2(2 * M, (uintptr_t)py2);

      vy0 = BBE_ZALIGN();
      vy1 = BBE_ZALIGN();
      vy2 = BBE_ZALIGN();
      vy3 = BBE_ZALIGN();
      px0 = (xb_vecNx16 *)x;
      px1 = (xb_vecNx16 *)XT_ADDX2(2 * N, (uintptr_t)px0);
      px2 = (xb_vecNx16 *)XT_ADDX4(2 * N, (uintptr_t)px0);
      px3 = (xb_vecNx16 *)XT_ADDX2(2 * N, (uintptr_t)px2);
      __Pragma("loop_count min=1");
      for (i = 0; i < (M >> 2); i++)
      {
        t0 = BBE_LV4X16_I(px0, 0);
        t1 = BBE_LV4X16_I(px0, 8);
        x0 = BBE_SELNX16I(t1, t0, BBE_SELI_INTERLEAVE_4_EVEN);

        t0 = BBE_LV4X16_I(px1, 0);
        t1 = BBE_LV4X16_I(px1, 8);
        x1 = BBE_SELNX16I(t1, t0, BBE_SELI_INTERLEAVE_4_EVEN);

        t0 = BBE_LV4X16_I(px2, 0);
        t1 = BBE_LV4X16_I(px2, 8);
        x2 = BBE_SELNX16I(t1, t0, BBE_SELI_INTERLEAVE_4_EVEN);

        t0 = BBE_LV4X16_I(px3, 0);
        t1 = BBE_LV4X16_I(px3, 8);
        x3 = BBE_SELNX16I(t1, t0, BBE_SELI_INTERLEAVE_4_EVEN);

        px0 = (xb_vecNx16 *)XT_ADDX8(2 * N, (uintptr_t)px0);
        px1 = (xb_vecNx16 *)XT_ADDX8(2 * N, (uintptr_t)px1);
        px2 = (xb_vecNx16 *)XT_ADDX8(2 * N, (uintptr_t)px2);
        px3 = (xb_vecNx16 *)XT_ADDX8(2 * N, (uintptr_t)px3);

        t0 = BBE_SELNX16I(x1, x0, BBE_SELI_PACK_8);
        t1 = BBE_SELNX16I(x3, x2, BBE_SELI_PACK_8);
 
        t0 = BBE_CONJSNX16C(t0);
        t1 = BBE_CONJSNX16C(t1);

        x0 = BBE_SELNX16(t1, t0, s0);
        x1 = BBE_SELNX16(t1, t0, s1);
        x2 = BBE_SELNX16(t1, t0, s2);
        x3 = BBE_SELNX16(t1, t0, s3);

        BBE_SAVNX16_XP(x0, vy0, py0, 16);
        BBE_SAVNX16_XP(x1, vy1, py1, 16);
        BBE_SAVNX16_XP(x2, vy2, py2, 16);
        BBE_SAVNX16_XP(x3, vy3, py3, 16);
      }
      BBE_SAPOS_FP(vy0, py0);
      BBE_SAPOS_FP(vy1, py1);
      BBE_SAPOS_FP(vy2, py2);
      BBE_SAPOS_FP(vy3, py3);
      y = (complex_fract16 *)XT_ADDX8(2 * M, (uintptr_t)y);
      x += 4;
    }
    x = (complex_fract16 *)XT_ADDX4(M*N - N, (uintptr_t)x);
  }
} /* cmattrannxmn() */

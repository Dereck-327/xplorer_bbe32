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

/* Streaming Order, MxN->NxM, S=MxN
   Restrictions:
     L must be a multiple of 16
*/
void mattrannxms ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int N, int M, int L )
{
  int i, j, k;
  int _M, _N;
  xb_vecNx16 x00;

  const xb_vecNx16 * restrict px = (const xb_vecNx16*)x;
  const xb_vecNx16 * restrict px1 = (const xb_vecNx16*)x;
  xb_vecNx16 * restrict py = (xb_vecNx16*)y;
  xb_vecNx16 * restrict py1 = (xb_vecNx16*)y;

  const int16_t * x0 = x;
  int16_t * y0 = y;

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % BBE_SIMD_WIDTH == 0);
  if (L <= 0 || N<=0 || M<=0 ) return;
  _M = M&~1;
  _N = N&~3;

  /* main square: (M&~1)x(N&~3) */
  for (i = 0; i<(M >> 1); i++)
  {
    for (j = 0; j<(N >> 2); j++)
    {
      px = (const xb_vecNx16 *)(x);
      py = (xb_vecNx16 *)(y);
      px1 = (const xb_vecNx16 *)XT_ADDX2((_M*N) / 2 * L, (uintptr_t)px);
      py1 = (xb_vecNx16 *)XT_ADDX2((_M / 2 * L), (uintptr_t)py);
      k = L >> LOG2_BBE_SIMD_WIDTH;
      do
      {
        x00 = BBE_LVNX16_X(px, 2 * L); BBE_SVNX16_X(x00, py, 2 * M*L);
        x00 = BBE_LVNX16_X(px, 2 * 2 * L); BBE_SVNX16_X(x00, py, 2 * 2 * M*L);
        x00 = BBE_LVNX16_X(px, 3 * 2 * L); BBE_SVNX16_X(x00, py, 3 * 2 * M*L);
        BBE_LVNX16_IP(x00, px, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_IP(x00, py, 2 * BBE_SIMD_WIDTH);

        x00 = BBE_LVNX16_X(px1, 2 * L); BBE_SVNX16_X(x00, py1, 2 * M*L);
        x00 = BBE_LVNX16_X(px1, 2 * 2 * L); BBE_SVNX16_X(x00, py1, 2 * 2 * M*L);
        x00 = BBE_LVNX16_X(px1, 3 * 2 * L); BBE_SVNX16_X(x00, py1, 3 * 2 * M*L);
        BBE_LVNX16_IP(x00, px1, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_IP(x00, py1, 2 * BBE_SIMD_WIDTH);
      } while (--k);
      x = (const int16_t*)XT_ADDX4((2 * L), (uintptr_t)x);
      y = (int16_t*)XT_ADDX4((2 * M*L), (uintptr_t)y);
    }
    x += (N & 3)*L;
    y += (1 - _N*M)*L;
  }
  /* low border (M&1)rows x (N&~3) columns */
  if (M & 1)
  {
    x = (int16_t*)XT_ADDX2((_M*N*L), (uintptr_t)x0);
    y = (int16_t*)XT_ADDX2((_M*L), (uintptr_t)y0);
    for (j = 0; j<N; j++)
    {
      px = (const xb_vecNx16 *)(x);
      py = (xb_vecNx16 *)(y);
      k = L >> LOG2_BBE_SIMD_WIDTH;
      do
      {
        BBE_LVNX16_IP(x00, px, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x00, py, 2 * BBE_SIMD_WIDTH);
      } while (--k);
      x = (const int16_t*)XT_ADDX2((L), (uintptr_t)x);
      y = (int16_t*)XT_ADDX2((M*L), (uintptr_t)y);
    }
  }
  /* right border: M rows x (N&3) columns */
  x = (const int16_t*)XT_ADDX2(_N*L, (uintptr_t)x0);
  y = (int16_t*)XT_ADDX2(_N*M*L, (uintptr_t)y0);
  for (j = 0; j<(N & 3); j++)
  {
    for (i = 0; i<M; i++)
    {
      px = (const xb_vecNx16 *)(x);
      py = (xb_vecNx16 *)(y);
      k = L >> LOG2_BBE_SIMD_WIDTH;
      do
      {
        BBE_LVNX16_IP(x00, px, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x00, py, 2 * BBE_SIMD_WIDTH);
      } while (--k);
      x = (const int16_t*)XT_ADDX2((N*L), (uintptr_t)x);
      y = (int16_t*)XT_ADDX2((L), (uintptr_t)y);
    }
    x = (const int16_t*)XT_ADDX2(((1 - M*N)*L), (uintptr_t)x);
  }
} /* mattrannxms() */

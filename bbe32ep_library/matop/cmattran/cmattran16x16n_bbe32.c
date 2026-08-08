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

/* Block Order, 16x16->16x16, S=256
   Restrictions:
     None
*/
void cmattran16x16n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L )
{
  int i, j, k;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
  const xb_vecNx16* restrict px;
  xb_vecNx16* restrict py;

  int wstride, zstride;

  static const int xsteps[] = { 2 * 2 * BBE_SIMD_WIDTH, -15 * 2 * 2 * BBE_SIMD_WIDTH + 2 * BBE_SIMD_WIDTH, 2 * 2 * BBE_SIMD_WIDTH, 2 * BBE_SIMD_WIDTH };

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  if (L <= 0) return;
  px = (const xb_vecNx16 *)(x);
  py = (xb_vecNx16 *)(y);
  __Pragma("loop_count min=1")
  __Pragma("ymemory(px)")
  for (i = 3, k = 0; k<4 * L; k++)
  {
    i = BBE_ADDMOD16U(i, ((4) << 16) | 1);
    wstride = xsteps[i];

    zstride = 2 * BBE_SIMD_WIDTH;
    j = XT_AND(i, 1);
    XT_MOVEQZ(zstride, -7 * 2 * 2 * BBE_SIMD_WIDTH + 2 * BBE_SIMD_WIDTH, j);

    BBE_LVNX16_IP(x0, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x6, px, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x7, px, wstride);

    BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x5, x4, x5, x4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x7, x6, x7, x6, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_2);

    x0 = BBE_CONJSNX16C(x0);
    x1 = BBE_CONJSNX16C(x1);
    x2 = BBE_CONJSNX16C(x2);
    x3 = BBE_CONJSNX16C(x3);
    x4 = BBE_CONJSNX16C(x4);
    x5 = BBE_CONJSNX16C(x5);
    x6 = BBE_CONJSNX16C(x6);
    x7 = BBE_CONJSNX16C(x7);

    BBE_SVNX16_IP(x0, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x4, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x5, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x6, py, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(x7, py, zstride);
  }
} /* cmattran16x16n() */

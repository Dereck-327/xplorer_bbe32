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

#if !(HAVE_VFPU)
DISCARD_FUN(void, cmattran8x8sf, ( complex_float * restrict y, 
                             const complex_float * restrict x, 
                                   int L ))
#else

#define sz_cf32 (int)sizeof(complex_float)

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

/* Streaming Order, Floating-Point, 8x8->8x8, S=64
   Restrictions:
     L must be a multiple of 4
*/
void cmattran8x8sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  const int M = 8;
  const int N = 8;
  int i, j, l;
  xb_vecN_4xcf32 vreg;

  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
        xb_vecN_4xcf32 * restrict py0;
        xb_vecN_4xcf32 * restrict py1;

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 4) == 0);

  for (i = 0; i<(M/4); i++)
  {
    for (j = 0; j<(N/2); j++)
    {
      px0 = (const xb_vecN_4xcf32 *)(x);
      py0 = (      xb_vecN_4xcf32 *)(y);
      px1 = (const xb_vecN_4xcf32 *)((complex_float *)px0+L);
      py1 = (      xb_vecN_4xcf32 *)((complex_float *)py0+L*M);

      for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
      {
        vreg = BBE_LVN_4XCF32_X(px0, 1*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py0, 1*L*sz_cf32);
        vreg = BBE_LVN_4XCF32_X(px0, 2*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py0, 2*L*sz_cf32);
        vreg = BBE_LVN_4XCF32_X(px0, 3*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py0, 3*L*sz_cf32);
        BBE_LVN_4XCF32_IP(vreg, px0, 2*BBE_SIMD_WIDTH); vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_IP(vreg, py0, 2*BBE_SIMD_WIDTH);

        vreg = BBE_LVN_4XCF32_X(px1, 1*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py1, 1*L*sz_cf32);
        vreg = BBE_LVN_4XCF32_X(px1, 2*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py1, 2*L*sz_cf32);
        vreg = BBE_LVN_4XCF32_X(px1, 3*N*L*sz_cf32);    vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_X (vreg, py1, 3*L*sz_cf32);
        BBE_LVN_4XCF32_IP(vreg, px1, 2*BBE_SIMD_WIDTH); vreg = BBE_CONJN_4XCF32(vreg);  BBE_SVN_4XCF32_IP(vreg, py1, 2*BBE_SIMD_WIDTH);
      }
      x += L*2;
      y += L*M*2;
    }
    x += L*N*3;
    y += L*4-L*M*N;
  }
} /* cmattran8x8sf() */
#endif

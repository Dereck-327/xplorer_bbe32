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
DISCARD_FUN(void, cmattran4x4nf, ( complex_float * restrict y, 
                             const complex_float * restrict x, 
                                   int L ))
#else

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

/* Block Order, Floating-Point, 4x4->4x4, S=16
   Restrictions:
     None
*/
void cmattran4x4nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
        xb_vecN_4xcf32 * restrict py0;
        xb_vecN_4xcf32 * restrict py1;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3;
  xb_vecN_4xcf32 Y0, Y1, Y2, Y3;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x+8);
  py0 = (      xb_vecN_4xcf32 *)(y);
  py1 = (      xb_vecN_4xcf32 *)(y+4);

  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X0, px0, 1*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 3*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px1, 1*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px1, 3*2*BBE_SIMD_WIDTH);

    /* Transpose and conjugate */
    BBE_DSELN_4XCF32I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_4);
    Y0 = BBE_CONJN_4XCF32(X0);
    Y1 = BBE_CONJN_4XCF32(X1);
    Y2 = BBE_CONJN_4XCF32(X2);
    Y3 = BBE_CONJN_4XCF32(X3);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y0, py0, 2*2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py1, 2*2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y2, py0, 2*2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y3, py1, 2*2*BBE_SIMD_WIDTH);
  }
} /* cmattran4x4nf() */
#endif

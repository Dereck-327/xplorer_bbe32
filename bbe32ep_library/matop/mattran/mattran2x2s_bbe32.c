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

/* Streaming Order, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 16
*/
void mattran2x2s ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L )
{
  xb_vecNx16 x00, x01, x02, x03;

  const xb_vecNx16 * restrict px0 = (const xb_vecNx16*)(x + L * 0);
  const xb_vecNx16 * restrict px1 = (const xb_vecNx16*)(x + L * 1);
  const xb_vecNx16 * restrict px2 = (const xb_vecNx16*)(x + L * 2);
  const xb_vecNx16 * restrict px3 = (const xb_vecNx16*)(x + L * 3);

  xb_vecNx16 * restrict py0 = (xb_vecNx16*)(y + L * 0);
  xb_vecNx16 * restrict py1 = (xb_vecNx16*)(y + L * 1);
  xb_vecNx16 * restrict py2 = (xb_vecNx16*)(y + L * 2);
  xb_vecNx16 * restrict py3 = (xb_vecNx16*)(y + L * 3);

  __Pragma("ymemory( py0 )");
  __Pragma("ymemory( py1 )");
  __Pragma("ymemory( py2 )");
  __Pragma("ymemory( py3 )");

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % BBE_SIMD_WIDTH == 0);
  if (L<=0) return;
  L >>= LOG2_BBE_SIMD_WIDTH;
  do
  {
    BBE_LVNX16_IP(x00, px0, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x01, px1, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x02, px2, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x03, px3, (2 * BBE_SIMD_WIDTH));

    BBE_SVNX16_IP(x00, py0, (2 * BBE_SIMD_WIDTH));
    BBE_SVNX16_IP(x01, py2, (2 * BBE_SIMD_WIDTH));
    BBE_SVNX16_IP(x02, py1, (2 * BBE_SIMD_WIDTH));
    BBE_SVNX16_IP(x03, py3, (2 * BBE_SIMD_WIDTH));
  } while (--L);
} /* mattran2x2s() */

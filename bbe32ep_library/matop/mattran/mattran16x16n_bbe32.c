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

/* Block Order, 16x16->16x16, S=256
   Restrictions:
     None
*/
void mattran16x16n ( int16_t * restrict y, 
               const int16_t * restrict x, 
               int L )
{
  int i;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;
  const xb_vecNx16* restrict px = (const xb_vecNx16 *)x;
  xb_vecNx16* restrict py = (xb_vecNx16 *)y;

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

  for (i = 0; i<L; i++)
  {
    BBE_LVNX16_IP(x0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x4, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x5, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x6, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x7, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x8, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x9, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xa, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xb, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xc, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xd, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xe, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(xf, px, 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x5, x4, x5, x4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x7, x6, x7, x6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x9, x8, x9, x8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xb, xa, xb, xa, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xd, xc, xd, xc, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xf, xe, xf, xe, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xa, x8, xa, x8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xb, x9, xb, x9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xe, xc, xe, xc, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xf, xd, xf, xd, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xc, x8, xc, x8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xd, x9, xd, x9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xe, xa, xe, xa, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xf, xb, xf, xb, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(x8, x0, x8, x0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x9, x1, x9, x1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xa, x2, xa, x2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xb, x3, xb, x3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xc, x4, xc, x4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xd, x5, xd, x5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xe, x6, xe, x6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(xf, x7, xf, x7, BBE_DSELI_DEINTERLEAVE_1);

    BBE_SVNX16_IP(x0, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x1, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x2, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x3, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x4, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x5, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x6, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x7, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x8, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(x9, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xa, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xb, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xc, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xd, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xe, py, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(xf, py, 2 * BBE_SIMD_WIDTH);
  }
} /* mattran16x16n() */

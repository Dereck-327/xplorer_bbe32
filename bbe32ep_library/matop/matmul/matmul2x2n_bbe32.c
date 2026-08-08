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
    Real Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, matmul2x2n,(int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q))
#else
/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 4
*/
void matmul2x2n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int l;

  vsaN  q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 X, Y, Z;
  xb_vecNx16 x0, y0, x1, y1;

  xb_vecNx40 r;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % 4 == 0);

  for (l = 0; l<L; l += (BBE_SIMD_WIDTH / 4))
  {
    /* Load input matrix X and Y */
    BBE_LVNX16_IP(X, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

    x0 = BBE_SHFLNX16I(X, BBE_SHFLI_DUPLICATE_1_EVEN);
    y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M1_STEP_1);

    /* Multiply input matrix X and Y */
    r = BBE_MULRNX16(x0, y0, q);

    x1 = BBE_SHFLNX16I(X, BBE_SHFLI_DUPLICATE_1_ODD);
    y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M1_STEP_2);

    /* Multiply input matrix X and Y */
    BBE_MULANX16(r, x1, y1);

    /* Pack and save rezults */
    Z = BBE_PACKVNX40(r, q);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* matmul2x2n() */
#endif

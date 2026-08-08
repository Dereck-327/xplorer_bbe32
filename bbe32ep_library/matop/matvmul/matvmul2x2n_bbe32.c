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

/* Block Order, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions: 
     L must be a multiple of 8 
*/
void matvmul2x2n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  int l;

  xb_vecNx16 X0, X1, Y, Y0, Y1, Z, Z0, Z1;

  xb_vecNx40 acc0, acc1;

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  vsaN  q = BBE_MOVVSA32(Q);

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L % (BBE_SIMD_WIDTH / 2) == 0);
  if (L<=0) return;
  if (Q == 15)
  {
    __Pragma("ymemory( px )");
    __Pragma("ymemory( py )");
    for (l = 0; l<L; l += 8)
    {
      BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Y1, Y0, Y, Y, BBE_DSELI_INTERLEAVE_2);

      Z0 = BBE_MULNX16JPACKQ(X0, Y0); Z1 = BBE_MULNX16JPACKQ(X1, Y1);

      Z = BBE_SELNX16I(Z1, Z0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);

      BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
    }
  }
  else
  {
    __Pragma("ymemory( px )");
    __Pragma("ymemory( py )");
    for (l = 0; l<L; l += 8)
    {
      BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Y1, Y0, Y, Y, BBE_DSELI_INTERLEAVE_2);

      acc0 = BBE_MULRNX16J(X0, Y0, q); acc1 = BBE_MULRNX16J(X1, Y1, q);

      Z0 = BBE_PACKVNX40(acc0, q); Z1 = BBE_PACKVNX40(acc1, q);
      Z = BBE_SELNX16I(Z1, Z0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);

      BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
    }
  }
} /* matvmul2x2n() */

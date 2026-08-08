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

/* Streaming Order, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 16
*/
void matvmul3x3s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
int l;

  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  xb_vecNx16 vX00, vX01, vX02, vX10, vX11, vX12, vX20, vX21, vX22;
  xb_vecNx16 vY0, vY1, vY2;

  xb_vecNx16 vZ0, vZ1, vZ2;

  xb_vecNx40 Z0, Z1, Z2;

  const vsaN q = BBE_MOVVSA32(Q);
  const int dOffset = L*sizeof(int16_t);

  /* check restrictions */
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT_ALIGN16(z);
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);

  if (L<=0) return; 
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += BBE_SIMD_WIDTH)
  {
    /* Loading 16 input matrices X */
    BBE_LVNX16_XP(vX00, pXrd, dOffset);
    BBE_LVNX16_XP(vX01, pXrd, dOffset);
    BBE_LVNX16_XP(vX02, pXrd, dOffset);

    BBE_LVNX16_XP(vX10, pXrd, dOffset);
    BBE_LVNX16_XP(vX11, pXrd, dOffset);
    BBE_LVNX16_XP(vX12, pXrd, dOffset);

    BBE_LVNX16_XP(vX20, pXrd, dOffset);
    BBE_LVNX16_XP(vX21, pXrd, dOffset);
    BBE_LVNX16_XP(vX22, pXrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);

    /* Loading 16 input vectors Y */
    BBE_LVNX16_XP(vY0, pYrd, dOffset);
    BBE_LVNX16_XP(vY1, pYrd, dOffset);
    BBE_LVNX16_XP(vY2, pYrd, -2 * dOffset + 2 * BBE_SIMD_WIDTH);

    Z0 = BBE_MULRNX16(vX00, vY0, q);
    BBE_MULANX16(Z0, vX01, vY1);
    BBE_MULANX16(Z0, vX02, vY2);

    Z1 = BBE_MULRNX16(vX10, vY0, q);
    BBE_MULANX16(Z1, vX11, vY1);
    BBE_MULANX16(Z1, vX12, vY2);

    Z2 = BBE_MULRNX16(vX20, vY0, q);
    BBE_MULANX16(Z2, vX21, vY1);
    BBE_MULANX16(Z2, vX22, vY2);

    /* Pack and save result */
    vZ0 = BBE_PACKVNX40(Z0, q);
    vZ1 = BBE_PACKVNX40(Z1, q);
    vZ2 = BBE_PACKVNX40(Z2, q);
    BBE_SVNX16_XP(vZ0, pZwr, dOffset);
    BBE_SVNX16_XP(vZ1, pZwr, dOffset);
    BBE_SVNX16_XP(vZ2, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);
  }
} /* matvmul3x3s() */

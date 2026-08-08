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

/* Streaming Order, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 16
*/
void matmul3x3s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int l;

  xb_vecNx16  vX0, vX1, vX2;
  xb_vecNx16  vY;
  xb_vecNx16  vZ;

  xb_vecNx40 ACC;

  const vsaN q = BBE_MOVVSA32(Q);
  const int dOffset = L*sizeof(int16_t);

  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);
  if (L <= 0) return;

  __Pragma("loop_count min=1")
    for (l = 0; l<L; l += BBE_SIMD_WIDTH)
    {
      /* Loading a row of 16 input matrices X */
      BBE_LVNX16_XP(vX0, pXrd, dOffset);
      BBE_LVNX16_XP(vX1, pXrd, dOffset);
      BBE_LVNX16_XP(vX2, pXrd, dOffset);

      /* Multiple pack and save */
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -8 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      /* Loading a row of 16 input matrices X */
      BBE_LVNX16_XP(vX0, pXrd, dOffset);
      BBE_LVNX16_XP(vX1, pXrd, dOffset);
      BBE_LVNX16_XP(vX2, pXrd, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -8 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      /* Loading a row of 16 input matrices X */
      BBE_LVNX16_XP(vX0, pXrd, dOffset);
      BBE_LVNX16_XP(vX1, pXrd, dOffset);
      BBE_LVNX16_XP(vX2, pXrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, dOffset);

      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      ACC = BBE_MULRNX16(vX0, vY, q);
      BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
      BBE_MULANX16(ACC, vX1, vY);
      BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
      BBE_MULANX16(ACC, vX2, vY);
      vZ = BBE_PACKVNX40(ACC, q);
      BBE_SVNX16_XP(vZ, pZwr, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    }
} /* matmul3x3s() */

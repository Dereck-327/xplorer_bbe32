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
    Real Matrix by Complex Matrix/Vector Multiply
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
Real Matrix by Complex Matrix/Vector Multiply 

Description: These functions perform pairwise multiplication of left-hand
real matrices by right-hand complex matrices or vectors. Both the block order
and streaming order are allowed for input/output matrix sequences.

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand real matrices
y[L*Sy]     Sequence of right-hand complex matrices or vectors
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Streaming Order, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul3x3s ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  int l;

  xb_vecNx16  vX0l, vX1l, vX2l;
  xb_vecNx16  vX0h, vX1h, vX2h;
  xb_vecNx16  vXX0, vXX1, vXX2;
  xb_vecNx16  vY;
  xb_vecNx16  vZ;
  const xb_vecNx16 vZero = 0;

  xb_vecNx40 ACC;

  const vsaN q = BBE_MOVVSA32(Q);
  const int dOffset0 = L*sizeof(int16_t);
  const int dOffset = 2 * dOffset0;


  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += BBE_SIMD_WIDTH)
  {
    /*Load input matrix X*/
    BBE_LVNX16_XP(vXX0, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX1, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX2, pXrd, -2 * dOffset0 + 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(vX0h, vX0l, vZero, vXX0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX1h, vX1l, vZero, vXX1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX2h, vX2l, vZero, vXX2, BBE_DSELI_INTERLEAVE_1);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);
  }
  pXrd = (const xb_vecNx16*)XT_ADDX4(L, (uintptr_t)pXrd);
  pYrd = (const xb_vecNx16*)XT_ADDX4(-L, (uintptr_t)pYrd);
  pZwr = (xb_vecNx16*)XT_ADDX8(L, (uintptr_t)pZwr);
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += BBE_SIMD_WIDTH)
  {
    /*Load input matrix X*/
    BBE_LVNX16_XP(vXX0, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX1, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX2, pXrd, -2 * dOffset0 + 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(vX0h, vX0l, vZero, vXX0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX1h, vX1l, vZero, vXX1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX2h, vX2l, vZero, vXX2, BBE_DSELI_INTERLEAVE_1);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);
  }
  pXrd = (const xb_vecNx16*)XT_ADDX4(L, (uintptr_t)pXrd);
  pYrd = (const xb_vecNx16*)XT_ADDX4(-L, (uintptr_t)pYrd);
  pZwr = (xb_vecNx16*)XT_ADDX8(L, (uintptr_t)pZwr);
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += BBE_SIMD_WIDTH)
  {
    /*Load input matrix X*/
    BBE_LVNX16_XP(vXX0, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX1, pXrd, dOffset0);
    BBE_LVNX16_XP(vXX2, pXrd, -2 * dOffset0 + 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(vX0h, vX0l, vZero, vXX0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX1h, vX1l, vZero, vXX1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(vX2h, vX2l, vZero, vXX2, BBE_DSELI_INTERLEAVE_1);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0l, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1l, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2l, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -5 * dOffset);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, dOffset);

    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    ACC = BBE_MULRNX16C(vX0h, vY, q);
    BBE_LVNX16_XP(vY, pYrd, 3 * dOffset);
    BBE_MULANX16C(ACC, vX1h, vY);
    BBE_LVNX16_XP(vY, pYrd, -8 * dOffset + 2 * BBE_SIMD_WIDTH);
    BBE_MULANX16C(ACC, vX2h, vY);
    /* Pack and save results */
    vZ = BBE_PACKVNX40(ACC, q);
    BBE_SVNX16_XP(vZ, pZwr, -2 * dOffset + 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatmul3x3s() */

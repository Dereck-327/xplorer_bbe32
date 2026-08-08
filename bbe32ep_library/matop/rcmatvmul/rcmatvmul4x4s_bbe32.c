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

/* Streaming Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul4x4s ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  int _L0 = -L * 16 + BBE_SIMD_WIDTH;
  int _L1 = L * 16 - 4 * BBE_SIMD_WIDTH;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 Y0_0, Y1_0, Y2_0, Y3_0;
  xb_vecNx16 Y0_1, Y1_1, Y2_1, Y3_1;
  xb_vecNx16 X0, X1, X2, X3;
  xb_vecNx16 X0_0, X1_0, X2_0, X3_0;
  xb_vecNx16 X0_1, X1_1, X2_1, X3_1;

  xb_vecNx16 Z0, Z1;

  xb_vecNx40 Acc0, Acc1;

  xb_vecNx16 zero;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH) == 0);
  if (L <= 0) return;
  zero = 0;

  __Pragma("ymemory( py )");
  __Pragma("loop_count min=2");
  for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
  {
    Y0_1 = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(Y0_0, py, 2 * 2 * L);
    Y1_1 = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(Y1_0, py, 2 * 2 * L);
    Y2_1 = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(Y2_0, py, 2 * 2 * L);
    Y3_1 = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(Y3_0, py, -3 * 2 * 2 * L);

    /*Multiply first string*/

    BBE_LVNX16_XP(X0, px, 2 * L);
    BBE_LVNX16_XP(X1, px, 2 * L);
    BBE_LVNX16_XP(X2, px, 2 * L);
    BBE_LVNX16_XP(X3, px, 2 * L);

    BBE_DSELNX16I(X0_1, X0_0, zero, X0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X1_1, X1_0, zero, X1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X2_1, X2_0, zero, X2, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3_1, X3_0, zero, X3, BBE_DSELI_INTERLEAVE_1);

    Acc0 = BBE_MULRNX16C(Y0_0, X0_0, q); Acc1 = BBE_MULRNX16C(Y0_1, X0_1, q);
    BBE_MULANX16C(Acc0, Y1_0, X1_0);     BBE_MULANX16C(Acc1, Y1_1, X1_1);
    BBE_MULANX16C(Acc0, Y2_0, X2_0);     BBE_MULANX16C(Acc1, Y2_1, X2_1);
    BBE_MULANX16C(Acc0, Y3_0, X3_0);     BBE_MULANX16C(Acc1, Y3_1, X3_1);

    Z0 = BBE_PACKVNX40(Acc0, q); Z1 = BBE_PACKVNX40(Acc1, q);

    BBE_SVNX16_I(Z1, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(Z0, pz, 2 * 2 * L);

    /*Multiply second string*/

    BBE_LVNX16_XP(X0, px, 2 * L);
    BBE_LVNX16_XP(X1, px, 2 * L);
    BBE_LVNX16_XP(X2, px, 2 * L);
    BBE_LVNX16_XP(X3, px, 2 * L);

    BBE_DSELNX16I(X0_1, X0_0, zero, X0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X1_1, X1_0, zero, X1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X2_1, X2_0, zero, X2, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3_1, X3_0, zero, X3, BBE_DSELI_INTERLEAVE_1);

    Acc0 = BBE_MULRNX16C(Y0_0, X0_0, q); Acc1 = BBE_MULRNX16C(Y0_1, X0_1, q);
    BBE_MULANX16C(Acc0, Y1_0, X1_0);     BBE_MULANX16C(Acc1, Y1_1, X1_1);
    BBE_MULANX16C(Acc0, Y2_0, X2_0);     BBE_MULANX16C(Acc1, Y2_1, X2_1);
    BBE_MULANX16C(Acc0, Y3_0, X3_0);     BBE_MULANX16C(Acc1, Y3_1, X3_1);

    Z0 = BBE_PACKVNX40(Acc0, q); Z1 = BBE_PACKVNX40(Acc1, q);

    BBE_SVNX16_I(Z1, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(Z0, pz, 2 * 2 * L);

    px = (const xb_vecNx16 *)XT_ADDX2(_L0*(l & 1), (int32_t)px);
    py = (const xb_vecNx16 *)XT_ADDX4(BBE_SIMD_WIDTH*(l & 1), (int32_t)py);
    pz = (xb_vecNx16 *)XT_SUB((int32_t)pz, _L1*(l & 1));
  }
} /* rcmatvmul4x4s() */

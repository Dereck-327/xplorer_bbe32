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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, rcmatvmul4x4n,(complex_fract16 * restrict z, 
            const int16_t * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
#else
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

/* Block Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void rcmatvmul4x4n ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  vsaN q = BBE_MOVVSA32(Q);

  static const int16_t ALIGN(32) sel0[16] = { 0, 16, 1, 17, 4, 18, 5, 19, 8, 20, 9, 21, 12, 22, 13, 23 };
  static const int16_t ALIGN(32) sel1[16] = { 2, 16, 3, 17, 6, 18, 7, 19, 10, 20, 11, 21, 14, 22, 15, 23 };

  const xb_vecNx16 *pSEL0 = (const xb_vecNx16 *)sel0;
  const xb_vecNx16 *pSEL1 = (const xb_vecNx16 *)sel1;

  xb_vecNx16 Sel0, Sel1;

  vselN vSel0, vSel1;

  xb_vecNx16 X0, X1, Y, X_, Y_, Z_;

  xb_vecNx16 zero = 0;

  xb_vecNx40 Acc;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % (BBE_SIMD_WIDTH / 8) == 0);

  if (L<=0) return;
  Sel0 = BBE_LVNX16_I(pSEL0, 0);
  vSel0 = BBE_MOVVSV(Sel0, 0);

  Sel1 = BBE_LVNX16_I(pSEL1, 0);
  vSel1 = BBE_MOVVSV(Sel1, 0);

  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  __Pragma("loop_count min=1");
  for (l = 0; l<L; l += 2)
  {
    /* Load input matrix X */
    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

    Y_ = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_LOW_HALF);
    X_ = BBE_SELNX16(zero, X0, vSel0);
    Acc = BBE_MULRNX16PC_0(X_, Y_, q);

    Y_ = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_LOW_HALF);
    X_ = BBE_SELNX16(zero, X0, vSel1);
    BBE_MULANX16PC_0(Acc, X_, Y_);

    Y_ = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_HIGH_HALF);
    X_ = BBE_SELNX16(zero, X1, vSel0);
    BBE_MULANX16PC_1(Acc, X_, Y_);

    Y_ = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_HIGH_HALF);
    X_ = BBE_SELNX16(zero, X1, vSel1);
    BBE_MULANX16PC_1(Acc, X_, Y_);

    /* Pack and save results */
    Acc = BBE_SHFLNX40I(Acc, BBE_W_SHFLI_DITLV_2);
    Z_ = BBE_PACKVNX40(Acc, q);
    BBE_SVNX16_IP(Z_, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul4x4n() */
#endif

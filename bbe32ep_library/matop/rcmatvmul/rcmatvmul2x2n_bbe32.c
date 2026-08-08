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
#if !(HAVE_MULPC && 1)
DISCARD_FUN(void, rcmatvmul2x2n,(complex_fract16 * restrict z, 
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

/* Block Order, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 4
*/
void rcmatvmul2x2n ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  xb_vecNx16 zero = 0;

  xb_vecNx16 X0, Y0, Z0, X_l0, X_h0;
  xb_vecNx16 X1, Y1, Z1, X_l1, X_h1;

  vsaN q = BBE_MOVVSA32(Q);

  static const int16_t ALIGN(32) sel[16] = { 0, 16, 1, 17, 4, 18, 5, 19, 8, 20, 9, 21, 12, 22, 13, 23 };

  const xb_vecNx16 *pSEL = (const xb_vecNx16 *)sel;

  xb_vecNx16 Sel;

  xb_vecNx40 Acc0, Acc1;

  vselN vSel;

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % (BBE_SIMD_WIDTH / 4) == 0);

  if (L<=0) return;
  Sel = BBE_LVNX16_I(pSEL, 0);
  vSel = BBE_MOVVSV(Sel, 0);

  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  for (l = 0; l<L; l += 8)
  {
    /*Load input matrix X*/
    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);

    /*Load input matrix Y*/
    BBE_LVNX16_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, py, 2 * BBE_SIMD_WIDTH);

    X_l1 = BBE_SELNX16(zero, X1, vSel);
    BBE_SELUNX16(X_l0, zero, X0, vSel, 2);
    X_h1 = BBE_SELNX16(zero, X1, vSel);
    BBE_SELUNX16(X_h0, zero, X0, vSel, 126);

    Acc0 = BBE_MULRNX16PC_0(X_l0, Y0, q);
    Acc1 = BBE_MULRNX16PC_0(X_l1, Y1, q);
    BBE_MULANX16PC_1(Acc0, X_h0, Y0);
    BBE_MULANX16PC_1(Acc1, X_h1, Y1);

    /* Pack and save results */
    Z0 = BBE_PACKVNX40(Acc0, q);
    Z1 = BBE_PACKVNX40(Acc1, q);
    BBE_SVNX16_IP(Z0, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Z1, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul2x2n() */
#endif

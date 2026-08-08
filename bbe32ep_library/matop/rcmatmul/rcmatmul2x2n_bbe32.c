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
DISCARD_FUN(void ,rcmatmul2x2n,(complex_fract16 * restrict z, 
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

/* Block Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void rcmatmul2x2n ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  int l;

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  vsaN q = BBE_MOVVSA32(Q);

  static const int16_t ALIGN(32) sel[16] = { 0, 1, 4, 5, 0, 1, 4, 5, 8, 9, 12, 13, 8, 9, 12, 13 };

  const xb_vecNx16 *pSEL = (const xb_vecNx16 *)sel;

  xb_vecNx16 Sel;

  vselN vSel0, vSel1;

  xb_vecNx16 X, X0, X1, Y, Y0, Y1, Z0, Z1;

  xb_vecNx16 zero = 0;

  xb_vecNx40 Acc0, Acc1;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L % (BBE_SIMD_WIDTH / 4) == 0);
  if (L<=0) return;
  Sel = BBE_LVNX16_I(pSEL, 0);
  vSel1 = vSel0 = BBE_MOVVSV(Sel, 0);

  BBE_SELUNX16(X, X0, X1, vSel1, 2);

  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  __Pragma("loop_count min=1");
  for (l = 0; l<L; l += 4)
  {
    /* Load input matrix X */
    BBE_LVNX16_IP(X, px, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X1, X0, zero, X, BBE_DSELI_INTERLEAVE_1);

    /* Load input matrix Y */
    BBE_LVNX16_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, py, 2 * BBE_SIMD_WIDTH);

    Y = BBE_SHFLNX16(Y0, vSel0);
    Acc0 = BBE_MULRNX16PC_0(Y, X0, q);

    Y = BBE_SHFLNX16(Y0, vSel1);
    BBE_MULANX16PC_1(Acc0, Y, X0);

    /* Pack and save results */
    Z0 = BBE_PACKVNX40(Acc0, q);
    BBE_SVNX16_IP(Z0, pz, 2 * BBE_SIMD_WIDTH);

    Y = BBE_SHFLNX16(Y1, vSel0);
    Acc1 = BBE_MULRNX16PC_0(Y, X1, q);

    Y = BBE_SHFLNX16(Y1, vSel1);
    BBE_MULANX16PC_1(Acc1, Y, X1);

    /* Pack and save results */
    Z1 = BBE_PACKVNX40(Acc1, q);
    BBE_SVNX16_IP(Z1, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatmul2x2n() */
#endif


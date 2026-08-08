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

/* Streaming Order, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul3x3s ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  xb_vecNx16 vXX, vXl, vXh;
  xb_vecNx16 vY0l, vY1l, vY2l;
  xb_vecNx16 vY0h, vY1h, vY2h;

  xb_vecNx16 vZl, vZh;

  const xb_vecNx16 vZero = 0;

  xb_vecNx40 Zl, Zh;

  const vsaN q = BBE_MOVVSA32(Q);
  const int dOffset0 = L*sizeof(int16_t);
  const int dOffset = 2 * dOffset0;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(L % (BBE_SIMD_WIDTH) == 0);
  NASSERT(Q >= 0 && Q <= 16);
  if (L <= 0) return;

  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += BBE_SIMD_WIDTH)
  {
    /* Load 8 input vectors Y */
    vY0h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY0l, pYrd, dOffset);
    vY1h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY1l, pYrd, dOffset);
    vY2h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY2l, pYrd, -dOffset);

    /* load a row of 8 X matrices */
    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    Zl = BBE_MULRNX16C(vXl, vY0l, q);
    Zh = BBE_MULRNX16C(vXh, vY0h, q);

    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY1l);
    BBE_MULANX16C(Zh, vXh, vY1h);

    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY2l);
    BBE_MULANX16C(Zh, vXh, vY2h);

    vZl = BBE_PACKVNX40(Zl, q);
    vZh = BBE_PACKVNX40(Zh, q);
    BBE_SVNX16_I(vZh, pZwr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(vZl, pZwr, dOffset);

    /* Reload 8 input vectors Y */
    vY1h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY1l, pYrd, dOffset);
    vY2h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY2l, pYrd, -dOffset);

    /* load a row of 8 X matrices */
    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    Zl = BBE_MULRNX16C(vXl, vY0l, q);
    Zh = BBE_MULRNX16C(vXh, vY0h, q);

    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY1l);
    BBE_MULANX16C(Zh, vXh, vY1h);

    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY2l);
    BBE_MULANX16C(Zh, vXh, vY2h);

    vZl = BBE_PACKVNX40(Zl, q);
    vZh = BBE_PACKVNX40(Zh, q);
    BBE_SVNX16_I(vZh, pZwr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(vZl, pZwr, dOffset);

    /* Reload 8 input vectors Y */
    vY1h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY1l, pYrd, dOffset);
    vY2h = BBE_LVNX16_I(pYrd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(vY2l, pYrd, -2 * dOffset + 4 * BBE_SIMD_WIDTH);

    /* load a row of 8 X matrices */
    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    Zl = BBE_MULRNX16C(vXl, vY0l, q);
    Zh = BBE_MULRNX16C(vXh, vY0h, q);

    BBE_LVNX16_XP(vXX, pXrd, dOffset0);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY1l);
    BBE_MULANX16C(Zh, vXh, vY1h);

    BBE_LVNX16_XP(vXX, pXrd, -8 * dOffset0 + 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(vXh, vXl, vZero, vXX, BBE_DSELI_INTERLEAVE_1);
    BBE_MULANX16C(Zl, vXl, vY2l);
    BBE_MULANX16C(Zh, vXh, vY2h);

    vZl = BBE_PACKVNX40(Zl, q);
    vZh = BBE_PACKVNX40(Zh, q);
    BBE_SVNX16_I(vZh, pZwr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(vZl, pZwr, -2 * dOffset + 4 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul3x3s() */

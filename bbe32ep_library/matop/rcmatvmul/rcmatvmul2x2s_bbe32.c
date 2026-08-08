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

/* Streaming Order, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul2x2s ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  xb_vecNx16 xx00, xx01, xx10, xx11, zero;
  xb_vecNx16 x00l, x00h, x01l, x01h, x10l, x10h, x11l, x11h;
  xb_vecNx16 y0l, y0h, y1l, y1h;
  xb_vecNx16 z0l, z0h, z1l, z1h;

  xb_vecNx40 Zl, Zh;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH) == 0);
  if (L <= 0) return;
  zero = 0;

  __Pragma("ymemory( px )");
  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    BBE_LVNX16_XP(xx00, px, 2 * L);
    BBE_LVNX16_XP(xx01, px, 2 * L);
    BBE_LVNX16_XP(xx10, px, 2 * L);
    BBE_LVNX16_XP(xx11, px, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

    y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y0l, py, 2 * 2 * L);
    y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1l, py, -2 * 2 * L + 2 * 2 * BBE_SIMD_WIDTH);

    x00l = BBE_SELNX16I(zero, xx00, BBE_SELI_INTERLEAVE_1_LO);
    x00h = BBE_SELNX16I(zero, xx00, BBE_SELI_INTERLEAVE_1_HI);
    x01l = BBE_SELNX16I(zero, xx01, BBE_SELI_INTERLEAVE_1_LO);
    x01h = BBE_SELNX16I(zero, xx01, BBE_SELI_INTERLEAVE_1_HI);
    x10l = BBE_SELNX16I(zero, xx10, BBE_SELI_INTERLEAVE_1_LO);
    x10h = BBE_SELNX16I(zero, xx10, BBE_SELI_INTERLEAVE_1_HI);
    x11l = BBE_SELNX16I(zero, xx11, BBE_SELI_INTERLEAVE_1_LO);
    x11h = BBE_SELNX16I(zero, xx11, BBE_SELI_INTERLEAVE_1_HI);

    Zl = BBE_MULRNX16C(x00l, y0l, q); Zh = BBE_MULRNX16C(x00h, y0h, q);
    BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
    z0l = BBE_PACKVNX40(Zl, q);  z0h = BBE_PACKVNX40(Zh, q);

    Zl = BBE_MULRNX16C(x10l, y0l, q); Zh = BBE_MULRNX16C(x10h, y0h, q);
    BBE_MULANX16C(Zl, x11l, y1l); BBE_MULANX16C(Zh, x11h, y1h);
    z1l = BBE_PACKVNX40(Zl, q);   z1h = BBE_PACKVNX40(Zh, q);

    BBE_SVNX16_I(z0h, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(z0l, pz, 2 * 2 * L);

    BBE_SVNX16_I(z1h, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(z1l, pz, -2 * 2 * L + 2 * 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul2x2s() */

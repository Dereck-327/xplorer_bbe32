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

/* Streaming Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul2x2s ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  int l;

  xb_vecNx16  xx00, xx01, xx11, xx10, zero;
  xb_vecNx16  x00, x01, x11, x10;
  xb_vecNx16  y00, y01, y11, y10;
  xb_vecNx16  z00, z01, z11, z10;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  const int offset01 = (2 * L)*sizeof(int16_t);

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L%BBE_SIMD_WIDTH == 0);

  zero = 0;
  if (L <= 0) return;
  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    xb_vecNx40 r0;

    BBE_LVNX16_XP(xx00, px, 2 * L);
    BBE_LVNX16_XP(xx01, px, 2 * L);
    BBE_LVNX16_XP(xx10, px, 2 * L);
    BBE_LVNX16_XP(xx11, px, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

    x00 = BBE_SELNX16I(zero, xx00, BBE_SELI_INTERLEAVE_1_LO);
    x01 = BBE_SELNX16I(zero, xx01, BBE_SELI_INTERLEAVE_1_LO);
    x10 = BBE_SELNX16I(zero, xx10, BBE_SELI_INTERLEAVE_1_LO);
    x11 = BBE_SELNX16I(zero, xx11, BBE_SELI_INTERLEAVE_1_LO);

    BBE_LVNX16_XP(y00, py, offset01);
    BBE_LVNX16_XP(y01, py, offset01);
    BBE_LVNX16_XP(y10, py, offset01);
    BBE_LVNX16_XP(y11, py, -3 * offset01 + 2 * BBE_SIMD_WIDTH);

    r0 = BBE_MULNX16C(x00, y00); BBE_MULANX16C(r0, x01, y10); z00 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x00, y01); BBE_MULANX16C(r0, y11, x01); z01 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x10, y00); BBE_MULANX16C(r0, x11, y10); z10 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x10, y01); BBE_MULANX16C(r0, x11, y11); z11 = BBE_PACKVNX40(r0, q);

    BBE_SVNX16_XP(z00, pz, offset01);
    BBE_SVNX16_XP(z01, pz, offset01);
    BBE_SVNX16_XP(z10, pz, offset01);
    BBE_SVNX16_XP(z11, pz, -3 * offset01 + 2 * BBE_SIMD_WIDTH);

    x00 = BBE_SELNX16I(zero, xx00, BBE_SELI_INTERLEAVE_1_HI);
    x01 = BBE_SELNX16I(zero, xx01, BBE_SELI_INTERLEAVE_1_HI);
    x10 = BBE_SELNX16I(zero, xx10, BBE_SELI_INTERLEAVE_1_HI);
    x11 = BBE_SELNX16I(zero, xx11, BBE_SELI_INTERLEAVE_1_HI);

    BBE_LVNX16_XP(y00, py, offset01);
    BBE_LVNX16_XP(y01, py, offset01);
    BBE_LVNX16_XP(y10, py, offset01);
    BBE_LVNX16_XP(y11, py, -3 * offset01 + 2 * BBE_SIMD_WIDTH);

    r0 = BBE_MULNX16C(x00, y00); BBE_MULANX16C(r0, x01, y10); z00 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x00, y01); BBE_MULANX16C(r0, y11, x01); z01 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x10, y00); BBE_MULANX16C(r0, x11, y10); z10 = BBE_PACKVNX40(r0, q);
    r0 = BBE_MULNX16C(x10, y01); BBE_MULANX16C(r0, x11, y11); z11 = BBE_PACKVNX40(r0, q);

    BBE_SVNX16_XP(z00, pz, offset01);
    BBE_SVNX16_XP(z01, pz, offset01);
    BBE_SVNX16_XP(z10, pz, offset01);
    BBE_SVNX16_XP(z11, pz, -3 * offset01 + 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatmul2x2s() */

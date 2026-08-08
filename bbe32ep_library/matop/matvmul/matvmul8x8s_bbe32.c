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

/* Streaming Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 16
*/
void matvmul8x8s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  vsaN  q = BBE_MOVVSA32(Q);
  const int L8 = L >> LOG2_BBE_SIMD_WIDTH;
  int i, j, xstride, ystride;
  xb_vecNx40 z0, z2;
  xb_vecNx16 x00, y00;

  const xb_vecNx16 * restrict px01 = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py01;
  xb_vecNx16 * restrict pz00;
  xb_vecNx16 * restrict pz02;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L%BBE_SIMD_WIDTH == 0);

  if (L<=0) return;
  L >>= 1;

  pz00 = (xb_vecNx16*)z;
  pz02 = (xb_vecNx16*)(z + 8 * L);
  py01 = (const xb_vecNx16*)y;

  __Pragma("ymemory( px01 )");
  __Pragma("loop_count min=4 factor=4");
  for (i = j = 0; j<L8 * 4; j++)
  {
    /* i=(i+1);if(i==L8) i=0; */
    i = BBE_ADDMOD16U(i, (L8 << 16) | 1);   
    xstride = ystride = (2 * BBE_SIMD_WIDTH - 4 * 7 * L);
    XT_MOVEQZ(xstride, 2 * BBE_SIMD_WIDTH, i);
    XT_MOVEQZ(ystride, (2 * BBE_SIMD_WIDTH - 4 * 8 * L), i);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    z2 = BBE_MULRNX16(x00, y00, q);

    BBE_LVNX16_XP(x00, px01, 4 * L);
    z0 = BBE_MULRNX16(x00, y00, q);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, 4 * L);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, 4 * L);
    BBE_MULANX16(z0, x00, y00);

    BBE_LVNX16_XP(y00, py01, ystride);
    x00 = BBE_LVNX16_X(px01, 4 * L * 32);
    BBE_MULANX16(z2, x00, y00);
    BBE_LVNX16_XP(x00, px01, xstride);
    BBE_MULANX16(z0, x00, y00);

    x00 = BBE_PACKVNX40(z2, q); BBE_SVNX16_IP(x00, pz02, BBE_SIMD_WIDTH * 2);
    x00 = BBE_PACKVNX40(z0, q); BBE_SVNX16_IP(x00, pz00, BBE_SIMD_WIDTH * 2);
  }
} /* matvmul8x8s() */

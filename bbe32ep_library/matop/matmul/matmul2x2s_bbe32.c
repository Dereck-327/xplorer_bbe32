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

/* Streaming Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void matmul2x2s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int l;

  xb_vecNx16  x00, x01, x11, x10;
  xb_vecNx16  y00, y01, y11, y10;
  xb_vecNx16  z00, z01, z11, z10;

  vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 * restrict px = (xb_vecNx16 *)x;
  xb_vecNx16 * restrict py = (xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);

  if (L <= 0) return;
  NASSERT(L % BBE_SIMD_WIDTH == 0);

  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  __Pragma("loop_count min=1");
  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    xb_vecNx40 r0;

    /* Load input matrix X */
    BBE_LVNX16_XP(x00, px, 2 * L);
    BBE_LVNX16_XP(x01, px, 2 * L);
    BBE_LVNX16_XP(x10, px, 2 * L);
    BBE_LVNX16_XP(x11, px, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVNX16_XP(y00, py, 2 * L);
    BBE_LVNX16_XP(y01, py, 2 * L);
    BBE_LVNX16_XP(y10, py, 2 * L);
    BBE_LVNX16_XP(y11, py, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

    /* Multiply input matrix X and Y */
    r0 = BBE_MULRNX16(x00, y00, q);
    BBE_MULANX16(r0, x01, y10);
    z00 = BBE_PACKVNX40(r0, q);

    r0 = BBE_MULRNX16(x00, y01, q);
    BBE_MULANX16(r0, y11, x01);
    z01 = BBE_PACKVNX40(r0, q);

    r0 = BBE_MULRNX16(x10, y00, q);
    BBE_MULANX16(r0, x11, y10);
    z10 = BBE_PACKVNX40(r0, q);

    r0 = BBE_MULRNX16(x10, y01, q);
    BBE_MULANX16(r0, x11, y11);
    z11 = BBE_PACKVNX40(r0, q);

    /* Save rezults */
    BBE_SVNX16_X(z01, pz, 1 * 2 * L);
    BBE_SVNX16_X(z10, pz, 2 * 2 * L);
    BBE_SVNX16_X(z11, pz, 3 * 2 * L);
    BBE_SVNX16_IP(z00, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* matmul2x2s() */

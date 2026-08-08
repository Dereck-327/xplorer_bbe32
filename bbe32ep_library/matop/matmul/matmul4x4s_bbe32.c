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

/* Streaming Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void matmul4x4s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int i, k, stride, ystride;

  xb_vecNx16 xx00, xx01, xx02, xx03;
  xb_vecNx16 y0_0, y1_0, y2_0, y3_0;
  xb_vecNx16 y0_1, y1_1, y2_1, y3_1;
  xb_vecNx16 y0_2, y1_2, y2_2, y3_2;
  xb_vecNx16 y0_3, y1_3, y2_3, y3_3;

  xb_vecNx16 z_0, z_1, z_2, z_3;

  xb_vecNx40 Z0, Z1, Z2, Z3;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

  if (L <= 0) return;
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L>0 && L%BBE_SIMD_WIDTH == 0);
  __Pragma("ymemory( py )");
  __Pragma("loop_count min=4");
  for (k = i = 0; i<L * 4; i += BBE_SIMD_WIDTH)
  {
    k = BBE_ADDMOD16U(k, (L << 16) | BBE_SIMD_WIDTH);
    stride = BBE_SIMD_WIDTH;
    XT_MOVNEZ(stride, -3 * L + BBE_SIMD_WIDTH, k);
    ystride = -15 * 2 * L + 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(ystride, -15 * 2 * L + 2 * BBE_SIMD_WIDTH - 2 * L, k);

    /* Load input matrix Y */
    BBE_LVNX16_XP(y0_0, py, 2 * L);
    BBE_LVNX16_XP(y0_1, py, 2 * L);
    BBE_LVNX16_XP(y0_2, py, 2 * L);
    BBE_LVNX16_XP(y0_3, py, 2 * L);

    BBE_LVNX16_XP(y1_0, py, 2 * L);
    BBE_LVNX16_XP(y1_1, py, 2 * L);
    BBE_LVNX16_XP(y1_2, py, 2 * L);
    BBE_LVNX16_XP(y1_3, py, 2 * L);

    BBE_LVNX16_XP(y2_0, py, 2 * L);
    BBE_LVNX16_XP(y2_1, py, 2 * L);
    BBE_LVNX16_XP(y2_2, py, 2 * L);
    BBE_LVNX16_XP(y2_3, py, 2 * L);

    BBE_LVNX16_XP(y3_0, py, 2 * L);
    BBE_LVNX16_XP(y3_1, py, 2 * L);
    BBE_LVNX16_XP(y3_2, py, 2 * L);
    BBE_LVNX16_XP(y3_3, py, ystride);

    /* Load input matrix X */
    BBE_LVNX16_XP(xx00, px, 2 * L);
    BBE_LVNX16_XP(xx01, px, 2 * L);
    BBE_LVNX16_XP(xx02, px, 2 * L);
    BBE_LVNX16_XP(xx03, px, 0 * L);

    /* Multipy input matrix X and Y */
    Z0 = BBE_MULRNX16(xx00, y0_0, q); BBE_MULANX16(Z0, xx01, y1_0);
    BBE_MULANX16(Z0, xx02, y2_0); BBE_MULANX16(Z0, xx03, y3_0);
    Z1 = BBE_MULRNX16(xx00, y0_1, q); BBE_MULANX16(Z1, xx01, y1_1);
    BBE_MULANX16(Z1, xx02, y2_1); BBE_MULANX16(Z1, xx03, y3_1);

    Z2 = BBE_MULRNX16(xx00, y0_2, q); BBE_MULANX16(Z2, xx01, y1_2);
    BBE_MULANX16(Z2, xx02, y2_2); BBE_MULANX16(Z2, xx03, y3_2);
    Z3 = BBE_MULRNX16(xx00, y0_3, q); BBE_MULANX16(Z3, xx01, y1_3);
    BBE_MULANX16(Z3, xx02, y2_3); BBE_MULANX16(Z3, xx03, y3_3);

    px = (const xb_vecNx16 *)XT_ADDX2(stride, (int32_t)px);

    /* Pack and save rezult */
    z_0 = BBE_PACKVNX40(Z0, q);
    z_1 = BBE_PACKVNX40(Z1, q);
    z_2 = BBE_PACKVNX40(Z2, q);
    z_3 = BBE_PACKVNX40(Z3, q);
    BBE_SVNX16_XP(z_0, pz, 2 * L);
    BBE_SVNX16_XP(z_1, pz, 2 * L);
    BBE_SVNX16_XP(z_2, pz, 2 * L);
    BBE_SVNX16_XP(z_3, pz, 0 * L);

    pz = (xb_vecNx16 *)XT_ADDX2(stride, (int32_t)pz);
  }
} /* matmul4x4s() */

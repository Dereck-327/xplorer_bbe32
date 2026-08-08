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
    Complex Matrix-Matrix/Matrix-Vector Multiply
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
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
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

/* Streaming Order, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul8x8s ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q )
{
  int i, m, k;
  xb_vecNx16 x0, x1, y0;
  xb_vecNx16 z0, z1;
  xb_vecNx40 Z0, Z1;
  vsaN q = BBE_MOVVSA32(Q);
  const xb_vecNx16 * px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * py = (const xb_vecNx16 *)y;
        xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);
  NASSERT(Q >= 0 && Q <= 16);
  if (L <= 0) return;

  px = (const xb_vecNx16 *)x;
  for (m = 0; m<4; m++)
  {
    int xstride;
    py = (const xb_vecNx16 *)y;
    __Pragma("ymemory( py )");
    __Pragma("loop_count min=8");
    for (k = i = 0; i<L * 8; i += BBE_SIMD_WIDTH / 2)
    {
      k = BBE_ADDMOD16U(k, (L << 16) | (BBE_SIMD_WIDTH / 2));
      xstride = -7 * 4 * L + 2 * BBE_SIMD_WIDTH;
      XT_MOVEQZ(xstride, -8 * 4 * L + 2 * BBE_SIMD_WIDTH, k);
      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      Z0 = BBE_MULNX16C(x0, y0);
      Z1 = BBE_MULNX16C(x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, 4 * L);
      BBE_LVNX16_XP(y0, py, 8 * 4 * L);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      x1 = BBE_LVNX16_X(px, 8 * 4 * L);
      BBE_LVNX16_XP(x0, px, xstride);
      BBE_LVNX16_XP(y0, py, -7 * 8 * 4 * L + 2 * BBE_SIMD_WIDTH);
      BBE_MULANX16C(Z0, x0, y0);
      BBE_MULANX16C(Z1, x1, y0);

      /*Pack and save results*/
      z0 = BBE_PACKVNX40(Z0, q);
      z1 = BBE_PACKVNX40(Z1, q);
      BBE_SVNX16_X(z1, pz, 8 * 4 * L);
      BBE_SVNX16_IP(z0, pz, 2 * BBE_SIMD_WIDTH);
    }
    /* jump over row */
    px = (const xb_vecNx16*)XT_ADDX2(32 * L, (uintptr_t)px);
    pz = (xb_vecNx16*)XT_ADDX8(4 * L, (uintptr_t)pz);
  }
} /* cmatmul8x8s() */

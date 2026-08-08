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

/* Block Order, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     None
*/
void matmul16x16n ( int16_t * restrict z, 
              const int16_t * restrict x, 
              const int16_t * restrict y, 
              int L, int Q )
{
  int k, i;
  vsaN q;
  int stride;

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;
  xb_vecNx16  x0, y0, x_sel;
  xb_vecNx40 acc;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;
  q = BBE_MOVVSA32(Q);
  __Pragma("loop_count min=16, factor=16");
  for (k = i = 0; i<L * 16; i++)
  {
    k = BBE_ADDMOD16U(k, (16 << 16) | 1);
    stride = -15 * 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(stride, 2 * BBE_SIMD_WIDTH, k);

    BBE_LVNX16_IP(x0, px, 2 * BBE_SIMD_WIDTH);
    x_sel = BBE_REPNX16(x0, 0);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(x_sel, y0, q);
    x_sel = BBE_REPNX16(x0, 1);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 2);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 3);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 4);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 5);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 6);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 7);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 8);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 9);  BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 10); BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 11); BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 12); BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 13); BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 14); BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, x_sel, y0);
    x_sel = BBE_REPNX16(x0, 15); BBE_LVNX16_XP(y0, py, stride); BBE_MULANX16(acc, x_sel, y0);
    x0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_IP(x0, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* matmul16x16n() */

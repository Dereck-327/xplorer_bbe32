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
    Matrix Hermitian Product
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
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     None
*/
void cmatherm8x8n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q )
{
  int k, i, repstep, ystride, zstride;
  vsaN q;

  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)x;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)y;
  xb_vecNx16  x0, y0, y1, y2, y3, y4, y5, y6, y7, x_sel;
  xb_vecNx40 acc;
  vselN rep0, rep1;

  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  if (L <= 0) return;

  q = BBE_MOVVSA32(Q);
  x0 = BBE_MOVVA16C(1 << 16); /* pair 1,0 */
  rep0 = BBE_MOVVSELNX16(x0, 0);
  x0 = BBE_MOVVA16C(0x30002); /* pair 3,2 */
  rep1 = BBE_MOVVSELNX16(x0, 0);

  __Pragma("loop_count min=4, factor=4");
  for (i = k = 0; k<L * 4; k++)
  {
    i = BBE_ADDMOD16U(i, (L << 16) | 1);    /*m=(m+1)&3;*/
    ystride = 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(ystride, -64 * 4 * L + 2 * BBE_SIMD_WIDTH, i);
    repstep = 0;
    XT_MOVEQZ(repstep, 4, i);
    zstride = 7 * 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(zstride, -64 * 4 * L + 9 * 2 * BBE_SIMD_WIDTH, i);

    BBE_LVNX16_IP(y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y2, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y3, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y4, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y5, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y6, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y7, py, ystride);

    x_sel = BBE_SHFLNX16(y0, rep0);  acc = BBE_MULRNX16J(y0, x_sel, q);
    x_sel = BBE_SHFLNX16(y1, rep0);  BBE_MULANX16J(acc, y1, x_sel);
    x_sel = BBE_SHFLNX16(y2, rep0);  BBE_MULANX16J(acc, y2, x_sel);
    x_sel = BBE_SHFLNX16(y3, rep0);  BBE_MULANX16J(acc, y3, x_sel);
    x_sel = BBE_SHFLNX16(y4, rep0);  BBE_MULANX16J(acc, y4, x_sel);
    x_sel = BBE_SHFLNX16(y5, rep0);  BBE_MULANX16J(acc, y5, x_sel);
    x_sel = BBE_SHFLNX16(y6, rep0);  BBE_MULANX16J(acc, y6, x_sel);
    BBE_SHFLUNX16(x_sel, y7, rep0, repstep);  BBE_MULANX16J(acc, y7, x_sel);
    x0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_IP(x0, pz, 2 * BBE_SIMD_WIDTH);

    x_sel = BBE_SHFLNX16(y0, rep1);  acc = BBE_MULRNX16J(y0, x_sel, q);
    x_sel = BBE_SHFLNX16(y1, rep1);  BBE_MULANX16J(acc, y1, x_sel);
    x_sel = BBE_SHFLNX16(y2, rep1);  BBE_MULANX16J(acc, y2, x_sel);
    x_sel = BBE_SHFLNX16(y3, rep1);  BBE_MULANX16J(acc, y3, x_sel);
    x_sel = BBE_SHFLNX16(y4, rep1);  BBE_MULANX16J(acc, y4, x_sel);
    x_sel = BBE_SHFLNX16(y5, rep1);  BBE_MULANX16J(acc, y5, x_sel);
    x_sel = BBE_SHFLNX16(y6, rep1);  BBE_MULANX16J(acc, y6, x_sel);
    BBE_SHFLUNX16(x_sel, y7, rep1, repstep); BBE_MULANX16J(acc, y7, x_sel);
    x0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(x0, pz, zstride);
  }
} /* cmatherm8x8n() */

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
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatherm2x2n,(complex_fract16* restrict y, 
               const complex_fract16 * restrict x, 
               int L, int Q))
#else
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

/* Block Order, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be even
*/
void cmatherm2x2n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q )
{
  int l;
  vsaN  q = BBE_MOVVSA32(Q);
  static const int16_t ALIGN(32) seli[BBE_SIMD_WIDTH] = { 0, 1, 0, 1, 2, 3, 2, 3, 8, 9, 8, 9, 10, 11, 10, 11 };
  const xb_vecNx16 *px = (const xb_vecNx16 *)x;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)y;
  xb_vecNx16 X0, X1, Y0, Y1, Z;
  xb_vecNx40 r;
  vselN sel0, sel1;
  /* check restrictions */
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  NASSERT(L % 2 == 0);
  if (L <= 0) return;
  Z = BBE_LVNX16_I((const xb_vecNx16*)seli, 0);
  sel0 = BBE_MOVVSELNX16(Z, 0);
  sel1 = sel0;
  BBE_SELUNX16(Z, Z, Z, sel1, 4);
  __Pragma("loop_count min=1");
  for (l = 0; l<L; l += 2)
  {
    BBE_LVNX16_IP(Z, px, 2 * BBE_SIMD_WIDTH);
    X0 = BBE_SHFLNX16(Z, sel0);
    Y0 = BBE_SHFLNX16I(Z, BBE_SHFLI_MMC2X2X2X2_M2_STEP_1);
    X1 = BBE_SHFLNX16(Z, sel1);
    Y1 = BBE_SHFLNX16I(Z, BBE_SHFLI_MMC2X2X2X2_M2_STEP_2);
    r = BBE_MULRNX16J(Y0, X0, q);
    BBE_MULANX16J(r, Y1, X1);
    /* Pack and save results */
    Z = BBE_PACKVNX40(r, q);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* cmatherm2x2n() */
#endif

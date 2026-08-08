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

#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatmul2x2n, (complex_fract16 * restrict z,
  const complex_fract16 * restrict x,
  const complex_fract16 * restrict y,
  int L, int Q))
#else

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

/* Block Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void cmatmul2x2n(complex_fract16 * restrict z,
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y,
             int L, int Q )
{
  int l;
  vsaN  q = BBE_MOVVSA32(Q);
  const xb_vecNx16 *px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;
  xb_vecNx16 X, Y, Z;
  xb_vecNx16 sel_x0, sel_y0, sel_x1, sel_y1;
  xb_vecNx40 r;

  /* check restrictions */
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2*BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L % 2 == 0);

  if (L<=0) return;
  if (Q == 15)
  {
    __Pragma("loop_count min=1");
    for (l = 0; l<L; l += 2)
    {
      /* Load input matrices X and Y */
      BBE_LVNX16_XP(X, px, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(Y, py, 2 * BBE_SIMD_WIDTH);
      /* Select input matrix for multiple */
      sel_x0 = BBE_SHFLNX16I(X, BBE_SHFLI_MMC2X2X2X2_M1_STEP_1);
      sel_y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M2_STEP_1);
      /* Multiple input matrix X and Y */
      r = BBE_MULNX16C(sel_x0, sel_y0);
      /* Select input matrix for multiple */
      sel_x1 = BBE_SHFLNX16I(X, BBE_SHFLI_MMC2X2X2X2_M1_STEP_2);
      sel_y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M2_STEP_2);
      /* Multiple input matrix X and Y */
      BBE_MULANX16C(r, sel_x1, sel_y1);
      /* Pack and save results */
      Z = BBE_PACKQNX40(r);
      BBE_SVNX16_XP(Z, pz, 2 * BBE_SIMD_WIDTH);
    }
  }
  else
  {
    __Pragma("loop_count min=1");
    for (l = 0; l<L; l += 2)
    {
      /* load input matrix X and Y */
      BBE_LVNX16_IP(X, px, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);
      /* Select input matrix for multiple */
      sel_x0 = BBE_SHFLNX16I(X, BBE_SHFLI_MMC2X2X2X2_M1_STEP_1);
      sel_y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M2_STEP_1);
      /* Multiple input matrix X and Y */
      r = BBE_MULRNX16C(sel_x0, sel_y0, q);
      /* Select input matrix for multiple */
      sel_x1 = BBE_SHFLNX16I(X, BBE_SHFLI_MMC2X2X2X2_M1_STEP_2);
      sel_y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_MMC2X2X2X2_M2_STEP_2);
      /* Multiple input matrix X and Y */
      BBE_MULANX16C(r, sel_x1, sel_y1);
      /* Pack and save results */
      Z = BBE_PACKVNX40(r, q);
      BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
    }
  }
} /* cmatmul2x2n() */
#endif

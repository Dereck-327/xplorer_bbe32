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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, matvmul4x4n,(int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q))
#else
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

/* Block Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 4
*/
void matvmul4x4n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  int l;

  xb_vecNx16  y_, x0, x1, x2, x3, t_;
  xb_vecNx16  y0123, y4567;
  xb_vecNx40 acc0, acc1;

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  vsaN  q = BBE_MOVVSA32(Q);

  static const int16_t ALIGN(32) sel[16] = { 0, 4, 8, 12, 2, 6, 10, 14, 16, 20, 24, 28, 18, 22, 26, 30 };

  const xb_vecNx16 *SEL = (const xb_vecNx16 *)sel;

  xb_vecNx16 Sel;

  vselN SelN;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % (BBE_SIMD_WIDTH / 4) == 0);

  if (L<=0) return;
  Sel = BBE_LVNX16_I(SEL, 0);
  SelN = BBE_MOVVSV(Sel, 0);

  __Pragma("ymemory( px )");
  for (l = 0; l<L; l += 4)
  {
    BBE_LVNX16_IP(x0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x2, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x3, px, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(y_, py, 2 * BBE_SIMD_WIDTH);

    y_ = BBE_CONJSNX16C(y_);
    y0123 = BBE_SHFLNX16I(y_, BBE_SHFLI_REP_0X4);
    y4567 = BBE_SHFLNX16I(y_, BBE_SHFLI_REP_1X4);

    acc0 = BBE_MULRNX16PC_0(x0, y0123, q);
    BBE_MULANX16PC_1(acc0, x1, y4567);

    y0123 = BBE_SHFLNX16I(y_, BBE_SHFLI_REP_2X4);
    y4567 = BBE_SHFLNX16I(y_, BBE_SHFLI_REP_3X4);

    acc1 = BBE_MULRNX16PC_0(x2, y0123, q);
    BBE_MULANX16PC_1(acc1, x3, y4567);

    x0 = BBE_PACKVNX40(acc0, q);
    x1 = BBE_PACKVNX40(acc1, q);

    t_ = BBE_SELNX16(x1, x0, SelN);

    BBE_SVNX16_IP(t_, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* matvmul4x4n() */
#endif

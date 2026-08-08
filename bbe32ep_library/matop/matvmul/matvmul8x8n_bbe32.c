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
DISCARD_FUN(void, matvmul8x8n,(int16_t * restrict z, 
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

/* Block Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
      L must be even
*/
void matvmul8x8n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  int l;

  vsaN  q = BBE_MOVVSA32(Q);

  xb_vecNx16 Y, Y0, Y1, Y2, Y3, X0, X1, X2, X3, X4, X5, X6, X7, X_0, X_1, X_2, X_3, X_4, X_5, X_6, X_7, Z, Z0, Z1;

  xb_vecNx40 Acc0, Acc1;

  const xb_vecNx16 * restrict X_ = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict Y_ = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict Z_ = (xb_vecNx16 *)z;

  static const int16_t ALIGN(32) sel[16] = { 0, 4, 8, 12, 2, 6, 10, 14, 16, 20, 24, 28, 18, 22, 26, 30 };

  const xb_vecNx16 *restrict SEL = (const xb_vecNx16 *)sel;

  xb_vecNx16 Sel;

  vselN SelN;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % (BBE_SIMD_WIDTH / 8) == 0);

  if (L<=0) return;
  Sel = BBE_LVNX16_I(SEL, 0);
  SelN = BBE_MOVVSV(Sel, 0);

  __Pragma("ymemory( X_ )");
  for (l = 0; l<L; l += 2)
  {
    BBE_LVNX16_IP(Y, Y_, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(X0, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, X_, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(X4, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X5, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X6, X_, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X7, X_, 2 * BBE_SIMD_WIDTH);

    Y = BBE_CONJSNX16C(Y);

    Y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_0X4);
    Y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_1X4);
    Y2 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_2X4);
    Y3 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_3X4);

    //Multipy first and second matrix
    X_0 = BBE_SELNX16I(X1, X0, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    X_4 = BBE_SELNX16I(X5, X4, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    Acc0 = BBE_MULRNX16PC_0(X_0, Y0, q);
    Acc1 = BBE_MULRNX16PC_0(X_4, Y2, q);

    X_1 = BBE_SELNX16I(X1, X0, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    X_5 = BBE_SELNX16I(X5, X4, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    BBE_MULANX16PC_0(Acc0, X_1, Y1);
    BBE_MULANX16PC_0(Acc1, X_5, Y3);

    X_2 = BBE_SELNX16I(X3, X2, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    X_6 = BBE_SELNX16I(X7, X6, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    BBE_MULANX16PC_1(Acc0, X_2, Y0);
    BBE_MULANX16PC_1(Acc1, X_6, Y2);

    X_3 = BBE_SELNX16I(X3, X2, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    X_7 = BBE_SELNX16I(X7, X6, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    BBE_MULANX16PC_1(Acc0, X_3, Y1);
    BBE_MULANX16PC_1(Acc1, X_7, Y3);

    //Pack and save rezult
    Z0 = BBE_PACKVNX40(Acc0, q); Z1 = BBE_PACKVNX40(Acc1, q);
    Z = BBE_SELNX16(Z1, Z0, SelN);
    BBE_SVNX16_IP(Z, Z_, 2 * BBE_SIMD_WIDTH);
  }
} /* matvmul8x8n() */
#endif

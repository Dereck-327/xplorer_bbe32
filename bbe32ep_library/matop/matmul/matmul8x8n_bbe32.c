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
DISCARD_FUN(void, matmul8x8n,(int16_t * restrict z, 
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

/* Block Order, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void matmul8x8n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int l;

  //static const int16_t ALIGN(32) sel0[16] = {0, 8, 16, 24, 1, 9, 17, 25,  2, 10, 18, 26,  3, 11, 19, 27};
  //static const int16_t ALIGN(32) sel1[16] = { 0, 4, 8, 12, 2, 6, 10, 14, 16, 20, 24, 28, 18, 22, 26, 30 };

  static const int16_t ALIGN(32) sel0[16] = { 0, 8, 16, 24, 1, 9, 17, 25, 2, 10, 18, 26, 3, 11, 19, 27 };
  static const int16_t ALIGN(32) sel1[16] = { 1, 5, 9, 13, 3, 7, 11, 15, 17, 21, 25, 29, 19, 23, 27, 31 };

  vsaN  q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          X_ = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          Y_ = (const xb_vecNx16 *)y;
        xb_vecNx16 * restrict Z_ = (      xb_vecNx16 *)z;

  xb_vecNx16 X, Y0, Y1, Y2, Y3, X0, X1, X2, X3, Y_0, Y_1, Y_2, Y_3, Z0, Z1, Z;

  const xb_vecNx16 *SEL0 = (const xb_vecNx16 *)sel0;
  const xb_vecNx16 *SEL1 = (const xb_vecNx16 *)sel1;

  xb_vecNx16 Sel0, Sel1;

  xb_vecNx40 Acc;

  vselN SelN0, SelN1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;

  Sel0  = BBE_LVNX16_I(SEL0, 0);
  SelN0 = BBE_MOVVSV(Sel0, 0);

  Sel1  = BBE_LVNX16_I(SEL1, 0);
  SelN1 = BBE_MOVVSV(Sel1, 0);
    
  __Pragma( "ymemory( Y_ )" );
  __Pragma( "ymemory( X_ )" );
  __Pragma( "loop_count min=2" );
  for (l=0; l<2*L; l++)
  {
      BBE_LVNX16_IP(Y_0, Y_, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Y_1, Y_, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Y_2, Y_, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Y_3, Y_, 2*BBE_SIMD_WIDTH);

      Y0 = BBE_SELNX16(Y_1, Y_0, SelN0);
      BBE_SELUNX16(Y2, Y_3, Y_2, SelN0, 4);
      Y1 = BBE_SELNX16(Y_1, Y_0, SelN0);
      BBE_SELUNX16(Y3, Y_3, Y_2, SelN0, 124);

      BBE_LVNX16_IP(X, X_, 2*BBE_SIMD_WIDTH);
      X = BBE_SHFLNX16I(X, BBE_SHFLI_SWAP_1);

      X0 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_0X4);
      X1 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_1X4);
      X2 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_2X4);
      X3 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_3X4);

      Acc = BBE_MULRNX16PC_0(Y0, X0, q);
      BBE_MULANX16PC_0(Acc, Y2, X1);
      BBE_MULANX16PC_1(Acc, Y1, X0);
      BBE_MULANX16PC_1(Acc, Y3, X1);

      Z0 = BBE_PACKVNX40(Acc, q);
        
      Acc = BBE_MULRNX16PC_0(Y0, X2, q);
      BBE_MULANX16PC_0(Acc, Y2, X3);
      BBE_MULANX16PC_1(Acc, Y1, X2);
      BBE_MULANX16PC_1(Acc, Y3, X3);

      Z1 = BBE_PACKVNX40(Acc, q);

      BBE_LVNX16_IP(X, X_, 2*BBE_SIMD_WIDTH);
      X = BBE_SHFLNX16I(X, BBE_SHFLI_SWAP_1);

      Z = BBE_SELNX16(Z1, Z0, SelN1);
      BBE_SVNX16_IP(Z, Z_, 2*BBE_SIMD_WIDTH);

      X0 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_0X4);
      X1 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_1X4);
      X2 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_2X4);
      X3 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_3X4);

      Acc = BBE_MULRNX16PC_0(Y0, X0,q);
      BBE_MULANX16PC_0(Acc, Y2, X1);
      BBE_MULANX16PC_1(Acc, Y1, X0);
      BBE_MULANX16PC_1(Acc, Y3, X1);

      Z0 = BBE_PACKVNX40(Acc, q);
        
      Acc = BBE_MULRNX16PC_0(Y0, X2, q);
      BBE_MULANX16PC_0(Acc, Y2, X3);
      BBE_MULANX16PC_1(Acc, Y1, X2);
      BBE_MULANX16PC_1(Acc, Y3, X3);

      Z1 = BBE_PACKVNX40(Acc, q);

      Z = BBE_SELNX16(Z1, Z0, SelN1);
      BBE_SVNX16_IP(Z, Z_, 2*BBE_SIMD_WIDTH); 

      Y_ = (const xb_vecNx16 *)XT_ADDX4((l&1)*2*BBE_SIMD_WIDTH, (int32_t)Y_);
      Y_ = (const xb_vecNx16 *)XT_ADDX4(-2*BBE_SIMD_WIDTH, (int32_t)Y_);
  }
} /* matmul8x8n() */
#endif

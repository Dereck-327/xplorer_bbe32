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
DISCARD_FUN(void, matmul4x4n,(int16_t * restrict z, 
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

/* Block Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void matmul4x4n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
  int l;

  static const int16_t ALIGN(32) sel0[16] = {0, 4, 8, 12, 1, 5,  9, 13,  2,  6, 10, 14,  3,  7, 11, 15};
  static const int16_t ALIGN(32) sel1[16] = {0, 4, 8, 12, 2, 6, 10, 14, 16, 20, 24, 28, 18, 22, 26, 30};

	vsaN  q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *SEL0 = (const xb_vecNx16 *)sel0;
  const xb_vecNx16 *SEL1 = (const xb_vecNx16 *)sel1;

  xb_vecNx16 Sel0, Sel1;

  vselN SelN0, SelN1;

  const xb_vecNx16 *          X_ = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          Y_ = (const xb_vecNx16 *)y;
        xb_vecNx16 * restrict Z_ = (      xb_vecNx16 *)z;

  xb_vecNx16 X, Y, X0, X1, X2, X3, Z, Z0, Z1;

  xb_vecNx40 Acc0, Acc1;

  NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);;
  NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);;
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);;
  NASSERT(Q>=0 && Q<=16);
  
  if (L<=0) return;

  Sel0  = BBE_LVNX16_I(SEL0, 0);
  SelN0 = BBE_MOVVSV(Sel0, 0);

  Sel1  = BBE_LVNX16_I(SEL1, 0);
  SelN1 = BBE_MOVVSV(Sel1, 0);
    
  __Pragma( "ymemory( Y_ )" );
  __Pragma( "ymemory( X_ )" );
  for (l=0; l<L; l++)
  {
      BBE_LVNX16_IP(X, X_, 2*BBE_SIMD_WIDTH);
      X = BBE_CONJSNX16C(X);

      X0 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_0X4);
      X1 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_1X4);
      X2 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_2X4);
      X3 = BBE_SHFLNX16I(X, BBE_SHFLI_REP_3X4);
 
      BBE_LVNX16_IP(Y, Y_, 2*BBE_SIMD_WIDTH);
      Y = BBE_SHFLNX16(Y, SelN0);

      Acc0 = BBE_MULRNX16PC_0(Y, X0, q);
      BBE_MULANX16PC_1(Acc0, Y, X1);
      Acc1 = BBE_MULRNX16PC_0(Y, X2, q);
      BBE_MULANX16PC_1(Acc1, Y, X3);

      Z0 = BBE_PACKVNX40(Acc0, q); 
      Z1 = BBE_PACKVNX40(Acc1, q);
      Z = BBE_SELNX16(Z1, Z0, SelN1);

      BBE_SVNX16_IP(Z, Z_, 2*BBE_SIMD_WIDTH);
  }
} /* matmul4x4n() */
#endif

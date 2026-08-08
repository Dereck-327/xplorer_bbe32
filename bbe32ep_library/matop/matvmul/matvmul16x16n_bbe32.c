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
DISCARD_FUN(void, matvmul16x16n,(int16_t * restrict z, 
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

/* Block Order, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     None
*/
void matvmul16x16n ( int16_t * restrict z, 
               const int16_t * restrict x, 
               const int16_t * restrict y, 
               int L, int Q )
{
  int l;

  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X_0, X_1, X_2, X_3, Y, Y0, Y1, Y2, Y3, Z, Z0, Z1;

  xb_vecNx40 Acc;

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
  if (L<=0) return;

  Sel = BBE_LVNX16_I(SEL, 0);
  SelN = BBE_MOVVSV(Sel, 0);

  BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

  BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X2, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X3, px, 2 * BBE_SIMD_WIDTH);

  __Pragma("ymemory( px )");
  for (l = 0; l<L - 1; l++)
  {
    BBE_LVNX16_IP(X4, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X5, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X6, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X7, px, 2 * BBE_SIMD_WIDTH);

    Y = BBE_CONJSNX16C(Y);

    Y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_0X4);
    Y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_1X4);
    Y2 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_2X4);
    Y3 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_3X4);

    /* Multiply 0..3 strings */
    BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_1, X_0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_3, X_2, X3, X2, BBE_DSELI_INTERLEAVE_4);

    Acc = BBE_MULRNX16PC_0(X_0, Y0, q);
    BBE_MULANX16PC_0(Acc, X_1, Y1);
    BBE_MULANX16PC_0(Acc, X_2, Y2);
    BBE_MULANX16PC_0(Acc, X_3, Y3);

    /* Multiply 4..7 strings */
    BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_1, X_0, X5, X4, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_3, X_2, X7, X6, BBE_DSELI_INTERLEAVE_4);

    BBE_MULANX16PC_1(Acc, X_0, Y0);
    BBE_MULANX16PC_1(Acc, X_1, Y1);
    BBE_MULANX16PC_1(Acc, X_2, Y2);
    BBE_MULANX16PC_1(Acc, X_3, Y3);

    /* Pack 0..7 strings */
    Z0 = BBE_PACKVNX40(Acc, q);

    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, px, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(X4, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X5, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X6, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X7, px, 2 * BBE_SIMD_WIDTH);

    /* Multiply 8..11 strings */
    BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_1, X_0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_3, X_2, X3, X2, BBE_DSELI_INTERLEAVE_4);

    Acc = BBE_MULRNX16PC_0(X_0, Y0, q);
    BBE_MULANX16PC_0(Acc, X_1, Y1);
    BBE_MULANX16PC_0(Acc, X_2, Y2);
    BBE_MULANX16PC_0(Acc, X_3, Y3);

    /* Multiply 12..15 strings */
    BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_1, X_0, X5, X4, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X_3, X_2, X7, X6, BBE_DSELI_INTERLEAVE_4);

    BBE_MULANX16PC_1(Acc, X_0, Y0);
    BBE_MULANX16PC_1(Acc, X_1, Y1);
    BBE_MULANX16PC_1(Acc, X_2, Y2);
    BBE_MULANX16PC_1(Acc, X_3, Y3);

    /* Pack 8..15 strings */
    Z1 = BBE_PACKVNX40(Acc, q);

    BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, px, 2 * BBE_SIMD_WIDTH);

    /* Save results */
    Z = BBE_SELNX16(Z1, Z0, SelN);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
  }

  BBE_LVNX16_IP(X4, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X5, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X6, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X7, px, 2 * BBE_SIMD_WIDTH);

  Y = BBE_CONJSNX16C(Y);

  Y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_0X4);
  Y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_1X4);
  Y2 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_2X4);
  Y3 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_3X4);

  /* Multiply 0..3 strings */
  BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_1, X_0, X1, X0, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_3, X_2, X3, X2, BBE_DSELI_INTERLEAVE_4);

  Acc = BBE_MULRNX16PC_0(X_0, Y0, q);
  BBE_MULANX16PC_0(Acc, X_1, Y1);
  BBE_MULANX16PC_0(Acc, X_2, Y2);
  BBE_MULANX16PC_0(Acc, X_3, Y3);

  /* Multiply 4..7 strings */
  BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_1, X_0, X5, X4, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_3, X_2, X7, X6, BBE_DSELI_INTERLEAVE_4);

  BBE_MULANX16PC_1(Acc, X_0, Y0);
  BBE_MULANX16PC_1(Acc, X_1, Y1);
  BBE_MULANX16PC_1(Acc, X_2, Y2);
  BBE_MULANX16PC_1(Acc, X_3, Y3);

  /* Pack 0..7 strings */
  Z0 = BBE_PACKVNX40(Acc, q);

  BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X2, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X3, px, 2 * BBE_SIMD_WIDTH);

  BBE_LVNX16_IP(X4, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X5, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X6, px, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(X7, px, 2 * BBE_SIMD_WIDTH);

  /* Multiply 8..11 strings */
  BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_1, X_0, X1, X0, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_3, X_2, X3, X2, BBE_DSELI_INTERLEAVE_4);

  Acc = BBE_MULRNX16PC_0(X_0, Y0, q);
  BBE_MULANX16PC_0(Acc, X_1, Y1);
  BBE_MULANX16PC_0(Acc, X_2, Y2);
  BBE_MULANX16PC_0(Acc, X_3, Y3);

  /* Multiply 12..15 strings */
  BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_1, X_0, X5, X4, BBE_DSELI_INTERLEAVE_4);
  BBE_DSELNX16I(X_3, X_2, X7, X6, BBE_DSELI_INTERLEAVE_4);

  BBE_MULANX16PC_1(Acc, X_0, Y0);
  BBE_MULANX16PC_1(Acc, X_1, Y1);
  BBE_MULANX16PC_1(Acc, X_2, Y2);
  BBE_MULANX16PC_1(Acc, X_3, Y3);

  /* Pack 8..15 strings */
  Z1 = BBE_PACKVNX40(Acc, q);

  /* Save results */
  Z = BBE_SELNX16(Z1, Z0, SelN);
  BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
} /* matvmul16x16n() */
#endif

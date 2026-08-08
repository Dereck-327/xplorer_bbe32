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

/* Streaming Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void matvmul4x4s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  int l;

  const xb_vecNx16 *          X_ = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          Y_ = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict Z_ = (xb_vecNx16 *)z;

  xb_vecNx16 X00, X01, X02, X03;
  xb_vecNx16 Y0, Y1, Y2, Y3;
  xb_vecNx16 z0, z1;

  xb_vecNx40 Z0, Z1;

  vsaN q = BBE_MOVVSA32(Q);

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

  NASSERT(Q >= 0 && Q <= 16);

  NASSERT(L%BBE_SIMD_WIDTH == 0);
  if (L <= 0) return;
  __Pragma("ymemory( X_ )");
  __Pragma("loop_count min=1");
  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    /* Load input matrix X */
    BBE_LVNX16_XP(X00, X_, 2 * L);
    BBE_LVNX16_XP(X01, X_, 2 * L);
    BBE_LVNX16_XP(X02, X_, 2 * L);
    BBE_LVNX16_XP(X03, X_, 2 * L);

    /* Load input matrix Y */
    BBE_LVNX16_XP(Y0, Y_, 2 * L); 
    BBE_LVNX16_XP(Y1, Y_, 2 * L);
    BBE_LVNX16_XP(Y2, Y_, 2 * L);
    BBE_LVNX16_XP(Y3, Y_, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

    /* Multiply input matrix X and Y */
    Z0 = BBE_MULRNX16(X00, Y0, q);
    BBE_MULANX16(Z0, X01, Y1);
    BBE_MULANX16(Z0, X02, Y2);
    BBE_MULANX16(Z0, X03, Y3);

    /* Load input matrix X */
    BBE_LVNX16_XP(X00, X_, 2 * L);
    BBE_LVNX16_XP(X01, X_, 2 * L);
    BBE_LVNX16_XP(X02, X_, 2 * L);
    BBE_LVNX16_XP(X03, X_, 2 * L);

    /* Multiply input matrix X and Y */
    Z1 = BBE_MULRNX16(X00, Y0, q);
    BBE_MULANX16(Z1, X01, Y1);
    BBE_MULANX16(Z1, X02, Y2);
    BBE_MULANX16(Z1, X03, Y3);

    /* Pack and save rezult */
    z0 = BBE_PACKVNX40(Z0, q);  z1 = BBE_PACKVNX40(Z1, q);
    BBE_SVNX16_XP(z0, Z_, 2 * L); BBE_SVNX16_XP(z1, Z_, 2 * L);

    /* Load input matrix X */
    BBE_LVNX16_XP(X00, X_, 2 * L);
    BBE_LVNX16_XP(X01, X_, 2 * L);
    BBE_LVNX16_XP(X02, X_, 2 * L);
    BBE_LVNX16_XP(X03, X_, 2 * L);

    /* Multiply input matrix X and Y */
    Z0 = BBE_MULRNX16(X00, Y0, q);
    BBE_MULANX16(Z0, X01, Y1);
    BBE_MULANX16(Z0, X02, Y2);
    BBE_MULANX16(Z0, X03, Y3);

    /* Load input matrix X */
    BBE_LVNX16_XP(X00, X_, 2 * L);
    BBE_LVNX16_XP(X01, X_, 2 * L);
    BBE_LVNX16_XP(X02, X_, 2 * L);
    BBE_LVNX16_XP(X03, X_, -15 * 2 * L + 2 * BBE_SIMD_WIDTH);

    /* Multiply input matrix X and Y */
    Z1 = BBE_MULRNX16(X00, Y0, q);
    BBE_MULANX16(Z1, X01, Y1);
    BBE_MULANX16(Z1, X02, Y2);
    BBE_MULANX16(Z1, X03, Y3);

    /* Pack and save rezult */
    z0 = BBE_PACKVNX40(Z0, q);  z1 = BBE_PACKVNX40(Z1, q);
    BBE_SVNX16_XP(z0, Z_, 2 * L); BBE_SVNX16_XP(z1, Z_, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);
  }
} /* matvmul4x4s() */

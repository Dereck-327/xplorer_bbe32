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
    Real Matrix by Complex Matrix/Vector Multiply
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
DISCARD_FUN(void, rcmatvmul8x8n,(complex_fract16 * restrict z, 
            const int16_t * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
#else
/*-------------------------------------------------------------------------
Real Matrix by Complex Matrix/Vector Multiply 

Description: These functions perform pairwise multiplication of left-hand
real matrices by right-hand complex matrices or vectors. Both the block order
and streaming order are allowed for input/output matrix sequences.

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand real matrices
y[L*Sy]     Sequence of right-hand complex matrices or vectors
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

/* Block Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     None
*/
void rcmatvmul8x8n ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 X01, X23, X45, X67, X0246_0, X0246_1, X1357_0, X1357_1;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, Y, Y0, Y1, Y2, Y3, Z;

  xb_vecNx40 Acc;

  vsaN  q = BBE_MOVVSA32(Q);

  const xb_vecNx16 zero = 0;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;
  __Pragma("ymemory( px )");
  __Pragma("ymemory( py )");
  __Pragma("loop_count min=1");
  for (l = 0; l<L; l++)
  {
    /* Load input vector Y */
    BBE_LVNX16_IP(Y, py, 2 * BBE_SIMD_WIDTH);

    Y0 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_0X4);
    Y1 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_1X4);
    Y2 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_2X4);
    Y3 = BBE_SHFLNX16I(Y, BBE_SHFLI_REP_3X4);

    /* load input matrix X */
    BBE_LVNX16_IP(X01, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X23, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X45, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X67, px, 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(X45, X01, X45, X01, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X67, X23, X67, X23, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X0246_1, X0246_0, X23, X01, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X1357_1, X1357_0, X67, X45, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(X1, X0, zero, X0246_0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, zero, X0246_1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X5, X4, zero, X1357_0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, zero, X1357_1, BBE_DSELI_INTERLEAVE_1);

    Acc = BBE_MULRNX16PC_0(Y0, X0, q);
    BBE_MULANX16PC_0(Acc, Y1, X1);
    BBE_MULANX16PC_0(Acc, Y2, X2);
    BBE_MULANX16PC_0(Acc, Y3, X3);
    BBE_MULANX16PC_1(Acc, Y0, X4);
    BBE_MULANX16PC_1(Acc, Y1, X5);
    BBE_MULANX16PC_1(Acc, Y2, X6);
    BBE_MULANX16PC_1(Acc, Y3, X7);

    /* Pack and save results */
    Z = BBE_PACKVNX40(Acc, q);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul8x8n() */
#endif

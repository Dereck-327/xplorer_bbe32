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
DISCARD_FUN(void, matvmul3x3n,(int16_t * restrict z, 
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

/* Block Order, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void matvmul3x3n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q )
{
  static
    const int16_t ALIGN(32) aSel[BBE_SIMD_WIDTH] = { 0x0000, 0x0104, 0x0208, 0x0f0c, 0x0302, 0x0406, 0x050a, 0x0f0e, 0x0610, 0x0714, 0x0818, 0x0f1c, 0x0f12, 0x0f16, 0x0f1a, 0x0f1e };
  int l;

  xb_vecNx16  vY, vX0, vX1, vX2, vX3, vZ;
  xb_vecNx16  vY0123, vY4567;
  xb_vecNx40 wvAcc0, wvAcc1;

  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  const vsaN  q = BBE_MOVVSA32(Q);

  xb_vecNx16 vTmp;

  vselN vsSelN, vsConv;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L % (BBE_SIMD_WIDTH / 4) == 0);

  if (L<=0) return;
  /* load selections */
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel, 0);
  vsSelN = BBE_MOVVSV(vTmp, 0);
  vsConv = BBE_MOVVSV(vTmp, 8);

  const vboolN bLoad = BBE_LTRNI(15);

  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += 4)
  {
    /* load 4 matrices X */
    BBE_LVNX16T_IP(vX0, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);
    BBE_LVNX16T_IP(vX1, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);
    BBE_LVNX16T_IP(vX2, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);
    BBE_LVNX16T_IP(vX3, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);

    /* Convert 3x3 to 4x4 */
    vX0 = BBE_SHFLNX16(vX0, vsConv);
    vX1 = BBE_SHFLNX16(vX1, vsConv);
    vX2 = BBE_SHFLNX16(vX2, vsConv);
    vX3 = BBE_SHFLNX16(vX3, vsConv);

    /* load 4 vectors (into 1 register) */
    BBE_LVNX16_IP(vY, pYrd, 2 * BBE_SIMD_WIDTH);

    vY = BBE_CONJSNX16C(vY);
    vY0123 = BBE_SHFLNX16I(vY, BBE_SHFLI_REP_0X4);
    vY4567 = BBE_SHFLNX16I(vY, BBE_SHFLI_REP_1X4);

    wvAcc0 = BBE_MULRNX16PC_0(vX0, vY0123, q);
    BBE_MULANX16PC_1(wvAcc0, vX1, vY4567);

    vY0123 = BBE_SHFLNX16I(vY, BBE_SHFLI_REP_2X4);
    vY4567 = BBE_SHFLNX16I(vY, BBE_SHFLI_REP_3X4);

    wvAcc1 = BBE_MULRNX16PC_0(vX2, vY0123, q);
    BBE_MULANX16PC_1(wvAcc1, vX3, vY4567);

    /* pack and store res */
    vX0 = BBE_PACKVNX40(wvAcc0, q);
    vX1 = BBE_PACKVNX40(wvAcc1, q);

    vZ = BBE_SELNX16(vX1, vX0, vsSelN);

    BBE_SVNX16_IP(vZ, pZwr, 2 * BBE_SIMD_WIDTH);
  }
} /* matvmul3x3n() */
#endif

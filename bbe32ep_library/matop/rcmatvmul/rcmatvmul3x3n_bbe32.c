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
DISCARD_FUN(void, rcmatvmul3x3n,(complex_fract16 * restrict z, 
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

/* Block Order, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void rcmatvmul3x3n ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  // packed selections
  static
    const int16_t ALIGN(32) aSel[BBE_SIMD_WIDTH] = { 0x0200, 0x0f0f, 0x0f01, 0x0f0f, 0x0503, 0x0f0f, 0x0f04, 0x0f0f, 0x0806, 0x0f0f, 0x0f07, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f };
  int l;

  const xb_vecNx16 *          pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  const vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 vTmp;

  vselN vSel0, vSel1;

  xb_vecNx16 vX0, vX1, vY, vXsel, vYsel, vZ;

  xb_vecNx40 wvZ;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % (BBE_SIMD_WIDTH / 8) == 0);

  if (L<=0) return;
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel, 0);
  vSel0 = BBE_MOVVSV(vTmp, 0);
  vSel1 = BBE_MOVVSV(vTmp, 8);

  const vboolN bLoad = BBE_LTRNI(15);

  __Pragma("loop_count min=1");
  for (l = 0; l<L; l += 2)
  {
    /* Load 2 input matrices X */
    BBE_LVNX16T_IP(vX0, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);
    BBE_LVNX16T_IP(vX1, pXrd, 2 * BBE_SIMD_WIDTH, bLoad);

    /* Load 2 input vectors vY */
    BBE_LVNX16_IP(vY, pYrd, 2 * BBE_SIMD_WIDTH);

    vYsel = BBE_SHFLNX16I(vY, BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_LOW_HALF);
    vXsel = BBE_SHFLNX16(vX0, vSel0);
    wvZ = BBE_MULRNX16PC_0(vXsel, vYsel, q);

    vYsel = BBE_SHFLNX16I(vY, BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_LOW_HALF);
    vXsel = BBE_SHFLNX16(vX0, vSel1);
    BBE_MULANX16PC_0(wvZ, vXsel, vYsel);

    vYsel = BBE_SHFLNX16I(vY, BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_HIGH_HALF);
    vXsel = BBE_SHFLNX16(vX1, vSel0);
    BBE_MULANX16PC_1(wvZ, vXsel, vYsel);

    vYsel = BBE_SHFLNX16I(vY, BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_HIGH_HALF);
    vXsel = BBE_SHFLNX16(vX1, vSel1);
    BBE_MULANX16PC_1(wvZ, vXsel, vYsel);

    /* Pack and save results */
    wvZ = BBE_SHFLNX40I(wvZ, BBE_W_SHFLI_DITLV_2);
    vZ = BBE_PACKVNX40(wvZ, q);
    BBE_SVNX16_IP(vZ, pZwr, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatvmul3x3n() */
#endif


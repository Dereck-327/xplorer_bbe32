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

/* Block Order, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     None
*/
void rcmatmul16x16n ( complex_fract16 * restrict z, 
                const int16_t * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q )
{
  int k, i, ystride;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py0 = (const xb_vecNx16 *)y;
  const xb_vecNx16 * restrict py1 = py0 + 1;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 x_, x0, x1, y0, z0, x_sel;
  xb_vecNx40 a0, a1;

  xb_vecNx16 zero = 0;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;

  for (i = k = 0; k<L * 16; k++)
  {
    i = BBE_ADDMOD16U(i, 0x100001);    /*i=(i+1)&15;*/
    ystride = -30 * 2 * BBE_SIMD_WIDTH;
    /* ystride= (i==0)? 2*2*BBE_SIMD_WIDTH:-30*2*BBE_SIMD_WIDTH; */
    XT_MOVEQZ(ystride, 2 * 2 * BBE_SIMD_WIDTH, i);

    BBE_LVNX16_IP(x_, px, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(x1, x0, zero, x_, BBE_DSELI_INTERLEAVE_1);

    x_sel = BBE_REPNX16C(x0, 0); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); a0 = BBE_MULRNX16C(x_sel, y0, q);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); a1 = BBE_MULRNX16C(x_sel, y0, q);
    x_sel = BBE_REPNX16C(x0, 1); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 2); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 3); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 4); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 5); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 6); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x0, 7); BBE_LVNX16_XP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_XP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);

    x_sel = BBE_REPNX16C(x1, 0); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 1); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 2); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 3); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 4); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 5); BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 6); BBE_LVNX16_IP(y0, py0, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_IP(y0, py1, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16C(a1, x_sel, y0);
    x_sel = BBE_REPNX16C(x1, 7); BBE_LVNX16_XP(y0, py0, ystride); BBE_MULANX16C(a0, x_sel, y0);
    BBE_LVNX16_XP(y0, py1, ystride); BBE_MULANX16C(a1, x_sel, y0);

    z0 = BBE_PACKVNX40(a0, q); BBE_SVNX16_IP(z0, pz, 2 * BBE_SIMD_WIDTH);
    z0 = BBE_PACKVNX40(a1, q); BBE_SVNX16_IP(z0, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatmul16x16n() */

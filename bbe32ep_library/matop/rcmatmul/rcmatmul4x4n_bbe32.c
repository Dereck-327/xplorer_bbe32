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
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, rcmatmul4x4n,(complex_fract16 * restrict z, 
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

/* Block Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void rcmatmul4x4n ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  int l;

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 X, X0, Y0, X1, Y1, Z;
  xb_vecNx16 sel_x0, sel_y0, sel_x1, sel_y1;

  xb_vecNx16 zero = 0;

  xb_vecNx40 r;

  vsaN  q = BBE_MOVVSA32(Q);

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L<=0) return;

  __Pragma("loop_count min=1");
  __Pragma("ymemory( py )");
  for (l = 0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVNX16_IP(X, px, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X1, X0, zero, X, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, py, 2 * BBE_SIMD_WIDTH);

    sel_x0 = BBE_SHFLNX16I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_1);
    sel_y0 = BBE_SHFLNX16I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    r = BBE_MULRNX16C(sel_x0, sel_y0, q);
    sel_x1 = BBE_SHFLNX16I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_2);
    sel_y1 = BBE_SHFLNX16I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
    BBE_MULANX16C(r, sel_x1, sel_y1);
    sel_x0 = BBE_SHFLNX16I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_3);
    sel_y0 = BBE_SHFLNX16I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    BBE_MULANX16C(r, sel_x0, sel_y0);
    sel_x1 = BBE_SHFLNX16I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_4);
    sel_y1 = BBE_SHFLNX16I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
    BBE_MULANX16C(r, sel_x1, sel_y1);

    /* Pack and save results (first 2 rows) */
    Z = BBE_PACKVNX40(r, q);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);

    sel_x0 = BBE_SHFLNX16I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_1);
    sel_y0 = BBE_SHFLNX16I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    r = BBE_MULRNX16C(sel_x0, sel_y0, q);
    sel_x1 = BBE_SHFLNX16I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_2);
    sel_y1 = BBE_SHFLNX16I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
    BBE_MULANX16C(r, sel_x1, sel_y1);
    sel_x0 = BBE_SHFLNX16I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_3);
    sel_y0 = BBE_SHFLNX16I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    BBE_MULANX16C(r, sel_x0, sel_y0);
    sel_x1 = BBE_SHFLNX16I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_4);
    sel_y1 = BBE_SHFLNX16I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
    BBE_MULANX16C(r, sel_x1, sel_y1);

    /* Pack and save results (second 2 rows) */
    Z = BBE_PACKVNX40(r, q);
    BBE_SVNX16_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
  }
} /* rcmatmul4x4n() */
#endif


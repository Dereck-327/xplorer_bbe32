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

#if !(HAVE_VFPU)
DISCARD_FUN(void, matmul4x4sf,( float32_t * restrict z, 
                          const float32_t * restrict x, 
                          const float32_t * restrict y, 
                          int L ))
#else

#define sz_f32 ((int)sizeof(float32_t))
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

/* Streaming Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 8
*/
void matmul4x4sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int stridexz, stridey;
  int k, modinc;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x03;
  xb_vecN_2xf32 y00, y01, y02, y03,
                y10, y11, y12, y13,
                y20, y21, y22, y23,
                y30, y31, y32, y33;
  xb_vecN_2xf32 z00, z01, z02, z03;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  px = (const xb_vecN_2xf32 *)(x);
  py = (const xb_vecN_2xf32 *)(y);
  pz = (      xb_vecN_2xf32 *)(z);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/2);

  __Pragma("loop_count factor=4");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1))*4; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridexz = 2*BBE_SIMD_WIDTH;
    XT_MOVNEZ(stridexz, -3*L*sz_f32+2*BBE_SIMD_WIDTH, k);
    stridey = -15*L*sz_f32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridey, -15*L*sz_f32+2*BBE_SIMD_WIDTH-L*sz_f32, k);

    /* Load 1 row of input matrix X */
    BBE_LVN_2XF32_XP(x00, px, L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px, L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px, L*sz_f32);
    x03 = BBE_LVN_2XF32_I(px, 0);
    px = (const xb_vecN_2xf32 *)((intptr_t)px+stridexz);

    /* Load input matrix Y */
    BBE_LVN_2XF32_XP(y00, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y01, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y02, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y03, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y10, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y11, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y12, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y13, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y20, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y21, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y22, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y23, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y30, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y31, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y32, py, L*sz_f32);
    BBE_LVN_2XF32_XP(y33, py, stridey);

    /* Compute 1 row of matrix Z */
    z00 = BBE_MULN_2XF32(x00, y00);
    z01 = BBE_MULN_2XF32(x00, y01);
    z02 = BBE_MULN_2XF32(x00, y02);
    z03 = BBE_MULN_2XF32(x00, y03);
    BBE_MULAN_2XF32(z00, x01, y10);
    BBE_MULAN_2XF32(z01, x01, y11);
    BBE_MULAN_2XF32(z02, x01, y12);
    BBE_MULAN_2XF32(z03, x01, y13);
    BBE_MULAN_2XF32(z00, x02, y20);
    BBE_MULAN_2XF32(z01, x02, y21);
    BBE_MULAN_2XF32(z02, x02, y22);
    BBE_MULAN_2XF32(z03, x02, y23);
    BBE_MULAN_2XF32(z00, x03, y30);
    BBE_MULAN_2XF32(z01, x03, y31);
    BBE_MULAN_2XF32(z02, x03, y32);
    BBE_MULAN_2XF32(z03, x03, y33);

    /* Save results (1 row) */
    BBE_SVN_2XF32_XP(z00, pz, L*sz_f32);
    BBE_SVN_2XF32_XP(z01, pz, L*sz_f32);
    BBE_SVN_2XF32_XP(z02, pz, L*sz_f32);
    BBE_SVN_2XF32_XP(z03, pz, stridexz);
  }
} /* matmul4x4sf() */
#endif

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
DISCARD_FUN(void, matmul3x3sf,( float32_t * restrict z, 
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

/* Streaming Order, Floating-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 8
*/
void matmul3x3sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict px2;
  const xb_vecN_2xf32 * restrict py0;
  const xb_vecN_2xf32 * restrict py1;
  const xb_vecN_2xf32 * restrict py2;
        xb_vecN_2xf32 * restrict pz0;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x10, x11, x12, x20, x21, x22;
  xb_vecN_2xf32 y00, y01, y02, y10, y11, y12, y20, y21, y22;
  xb_vecN_2xf32 z00, z01, z02, z10, z11, z12, z20, z21, z22;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  px0 = (const xb_vecN_2xf32 *)(x+L*0);
  px1 = (const xb_vecN_2xf32 *)(x+L*1);
  px2 = (const xb_vecN_2xf32 *)(x+L*2);
  py0 = (const xb_vecN_2xf32 *)(y+L*0);
  py1 = (const xb_vecN_2xf32 *)(y+L*1);
  py2 = (const xb_vecN_2xf32 *)(y+L*2);
  pz0 = (      xb_vecN_2xf32 *)(z);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
  {
    /* Load input matrix X */
    BBE_LVN_2XF32_XP(x00, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px0, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x21, px1, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x22, px2, -6*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_2XF32_XP(y00, py0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y01, py1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y02, py2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y10, py0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y11, py1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y12, py2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y20, py0, -6*L*sz_f32);
    BBE_LVN_2XF32_XP(y21, py1, -6*L*sz_f32);
    BBE_LVN_2XF32_XP(y22, py2, -6*L*sz_f32);

    /* Compute 2 rows of matrix Z */
    z00 = BBE_MULN_2XF32(x00, y00);
    BBE_MULAN_2XF32(z00, x01, y10);
    BBE_MULAN_2XF32(z00, x02, y20);
    z01 = BBE_MULN_2XF32(x00, y01);
    BBE_MULAN_2XF32(z01, x01, y11);
    BBE_MULAN_2XF32(z01, x02, y21);
    z02 = BBE_MULN_2XF32(x00, y02);
    BBE_MULAN_2XF32(z02, x01, y12);
    BBE_MULAN_2XF32(z02, x02, y22);

    z10 = BBE_MULN_2XF32(x10, y00);
    BBE_MULAN_2XF32(z10, x11, y10);
    BBE_MULAN_2XF32(z10, x12, y20);
    z11 = BBE_MULN_2XF32(x10, y01);
    BBE_MULAN_2XF32(z11, x11, y11);
    BBE_MULAN_2XF32(z11, x12, y21);
    z12 = BBE_MULN_2XF32(x10, y02);
    BBE_MULAN_2XF32(z12, x11, y12);
    BBE_MULAN_2XF32(z12, x12, y22);

    /* Reload matrix Y for better loop scheduling */
    BBE_LVN_2XF32_XP(y00, py0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y01, py1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y02, py2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y10, py0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y11, py1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y12, py2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(y20, py0, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y21, py1, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y22, py2, -6*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Compute 3rd row of matrix Z */
    z20 = BBE_MULN_2XF32(x20, y00);
    BBE_MULAN_2XF32(z20, x21, y10);
    BBE_MULAN_2XF32(z20, x22, y20);
    z21 = BBE_MULN_2XF32(x20, y01);
    BBE_MULAN_2XF32(z21, x21, y11);
    BBE_MULAN_2XF32(z21, x22, y21);
    z22 = BBE_MULN_2XF32(x20, y02);
    BBE_MULAN_2XF32(z22, x21, y12);
    BBE_MULAN_2XF32(z22, x22, y22);

    /* Save results */
    BBE_SVN_2XF32_XP(z00, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z01, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z02, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z10, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z11, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z12, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z20, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z21, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z22, pz0, -8*L*sz_f32+2*BBE_SIMD_WIDTH);
  }
} /* matmul3x3sf() */
#endif

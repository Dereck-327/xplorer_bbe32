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
DISCARD_FUN(void, matvmul8x8sf,( float32_t * restrict z, 
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

/* Streaming Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 8
*/
void matvmul8x8sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L )
{
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int stridex, stridey, stridez;
  int k, modinc;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x03, x04, x05, x06, x07,
                x10, x11, x12, x13, x14, x15, x16, x17,
                x20, x21, x22, x23, x24, x25, x26, x27,
                x30, x31, x32, x33, x34, x35, x36, x37;
  xb_vecN_2xf32 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecN_2xf32 z0, z1, z2, z3;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/2);

  px = (const xb_vecN_2xf32 *)(x);
  py = (const xb_vecN_2xf32 *)(y);
  pz = (      xb_vecN_2xf32 *)(z);

  __Pragma("loop_count factor=2");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1))*2; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridex = -31*L*sz_f32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, 2*BBE_SIMD_WIDTH, k);
    stridey = L*sz_f32;
    XT_MOVEQZ(stridey, k, k);
    stridez = -3*L*sz_f32;
    XT_MOVEQZ(stridez, k, k);

    /* Load 4 rows of input matrix X */
    BBE_LVN_2XF32_XP(x00, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x30, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x21, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x31, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x22, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x32, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x03, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x13, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x23, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x33, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x04, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x14, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x24, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x34, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x05, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x15, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x25, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x35, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x06, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x16, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x26, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x36, px, -23*L*sz_f32);
    BBE_LVN_2XF32_XP(x07, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x17, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x27, px,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(x37, px,   stridex);

    /* Load input matrix Y */
    BBE_LVN_2XF32_XP(y0, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y1, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y2, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y3, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y4, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y5, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y6, py,    L*sz_f32);
    BBE_LVN_2XF32_XP(y7, py, -8*L*sz_f32+2*BBE_SIMD_WIDTH);
    py = (const xb_vecN_2xf32 *)((intptr_t)py + stridey);

    /* Compute 4 rows of matrix Z */
    z0 = BBE_MULN_2XF32(x00, y0);
    z1 = BBE_MULN_2XF32(x10, y0);
    z2 = BBE_MULN_2XF32(x20, y0);
    z3 = BBE_MULN_2XF32(x30, y0);
    BBE_MULAN_2XF32(z0, x01, y1);
    BBE_MULAN_2XF32(z1, x11, y1);
    BBE_MULAN_2XF32(z2, x21, y1);
    BBE_MULAN_2XF32(z3, x31, y1);
    BBE_MULAN_2XF32(z0, x02, y2);
    BBE_MULAN_2XF32(z1, x12, y2);
    BBE_MULAN_2XF32(z2, x22, y2);
    BBE_MULAN_2XF32(z3, x32, y2);
    BBE_MULAN_2XF32(z0, x03, y3);
    BBE_MULAN_2XF32(z1, x13, y3);
    BBE_MULAN_2XF32(z2, x23, y3);
    BBE_MULAN_2XF32(z3, x33, y3);
    BBE_MULAN_2XF32(z0, x04, y4);
    BBE_MULAN_2XF32(z1, x14, y4);
    BBE_MULAN_2XF32(z2, x24, y4);
    BBE_MULAN_2XF32(z3, x34, y4);
    BBE_MULAN_2XF32(z0, x05, y5);
    BBE_MULAN_2XF32(z1, x15, y5);
    BBE_MULAN_2XF32(z2, x25, y5);
    BBE_MULAN_2XF32(z3, x35, y5);
    BBE_MULAN_2XF32(z0, x06, y6);
    BBE_MULAN_2XF32(z1, x16, y6);
    BBE_MULAN_2XF32(z2, x26, y6);
    BBE_MULAN_2XF32(z3, x36, y6);
    BBE_MULAN_2XF32(z0, x07, y7);
    BBE_MULAN_2XF32(z1, x17, y7);
    BBE_MULAN_2XF32(z2, x27, y7);
    BBE_MULAN_2XF32(z3, x37, y7);

    /* Save results (4 values) */
    BBE_SVN_2XF32_XP(z0, pz, L*sz_f32);
    BBE_SVN_2XF32_XP(z1, pz, L*sz_f32);
    BBE_SVN_2XF32_XP(z2, pz, L*sz_f32);
    BBE_SVN_2XF32_IP(z3, pz, 2*BBE_SIMD_WIDTH);
    pz = (xb_vecN_2xf32 *)((intptr_t)pz + stridez);
  }
} /* matvmul8x8sf() */
#endif

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
DISCARD_FUN(void, matmul8x8sf,( float32_t * restrict z, 
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

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 8
*/
void matmul8x8sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int stridex;
  int k, modinc;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x03, x04, x05, x06, x07,
                x10, x11, x12, x13, x14, x15, x16, x17,
                x20, x21, x22, x23, x24, x25, x26, x27,
                x30, x31, x32, x33, x34, x35, x36, x37;
  xb_vecN_2xf32 y00, y10, y20, y30, y40, y50, y60, y70;
  xb_vecN_2xf32 z00, z10, z20, z30;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/2);

  /*
   * Compute part of matrices (first 4 rows)
   */
  px = (const xb_vecN_2xf32 *)(x);
  py = (const xb_vecN_2xf32 *)(y);
  pz = (      xb_vecN_2xf32 *)(z);

  __Pragma("loop_count factor=4");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1))*8; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridex = -31*L*sz_f32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, -31*L*sz_f32+2*BBE_SIMD_WIDTH-L*sz_f32, k);

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

    /* Load 1 column of input matrix Y */
    BBE_LVN_2XF32_XP(y00, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y10, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y20, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y30, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y40, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y50, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y60, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y70, py, -56*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Compute 4 values for 4 rows of matrix Z */
    z00 = BBE_MULN_2XF32(x00, y00);
    z10 = BBE_MULN_2XF32(x10, y00);
    z20 = BBE_MULN_2XF32(x20, y00);
    z30 = BBE_MULN_2XF32(x30, y00);
    BBE_MULAN_2XF32(z00, x01, y10);
    BBE_MULAN_2XF32(z10, x11, y10);
    BBE_MULAN_2XF32(z20, x21, y10);
    BBE_MULAN_2XF32(z30, x31, y10);
    BBE_MULAN_2XF32(z00, x02, y20);
    BBE_MULAN_2XF32(z10, x12, y20);
    BBE_MULAN_2XF32(z20, x22, y20);
    BBE_MULAN_2XF32(z30, x32, y20);
    BBE_MULAN_2XF32(z00, x03, y30);
    BBE_MULAN_2XF32(z10, x13, y30);
    BBE_MULAN_2XF32(z20, x23, y30);
    BBE_MULAN_2XF32(z30, x33, y30);
    BBE_MULAN_2XF32(z00, x04, y40);
    BBE_MULAN_2XF32(z10, x14, y40);
    BBE_MULAN_2XF32(z20, x24, y40);
    BBE_MULAN_2XF32(z30, x34, y40);
    BBE_MULAN_2XF32(z00, x05, y50);
    BBE_MULAN_2XF32(z10, x15, y50);
    BBE_MULAN_2XF32(z20, x25, y50);
    BBE_MULAN_2XF32(z30, x35, y50);
    BBE_MULAN_2XF32(z00, x06, y60);
    BBE_MULAN_2XF32(z10, x16, y60);
    BBE_MULAN_2XF32(z20, x26, y60);
    BBE_MULAN_2XF32(z30, x36, y60);
    BBE_MULAN_2XF32(z00, x07, y70);
    BBE_MULAN_2XF32(z10, x17, y70);
    BBE_MULAN_2XF32(z20, x27, y70);
    BBE_MULAN_2XF32(z30, x37, y70);

    /* Save results (by one values for 4 rows) */
    BBE_SVN_2XF32_XP(z00, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z10, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z20, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z30, pz, -24*L*sz_f32+2*BBE_SIMD_WIDTH);
  }
  
  /*
   * Compute next part of matrices (second 4 rows)
   */
  px = (const xb_vecN_2xf32 *)(x+L*8*4);
  py = (const xb_vecN_2xf32 *)(y);
  pz = (      xb_vecN_2xf32 *)(z+L*8*4);

  __Pragma("loop_count factor=4");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1))*8; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridex = -31*L*sz_f32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, -31*L*sz_f32+2*BBE_SIMD_WIDTH-L*sz_f32, k);

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

    /* Load 1 column of input matrix Y */
    BBE_LVN_2XF32_XP(y00, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y10, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y20, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y30, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y40, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y50, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y60, py,   8*L*sz_f32);
    BBE_LVN_2XF32_XP(y70, py, -56*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Compute 4 values for 4 rows of matrix Z */
    z00 = BBE_MULN_2XF32(x00, y00);
    z10 = BBE_MULN_2XF32(x10, y00);
    z20 = BBE_MULN_2XF32(x20, y00);
    z30 = BBE_MULN_2XF32(x30, y00);
    BBE_MULAN_2XF32(z00, x01, y10);
    BBE_MULAN_2XF32(z10, x11, y10);
    BBE_MULAN_2XF32(z20, x21, y10);
    BBE_MULAN_2XF32(z30, x31, y10);
    BBE_MULAN_2XF32(z00, x02, y20);
    BBE_MULAN_2XF32(z10, x12, y20);
    BBE_MULAN_2XF32(z20, x22, y20);
    BBE_MULAN_2XF32(z30, x32, y20);
    BBE_MULAN_2XF32(z00, x03, y30);
    BBE_MULAN_2XF32(z10, x13, y30);
    BBE_MULAN_2XF32(z20, x23, y30);
    BBE_MULAN_2XF32(z30, x33, y30);
    BBE_MULAN_2XF32(z00, x04, y40);
    BBE_MULAN_2XF32(z10, x14, y40);
    BBE_MULAN_2XF32(z20, x24, y40);
    BBE_MULAN_2XF32(z30, x34, y40);
    BBE_MULAN_2XF32(z00, x05, y50);
    BBE_MULAN_2XF32(z10, x15, y50);
    BBE_MULAN_2XF32(z20, x25, y50);
    BBE_MULAN_2XF32(z30, x35, y50);
    BBE_MULAN_2XF32(z00, x06, y60);
    BBE_MULAN_2XF32(z10, x16, y60);
    BBE_MULAN_2XF32(z20, x26, y60);
    BBE_MULAN_2XF32(z30, x36, y60);
    BBE_MULAN_2XF32(z00, x07, y70);
    BBE_MULAN_2XF32(z10, x17, y70);
    BBE_MULAN_2XF32(z20, x27, y70);
    BBE_MULAN_2XF32(z30, x37, y70);

    /* Save results (by one values for 4 rows) */
    BBE_SVN_2XF32_XP(z00, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z10, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z20, pz,   8*L*sz_f32);
    BBE_SVN_2XF32_XP(z30, pz, -24*L*sz_f32+2*BBE_SIMD_WIDTH);
  }

} /* matmul8x8sf() */
#endif

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
DISCARD_FUN(void, matmul8x8nf,( float32_t * restrict z, 
                          const float32_t * restrict x, 
                          const float32_t * restrict y, 
                          int L ))
#else

#ifndef BBE_MOVN_2XF32_FROMN_4XCF32
#define BBE_MOVN_2XF32_FROMN_4XCF32(a) BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a))
#endif
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

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void matmul8x8nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
#if 0
  #define BBE_MOVN_2XF32_FROMN_4XCF32(a) BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a))
  int k, i, ystride;

  const xb_vecN_4xcf32 * restrict px = (const xb_vecN_4xcf32 *)x;
  const xb_vecN_2xf32  * restrict py = (const xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pz = (      xb_vecN_2xf32  *)z;
  xb_vecN_2xf32 x_sel01, x_sel23, x_sel45, x_sel67;
  xb_vecN_2xf32 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecN_2xf32 acc;
  xb_vecN_4xcf32 tmp;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  if (L<=0) return;

  __Pragma("loop_count min=2, factor=2");
  for (i = k = 0; k < L*2; k++)
  {
    /* i=(i+4)&7; */
    i = BBE_ADDMOD16U(i, 0x080004);   
    ystride = -7 * 2 * BBE_SIMD_WIDTH;
    /* ystride= (i==0)? 2*BBE_SIMD_WIDTH:-7*2*BBE_SIMD_WIDTH; */
    XT_MOVEQZ(ystride, 2 * BBE_SIMD_WIDTH, i);

    BBE_LVN_2XF32_IP(y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y2, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y3, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y4, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y5, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y6, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y7, py, ystride);

    BBE_LVN_4XCF32_IP(tmp, px, 2 * BBE_SIMD_WIDTH);
    x_sel01 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 0));
    x_sel23 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 1));
    x_sel45 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 2));
    x_sel67 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 3));
    acc = BBE_MULMN_2XF32( x_sel01, y0, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel01, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel23, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel23, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel45, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel45, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel67, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel67, y7, 0, 0xE);
    BBE_SVN_2XF32_IP(acc, pz, 2 * BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(tmp, px, 2 * BBE_SIMD_WIDTH);
    x_sel01 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 0));
    x_sel23 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 1));
    x_sel45 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 2));
    x_sel67 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 3));
    acc = BBE_MULMN_2XF32( x_sel01, y0, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel01, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel23, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel23, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel45, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel45, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel67, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel67, y7, 0, 0xE);
    BBE_SVN_2XF32_IP(acc, pz, 2 * BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(tmp, px, 2 * BBE_SIMD_WIDTH);
    x_sel01 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 0));
    x_sel23 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 1));
    x_sel45 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 2));
    x_sel67 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 3));
    acc = BBE_MULMN_2XF32( x_sel01, y0, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel01, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel23, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel23, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel45, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel45, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel67, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel67, y7, 0, 0xE);
    BBE_SVN_2XF32_IP(acc, pz, 2 * BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(tmp, px, 2 * BBE_SIMD_WIDTH);
    x_sel01 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 0));
    x_sel23 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 1));
    x_sel45 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 2));
    x_sel67 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp, 3));
    acc = BBE_MULMN_2XF32( x_sel01, y0, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel01, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel23, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel23, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel45, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel45, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc, x_sel67, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc, x_sel67, y7, 0, 0xE);
    BBE_SVN_2XF32_IP(acc, pz, 2 * BBE_SIMD_WIDTH);
  }
#else
  int k, i, ystride;

  const xb_vecN_4xcf32 * restrict px0 = (const xb_vecN_4xcf32 *)x;
  const xb_vecN_4xcf32 * restrict px1 = (const xb_vecN_4xcf32 *)x;
  const xb_vecN_2xf32  * restrict py  = (const xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pz  = (      xb_vecN_2xf32  *)z;
  xb_vecN_2xf32 x_sel001, x_sel023, x_sel045, x_sel067;
  xb_vecN_2xf32 x_sel101, x_sel123, x_sel145, x_sel167;
  xb_vecN_2xf32 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecN_2xf32 acc00, acc01, acc10, acc11;
  xb_vecN_2xf32 z0, z1;
  xb_vecN_4xcf32 tmp0, tmp1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  if (L<=0) return;

  __Pragma("loop_count min=4, factor=4");
  for (i = k = 0; k < L*4; k++)
  {
    /* i=(i+2)&7; */
    i = BBE_ADDMOD16U(i, 0x080002);   
    ystride = -8 * 2 * BBE_SIMD_WIDTH;
    /* ystride= (i==0)? 0 : -8*2*BBE_SIMD_WIDTH; */
    XT_MOVEQZ(ystride, i, i);

    BBE_LVN_2XF32_XP(y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y1, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y2, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y3, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y4, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y5, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y6, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(y7, py, 2 * BBE_SIMD_WIDTH);
    py = (const xb_vecN_2xf32 *)((intptr_t)py + ystride);

    BBE_LVN_4XCF32_IP(tmp0, px0, 2 * BBE_SIMD_WIDTH);
    x_sel001 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
    x_sel023 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 1));
    BBE_LVN_4XCF32_IP(tmp0, px0, 2 * BBE_SIMD_WIDTH);
    x_sel101 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
    x_sel123 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 1));
    BBE_LVN_4XCF32_IP(tmp1, px1, 2 * BBE_SIMD_WIDTH);
    x_sel045 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 2));
    x_sel067 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 3));
    BBE_LVN_4XCF32_IP(tmp1, px1, 2 * BBE_SIMD_WIDTH);
    x_sel145 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 2));
    x_sel167 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 3));

    acc00 = BBE_MULMN_2XF32( x_sel001, y0, 0, 0x4);
    acc01 = BBE_MULMN_2XF32( x_sel001, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc00, x_sel023, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc01, x_sel023, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc00, x_sel045, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc01, x_sel045, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc00, x_sel067, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc01, x_sel067, y7, 0, 0xE);

    acc10 = BBE_MULMN_2XF32( x_sel101, y0, 0, 0x4);
    acc11 = BBE_MULMN_2XF32( x_sel101, y1, 0, 0xE);
    BBE_MULMASN_2XF32(acc10, x_sel123, y2, 0, 0x4);
    BBE_MULMASN_2XF32(acc11, x_sel123, y3, 0, 0xE);
    BBE_MULMASN_2XF32(acc10, x_sel145, y4, 0, 0x4);
    BBE_MULMASN_2XF32(acc11, x_sel145, y5, 0, 0xE);
    BBE_MULMASN_2XF32(acc10, x_sel167, y6, 0, 0x4);
    BBE_MULMASN_2XF32(acc11, x_sel167, y7, 0, 0xE);

    z0 = BBE_ADDN_2XF32(acc00, acc01);
    z1 = BBE_ADDN_2XF32(acc10, acc11);

    BBE_SVN_2XF32_IP(z0, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(z1, pz, 2 * BBE_SIMD_WIDTH);
  }
#endif
} /* matmul8x8nf() */
#endif

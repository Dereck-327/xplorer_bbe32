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
    Complex Matrix-Matrix/Matrix-Vector Multiply
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
DISCARD_FUN(void, cmatmul8x8nf,( complex_float * restrict z, 
                           const complex_float * restrict x, 
                           const complex_float * restrict y, 
                           int L ))
#else
/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
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

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void cmatmul8x8nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L )
{
#if 0
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;
  int i, ystride;

  xb_vecN_4xcf32 X0;
  xb_vecN_4xcf32 Y00, Y01, Y10, Y11, Y20, Y21, Y30, Y31;
  xb_vecN_4xcf32 x00, x01, x02, x03;
  xb_vecN_4xcf32 Z0, Z1;
  xb_vecN_4xcf32 z00, z01, z10, z11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x);
  py  = (const xb_vecN_4xcf32 *)(y);
  pz  = (      xb_vecN_4xcf32 *)(z);
  i = 0;

  __Pragma("loop_count factor=2");
  for (l=0; l<L*8; l++)
  {
    /* i=(i+1)&7; */
    i = BBE_ADDMOD16U(i, 0x080001); 
    /* ystride= (i==0)? 0 : -4*2*BBE_SIMD_WIDTH; */
    ystride = -16 * 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(ystride, i, i);

    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(Y00, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y01, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y10, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y11, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y20, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y21, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y30, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y31, py, 2*BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    x00 = BBE_REPN_4XCF32(X0, 0);
    x01 = BBE_REPN_4XCF32(X0, 1);
    //BBE_LVN_4XCF32_IP(X0, px1, 2*BBE_SIMD_WIDTH);
    x02 = BBE_REPN_4XCF32(X0, 2);
    x03 = BBE_REPN_4XCF32(X0, 3);

    z00 = BBE_MULMN_4XCF32( x00, Y00, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x00, Y00, 1, 0xB);
    z10 = BBE_MULMN_4XCF32( x00, Y01, 0, 0x4);
    z11 = BBE_MULMN_4XCF32( x00, Y01, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, Y10, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, Y10, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x01, Y11, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x01, Y11, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, Y20, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, Y20, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x02, Y21, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x02, Y21, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x03, Y30, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x03, Y30, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x03, Y31, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x03, Y31, 1, 0xB);

    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(Y00, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y01, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y10, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y11, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y20, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y21, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y30, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y31, py, 2*BBE_SIMD_WIDTH);
    py = (const xb_vecN_4xcf32 *)((intptr_t)py + ystride);

    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    x00 = BBE_REPN_4XCF32(X0, 0);
    x01 = BBE_REPN_4XCF32(X0, 1);
    //BBE_LVN_4XCF32_IP(X0, px1, 2*BBE_SIMD_WIDTH);
    x02 = BBE_REPN_4XCF32(X0, 2);
    x03 = BBE_REPN_4XCF32(X0, 3);

    BBE_MULMASN_4XCF32(z00, x00, Y00, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x00, Y00, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x00, Y01, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x00, Y01, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, Y10, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, Y10, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x01, Y11, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x01, Y11, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, Y20, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, Y20, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x02, Y21, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x02, Y21, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x03, Y30, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x03, Y30, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x03, Y31, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x03, Y31, 1, 0xB);
    Z0 = BBE_ADDN_4XCF32(z00, z01);
    Z1 = BBE_ADDN_4XCF32(z10, z11);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z1, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;
  int i, j, rowstride, colstride;

  xb_vecN_4xcf32 X0, X1;
  xb_vecN_4xcf32 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  xb_vecN_4xcf32 x0, x1, x2, x3, x4, x5, x6, x7;
  xb_vecN_4xcf32 Z0;
  xb_vecN_4xcf32 z00, z01, z02, z03;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_4xcf32 *)(x);
  py = (const xb_vecN_4xcf32 *)(y);
  pz = (      xb_vecN_4xcf32 *)(z);
  i = j = 0;

  __Pragma("loop_count factor=4");
  for (l=0; l<L*16; l++)
  {
    j = XT_XOR(j, 1);
    /* i=(i+1)%16; */
    i = BBE_ADDMOD16U(i, 0x100001);
    /* rowstride= (j==0)? 0 : 2*2*BBE_SIMD_WIDTH; */
    rowstride = 2*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(rowstride, j, j);
    /* colstride= (i==0)? 0 : -16*2*BBE_SIMD_WIDTH; */
    colstride = -16*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(colstride, i, i);

    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(Y0, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y1, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y2, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y3, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y4, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y5, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y6, py, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y7, py, 1*2*BBE_SIMD_WIDTH);
    py = (const xb_vecN_4xcf32 *)((intptr_t)py + rowstride + colstride);

    BBE_LVN_4XCF32_IP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px, 2*BBE_SIMD_WIDTH);
    px = (const xb_vecN_4xcf32 *)((intptr_t)px - rowstride);
    x0 = BBE_REPN_4XCF32(X0, 0);
    x1 = BBE_REPN_4XCF32(X0, 1);
    x2 = BBE_REPN_4XCF32(X0, 2);
    x3 = BBE_REPN_4XCF32(X0, 3);
    x4 = BBE_REPN_4XCF32(X1, 0);
    x5 = BBE_REPN_4XCF32(X1, 1);
    x6 = BBE_REPN_4XCF32(X1, 2);
    x7 = BBE_REPN_4XCF32(X1, 3);

    z00 = BBE_MULMN_4XCF32( x0, Y0, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x0, Y0, 1, 0xB);
    z02 = BBE_MULMN_4XCF32( x1, Y1, 0, 0x4);
    z03 = BBE_MULMN_4XCF32( x1, Y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x2, Y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x2, Y2, 1, 0xB);
    BBE_MULMASN_4XCF32(z02, x3, Y3, 0, 0x4);
    BBE_MULMASN_4XCF32(z03, x3, Y3, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x4, Y4, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x4, Y4, 1, 0xB);
    BBE_MULMASN_4XCF32(z02, x5, Y5, 0, 0x4);
    BBE_MULMASN_4XCF32(z03, x5, Y5, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x6, Y6, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x6, Y6, 1, 0xB);
    BBE_MULMASN_4XCF32(z02, x7, Y7, 0, 0x4);
    BBE_MULMASN_4XCF32(z03, x7, Y7, 1, 0xB);
    z00 = BBE_ADDN_4XCF32(z00, z01);
    z02 = BBE_ADDN_4XCF32(z02, z03);
    Z0 = BBE_ADDN_4XCF32(z00, z02);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatmul8x8nf() */
#endif

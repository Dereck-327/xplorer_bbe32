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
DISCARD_FUN(void, cmatvmul8x8sf,( complex_float * restrict z, 
                            const complex_float * restrict x, 
                            const complex_float * restrict y, 
                            int L ))
#else

#define sz_cf32 ((int)sizeof(complex_float))
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

/* Streaming Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul8x8sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int stridex, stridey, stridez;
  int k, modinc;
  int l;

  xb_vecN_4xcf32 x00, x01, x02, x03, x04, x05, x06, x07,
                 x10, x11, x12, x13, x14, x15, x16, x17;
  xb_vecN_4xcf32 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecN_4xcf32 z00, z01, z10, z11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  px = (const xb_vecN_4xcf32 *)(x);
  py = (const xb_vecN_4xcf32 *)(y);
  pz = (      xb_vecN_4xcf32 *)(z);

  __Pragma("loop_count factor=2");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*4; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridex = -15*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, 2*BBE_SIMD_WIDTH, k);
    stridey = L*sz_cf32;
    XT_MOVEQZ(stridey, k, k);
    stridez = -L*sz_cf32;
    XT_MOVEQZ(stridez, k, k);

    /* Load 2 rows of input matrix X */
    BBE_LVN_4XCF32_XP(x00, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x10, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x11, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x02, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x12, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x03, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x13, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x04, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x14, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x05, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x15, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x06, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x16, px,  -7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x07, px,   8*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x17, px,   stridex);

    /* Load input matrix Y */
    BBE_LVN_4XCF32_XP(y0, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y1, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y2, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y3, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y4, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y5, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y6, py,    L*sz_cf32);
    BBE_LVN_4XCF32_XP(y7, py, -8*L*sz_cf32+2*BBE_SIMD_WIDTH);
    py = (const xb_vecN_4xcf32 *)((intptr_t)py + stridey);

    /* Compute 2 rows of matrix Z */
    z00 = BBE_MULMN_4XCF32( x00, y0, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x00, y0, 1, 0xB);
    z10 = BBE_MULMN_4XCF32( x10, y0, 0, 0x4);
    z11 = BBE_MULMN_4XCF32( x10, y0, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x11, y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x11, y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, y2, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x12, y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x12, y2, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x03, y3, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x03, y3, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x13, y3, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x13, y3, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x04, y4, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x04, y4, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x14, y4, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x14, y4, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x05, y5, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x05, y5, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x15, y5, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x15, y5, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x06, y6, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x06, y6, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x16, y6, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x16, y6, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x07, y7, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x07, y7, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x17, y7, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x17, y7, 1, 0xB);
    /* Make dummy operations for better loop scheduling */
    z00 = BBE_ADDN_4XCF32(z00, 0.0f);
    z10 = BBE_ADDN_4XCF32(z10, 0.0f);
    z01 = BBE_ADDN_4XCF32(z01, 0.0f);
    z11 = BBE_ADDN_4XCF32(z11, 0.0f);

    z00 = BBE_ADDN_4XCF32(z00, z01);
    z10 = BBE_ADDN_4XCF32(z10, z11);

    /* Save results (2 values) */
    BBE_SVN_4XCF32_XP(z00, pz, L*sz_cf32);
    BBE_SVN_4XCF32_IP(z10, pz, 2*BBE_SIMD_WIDTH);
    pz = (xb_vecN_4xcf32 *)((intptr_t)pz + stridez);
  }
} /* cmatvmul8x8sf() */
#endif

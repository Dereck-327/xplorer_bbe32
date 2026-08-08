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
DISCARD_FUN(void, cmatmul3x3sf, ( complex_float * restrict z, 
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

/* Streaming Order, Floating-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul3x3sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L )
{
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int stridexz, stridey;
  int k, modinc;
  int l;

  xb_vecN_4xcf32 x00, x01, x02;
  xb_vecN_4xcf32 y00, y01, y02,
                 y10, y11, y12,
                 y20, y21, y22;
  xb_vecN_4xcf32 z00, z01, z10, z11, z20, z21;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  px = (const xb_vecN_4xcf32 *)(x);
  py = (const xb_vecN_4xcf32 *)(y);
  pz = (      xb_vecN_4xcf32 *)(z);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*3; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridexz = 2*BBE_SIMD_WIDTH;
    XT_MOVNEZ(stridexz, -2*L*sz_cf32+2*BBE_SIMD_WIDTH, k);
    stridey = -8*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridey, -8*L*sz_cf32+2*BBE_SIMD_WIDTH-L*sz_cf32, k);

    /* Load one row of input matrix X */
    BBE_LVN_4XCF32_XP(x00, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px, L*sz_cf32);
    x02 = BBE_LVN_4XCF32_I(px, 0);
    px = (const xb_vecN_4xcf32 *)((intptr_t)px+stridexz);

    /* Load input matrix Y */
    BBE_LVN_4XCF32_XP(y00, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y01, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y02, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y10, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y11, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y12, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y20, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y21, py, L*sz_cf32);
    BBE_LVN_4XCF32_XP(y22, py, stridey);

    /* Compute one row of matrix Z */
    z00 = BBE_MULMN_4XCF32( x00, y00, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x00, y00, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, y10, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, y10, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, y20, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, y20, 1, 0xB);
    z00 = BBE_ADDN_4XCF32(z00, z01);

    z10 = BBE_MULMN_4XCF32( x00, y01, 0, 0x4);
    z11 = BBE_MULMN_4XCF32( x00, y01, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x01, y11, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x01, y11, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x02, y21, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x02, y21, 1, 0xB);
    z10 = BBE_ADDN_4XCF32(z10, z11);

    z20 = BBE_MULMN_4XCF32( x00, y02, 0, 0x4);
    z21 = BBE_MULMN_4XCF32( x00, y02, 1, 0xB);
    BBE_MULMASN_4XCF32(z20, x01, y12, 0, 0x4);
    BBE_MULMASN_4XCF32(z21, x01, y12, 1, 0xB);
    BBE_MULMASN_4XCF32(z20, x02, y22, 0, 0x4);
    BBE_MULMASN_4XCF32(z21, x02, y22, 1, 0xB);
    z20 = BBE_ADDN_4XCF32(z20, z21);

    /* Save results */
    BBE_SVN_4XCF32_XP(z00, pz, L*sz_cf32);
    BBE_SVN_4XCF32_XP(z10, pz, L*sz_cf32);
    BBE_SVN_4XCF32_XP(z20, pz, stridexz);
  }
} /* cmatmul3x3sf() */
#endif

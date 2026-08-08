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
DISCARD_FUN(void, cmatvmul4x4sf,( complex_float * restrict z, 
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

/* Streaming Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul4x4sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict px2;
  const xb_vecN_4xcf32 * restrict px3;
  const xb_vecN_4xcf32 * restrict py0;
  const xb_vecN_4xcf32 * restrict py1;
  const xb_vecN_4xcf32 * restrict py2;
  const xb_vecN_4xcf32 * restrict py3;
        xb_vecN_4xcf32 * restrict pz0;
  int l;

  xb_vecN_4xcf32 x00, x01, x02, x03,
                 x10, x11, x12, x13,
                 x20, x21, x22, x23,
                 x30, x31, x32, x33;
  xb_vecN_4xcf32 y0, y1, y2, y3;
  xb_vecN_4xcf32 z0, z1, z2, z3;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  px0 = (const xb_vecN_4xcf32 *)(x+L*0);
  px1 = (const xb_vecN_4xcf32 *)(x+L*1);
  px2 = (const xb_vecN_4xcf32 *)(x+L*2);
  px3 = (const xb_vecN_4xcf32 *)(x+L*3);
  py0 = (const xb_vecN_4xcf32 *)(y+L*0);
  py1 = (const xb_vecN_4xcf32 *)(y+L*1);
  py2 = (const xb_vecN_4xcf32 *)(y+L*2);
  py3 = (const xb_vecN_4xcf32 *)(y+L*3);
  pz0 = (      xb_vecN_4xcf32 *)(z);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_XP(x00, px0,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px1,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x02, px2,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x03, px3,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x10, px0,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x11, px1,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x12, px2,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x13, px3,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x20, px0,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x21, px1,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x22, px2,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x23, px3,   4*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x30, px0, -12*L*sz_cf32+2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(x31, px1, -12*L*sz_cf32+2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(x32, px2, -12*L*sz_cf32+2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(x33, px3, -12*L*sz_cf32+2*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_4XCF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(y3, py3, 2*BBE_SIMD_WIDTH);

    /* Compute matrix Z */
    z0 = BBE_MULN_4XCF32(x00, y0);
    z1 = BBE_MULN_4XCF32(x10, y0);
    z2 = BBE_MULN_4XCF32(x20, y0);
    z3 = BBE_MULN_4XCF32(x30, y0);
    BBE_MULAN_4XCF32(z0, x01, y1);
    BBE_MULAN_4XCF32(z1, x11, y1);
    BBE_MULAN_4XCF32(z2, x21, y1);
    BBE_MULAN_4XCF32(z3, x31, y1);
    BBE_MULAN_4XCF32(z0, x02, y2);
    BBE_MULAN_4XCF32(z1, x12, y2);
    BBE_MULAN_4XCF32(z2, x22, y2);
    BBE_MULAN_4XCF32(z3, x32, y2);
    BBE_MULAN_4XCF32(z0, x03, y3);
    BBE_MULAN_4XCF32(z1, x13, y3);
    BBE_MULAN_4XCF32(z2, x23, y3);
    BBE_MULAN_4XCF32(z3, x33, y3);

    /* Save results */
    BBE_SVN_4XCF32_XP(z0, pz0,    L*sz_cf32);
    BBE_SVN_4XCF32_XP(z1, pz0,    L*sz_cf32);
    BBE_SVN_4XCF32_XP(z2, pz0,    L*sz_cf32);
    BBE_SVN_4XCF32_XP(z3, pz0, -3*L*sz_cf32+2*BBE_SIMD_WIDTH);
  }
} /* cmatvmul4x4sf() */
#endif

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
    Matrix Hermitian Product
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
DISCARD_FUN(void, cmatherm4x4sf,( complex_float * restrict y, 
                            const complex_float * restrict x, 
                            int L ))
#else

#define sz_cf32 ((int)sizeof(complex_float))
/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Streaming Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm4x4sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
#if 0
  const xb_vecN_4xcf32 * restrict pxc;
  const xb_vecN_4xcf32 * restrict px;
        xb_vecN_4xcf32 * restrict py;
  int stridex, stridey;
  int k, modinc;
  int l;

  xb_vecN_4xcf32 xc0, xc1, xc2, xc3;
  xb_vecN_4xcf32 x00, x01, x02, x03,
                 x10, x11, x12, x13,
                 x20, x21, x22, x23,
                 x30, x31, x32, x33;
  xb_vecN_4xcf32 y00, y01, y02, y03;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  pxc = (const xb_vecN_4xcf32 *)(x);
  px  = (const xb_vecN_4xcf32 *)(x);
  py  = (      xb_vecN_4xcf32 *)(y);

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  __Pragma("loop_count factor=4");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*4; l++)
  {
    k = BBE_ADDMOD16U(k, modinc);
    stridey = 2*BBE_SIMD_WIDTH;
    XT_MOVNEZ(stridey, -3*L*sz_cf32+2*BBE_SIMD_WIDTH, k);
    stridex = -15*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, -15*L*sz_cf32+2*BBE_SIMD_WIDTH-L*sz_cf32, k);

    /* Load 1 column of input matrix X */
    BBE_LVN_4XCF32_XP(xc0, pxc, L*4*sz_cf32);
    BBE_LVN_4XCF32_XP(xc1, pxc, L*4*sz_cf32);
    BBE_LVN_4XCF32_XP(xc2, pxc, L*4*sz_cf32);
    BBE_LVN_4XCF32_XP(xc3, pxc, -3*L*4*sz_cf32+2*BBE_SIMD_WIDTH);

    /* Load input matrix x */
    BBE_LVN_4XCF32_XP(x00, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x02, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x03, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x10, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x11, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x12, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x13, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x20, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x21, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x22, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x23, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x30, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x31, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x32, px, L*sz_cf32);
    BBE_LVN_4XCF32_XP(x33, px, stridex);

    /* Compute 1 row of matrix Z */
    y00 = BBE_MULJN_4XCF32(x00, xc0);
    y01 = BBE_MULJN_4XCF32(x01, xc0);
    y02 = BBE_MULJN_4XCF32(x02, xc0);
    y03 = BBE_MULJN_4XCF32(x03, xc0);
    BBE_MULJAN_4XCF32(y00, x10, xc1);
    BBE_MULJAN_4XCF32(y01, x11, xc1);
    BBE_MULJAN_4XCF32(y02, x12, xc1);
    BBE_MULJAN_4XCF32(y03, x13, xc1);
    BBE_MULJAN_4XCF32(y00, x20, xc2);
    BBE_MULJAN_4XCF32(y01, x21, xc2);
    BBE_MULJAN_4XCF32(y02, x22, xc2);
    BBE_MULJAN_4XCF32(y03, x23, xc2);
    BBE_MULJAN_4XCF32(y00, x30, xc3);
    BBE_MULJAN_4XCF32(y01, x31, xc3);
    BBE_MULJAN_4XCF32(y02, x32, xc3);
    BBE_MULJAN_4XCF32(y03, x33, xc3);

    /* Save results (1 row) */
    BBE_SVN_4XCF32_XP(y00, py, L*sz_cf32);
    BBE_SVN_4XCF32_XP(y01, py, L*sz_cf32);
    BBE_SVN_4XCF32_XP(y02, py, L*sz_cf32);
    BBE_SVN_4XCF32_XP(y03, py, stridey);
  }
#else
  const xb_vecN_4xcf32 * restrict px;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 x00, x01, x02, x03,
                 x10, x11, x12, x13,
                 x20, x21, x22, x23,
                 x30, x31, x32, x33;
  xb_vecN_4xcf32 y00, y01, y02, y03,
                 y10, y11, y12, y13,
                 y20, y21, y22, y23,
                 y30, y31, y32, y33;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  /* Compute half of matrix */

  px  = (const xb_vecN_4xcf32 *)(x);
  py  = (      xb_vecN_4xcf32 *)(y);

  __Pragma("loop_count min=1");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
  {
    /* Load input matrix x */
    BBE_LVN_4XCF32_XP(x00, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x02, px,   2*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x10, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x11, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x12, px,   2*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x20, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x21, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x22, px,   2*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x30, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x31, px,     L*sz_cf32);
    BBE_LVN_4XCF32_XP(x32, px, -14*L*sz_cf32+2*BBE_SIMD_WIDTH);

    /* Compute half of matrix Y=X'*X */
    y00 = BBE_MULJN_4XCF32(x00, x00);
    y01 = BBE_MULJN_4XCF32(x01, x00);
    y11 = BBE_MULJN_4XCF32(x01, x01);
    y02 = BBE_MULJN_4XCF32(x02, x00);
    y12 = BBE_MULJN_4XCF32(x02, x01);
    BBE_MULJAN_4XCF32(y00, x10, x10);
    BBE_MULJAN_4XCF32(y01, x11, x10);
    BBE_MULJAN_4XCF32(y11, x11, x11);
    BBE_MULJAN_4XCF32(y02, x12, x10);
    BBE_MULJAN_4XCF32(y12, x12, x11);
    BBE_MULJAN_4XCF32(y00, x20, x20);
    BBE_MULJAN_4XCF32(y01, x21, x20);
    BBE_MULJAN_4XCF32(y11, x21, x21);
    BBE_MULJAN_4XCF32(y02, x22, x20);
    BBE_MULJAN_4XCF32(y12, x22, x21);
    BBE_MULJAN_4XCF32(y00, x30, x30);
    BBE_MULJAN_4XCF32(y01, x31, x30);
    BBE_MULJAN_4XCF32(y11, x31, x31);
    BBE_MULJAN_4XCF32(y02, x32, x30);
    BBE_MULJAN_4XCF32(y12, x32, x31);

    y10 = BBE_CONJN_4XCF32(y01);
    y20 = BBE_CONJN_4XCF32(y02);
    y21 = BBE_CONJN_4XCF32(y12);

    /* Save results (half) */
    BBE_SVN_4XCF32_XP(y00, py,    L*sz_cf32);
    BBE_SVN_4XCF32_XP(y01, py,  4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y11, py, -3*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y02, py,  4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y12, py, -2*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y10, py,  4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y20, py,    L*sz_cf32);
    BBE_SVN_4XCF32_XP(y21, py, -9*L*sz_cf32+2*BBE_SIMD_WIDTH);
  }

  /* Compute next half of matrix */

  px  = (const xb_vecN_4xcf32 *)(x+3*L);
  py  = (      xb_vecN_4xcf32 *)(y+3*L);

  __Pragma("loop_count min=1");
  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
  {
    /* Load input matrix x */
    BBE_LVN_4XCF32_XP(x03, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x02, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x01, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x00, px,  7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x13, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x12, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x11, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x10, px,  7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x23, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x22, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x21, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x20, px,  7*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x33, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x32, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x31, px, -1*L*sz_cf32);
    BBE_LVN_4XCF32_XP(x30, px, -9*L*sz_cf32+2*BBE_SIMD_WIDTH);

    /* Compute half of matrix Y=X'*X */
    y03 = BBE_MULJN_4XCF32(x03, x00);
    y13 = BBE_MULJN_4XCF32(x03, x01);
    y23 = BBE_MULJN_4XCF32(x03, x02);
    y33 = BBE_MULJN_4XCF32(x03, x03);
    y22 = BBE_MULJN_4XCF32(x02, x02);
    BBE_MULJAN_4XCF32(y03, x13, x10);
    BBE_MULJAN_4XCF32(y13, x13, x11);
    BBE_MULJAN_4XCF32(y23, x13, x12);
    BBE_MULJAN_4XCF32(y33, x13, x13);
    BBE_MULJAN_4XCF32(y22, x12, x12);
    BBE_MULJAN_4XCF32(y03, x23, x20);
    BBE_MULJAN_4XCF32(y13, x23, x21);
    BBE_MULJAN_4XCF32(y23, x23, x22);
    BBE_MULJAN_4XCF32(y33, x23, x23);
    BBE_MULJAN_4XCF32(y22, x22, x22);
    BBE_MULJAN_4XCF32(y03, x33, x30);
    BBE_MULJAN_4XCF32(y13, x33, x31);
    BBE_MULJAN_4XCF32(y23, x33, x32);
    BBE_MULJAN_4XCF32(y33, x33, x33);
    BBE_MULJAN_4XCF32(y22, x32, x32);

    y30 = BBE_CONJN_4XCF32(y03);
    y31 = BBE_CONJN_4XCF32(y13);
    y32 = BBE_CONJN_4XCF32(y23);

    /* Save results (half) */
    BBE_SVN_4XCF32_XP(y03, py,   4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y13, py,   4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y23, py,     L*sz_cf32);
    BBE_SVN_4XCF32_XP(y30, py,     L*sz_cf32);
    BBE_SVN_4XCF32_XP(y31, py,     L*sz_cf32);
    BBE_SVN_4XCF32_XP(y32, py,  -4*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y22, py,   5*L*sz_cf32);
    BBE_SVN_4XCF32_XP(y33, py, -12*L*sz_cf32+2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatherm4x4sf() */
#endif

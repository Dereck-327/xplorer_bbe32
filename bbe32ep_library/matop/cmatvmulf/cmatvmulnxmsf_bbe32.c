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
DISCARD_FUN(void, cmatvmulnxmsf,( complex_float * restrict z, 
                            const complex_float * restrict x, 
                            const complex_float * restrict y, 
                            int N, int M, int L ))
#else

#include <string.h>
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

/* Streaming Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmulnxmsf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int N, int M, int L )
{
  const xb_vecN_4xcf32 * restrict inX;
  const xb_vecN_4xcf32 * restrict inX_;
  const xb_vecN_4xcf32 * restrict inY;
        xb_vecN_4xcf32 * restrict out;

  xb_vecN_4xcf32 regX0, regX1, regX2, regX3, regY, regY_;
  xb_vecN_4xcf32 regout0, regout1, regout2, regout3;
  int m, n;
  int strideX, strideZ, k, modinc;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  if (M<=0 || L<=0) return;
  if (N <= 0)
  {
    memset(z, 0, M*L*sizeof(complex_float));
    return;
  }

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  WUR_CBEGIN((uintptr_t)(y));
  WUR_CEND  ((uintptr_t)(y+L));
  inX = (const xb_vecN_4xcf32 *)(x);
  inY = (const xb_vecN_4xcf32 *)(y);
  out = (      xb_vecN_4xcf32 *)(z);

  /* Compute matrices by 4 rows per iteration */
  for (m = 0; m < (M>>2)*(L >> (LOG2_BBE_SIMD_WIDTH-2)); m++)
  {
    regout0 = BBE_ZERON_4XCF32();
    regout1 = BBE_ZERON_4XCF32();
    regout2 = BBE_ZERON_4XCF32();
    regout3 = BBE_ZERON_4XCF32();
    /* compute by 4 values for each 4 rows *
     * in the innermost loop               */
    __Pragma("loop_count min=1");
    for (n = 0; n < N; n++)
    {
      /* load input data from x for each 4 rows */
      BBE_LVN_4XCF32_XP(regX0, inX,    sz_cf32*L*N);
      BBE_LVN_4XCF32_XP(regX1, inX,    sz_cf32*L*N);
      BBE_LVN_4XCF32_XP(regX2, inX,    sz_cf32*L*N);
      BBE_LVN_4XCF32_XP(regX3, inX, -3*sz_cf32*L*N+sz_cf32*L);
      /* load input data from y for all 4 rows */
      regY_ = BBE_LVN_4XCF32_X(inY, 0);
      BBE_LVN_4XCF32_XP(regY, inY, sz_cf32*L);
      /* perform multiplication */
      BBE_MULAN_4XCF32(regout0, regX0, regY_);
      BBE_MULAN_4XCF32(regout1, regX1, regY_);
      BBE_MULAN_4XCF32(regout2, regX2, regY);
      BBE_MULAN_4XCF32(regout3, regX3, regY);
    }

    /* move pointers to the next matrices or next row */
    k = BBE_ADDMOD16U(k, modinc);
    strideX = -N*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(strideX, 3*N*L*sz_cf32-L*sz_cf32+2*BBE_SIMD_WIDTH, k);
    strideZ = -3*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(strideZ, 2*BBE_SIMD_WIDTH, k);

    inX = (const xb_vecN_4xcf32 *)((intptr_t)inX + strideX);
    inY = (const xb_vecN_4xcf32 *)XT_ADDX8(-N*L, (uintptr_t)inY);
    BBE_LVN_4XCF32_IC(regY, inY);/* make speculative load to move pointer using circular addressing */

    /* Save result */
    BBE_SVN_4XCF32_XP(regout0, out, sz_cf32*L);
    BBE_SVN_4XCF32_XP(regout1, out, sz_cf32*L);
    BBE_SVN_4XCF32_XP(regout2, out, sz_cf32*L);
    BBE_SVN_4XCF32_XP(regout3, out, strideZ);
  }

  /* compute last (M%4) rows for L matrices */
  if (M&2)
  {
    int _M = M & (~3);

    inX = (const xb_vecN_4xcf32 *)(x + _M*N*L);
    inX_= (const xb_vecN_4xcf32 *)((complex_float *)inX + L*N);
    out = (      xb_vecN_4xcf32 *)(z + _M*L);

    __Pragma("loop_count min=1");
    for (m = 0; m < (L >> (LOG2_BBE_SIMD_WIDTH-2)); m++)
    {
      regout0 = regout1 = regout2 = regout3 = BBE_ZERON_4XCF32();

      __Pragma("loop_count min=1");
      for (n = 0; n < N; n++)
      {
        BBE_LVN_4XCF32_XP(regX0, inX,  sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX1, inX_, sz_cf32*L);
        BBE_LVN_4XCF32_XP(regY , inY,  sz_cf32*L);
        BBE_MULMASN_4XCF32(regout0, regX0, regY, 0, 0x4);
        BBE_MULMASN_4XCF32(regout1, regX0, regY, 1, 0xB);
        BBE_MULMASN_4XCF32(regout2, regX1, regY, 0, 0x4);
        BBE_MULMASN_4XCF32(regout3, regX1, regY, 1, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout1);
      regout2 = BBE_ADDN_4XCF32(regout2, regout3);
      /* move pointers to the next matrices */
      strideX = -N*L*sz_cf32+2*BBE_SIMD_WIDTH;
      strideZ = -L*sz_cf32+2*BBE_SIMD_WIDTH;
      inX  = (const xb_vecN_4xcf32 *)((intptr_t)inX  + strideX);
      inX_ = (const xb_vecN_4xcf32 *)((intptr_t)inX_ + strideX);
      inY  = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4-N*L, (uintptr_t)inY);

      /* Save result */
      BBE_SVN_4XCF32_XP(regout0, out, sz_cf32*L);
      BBE_SVN_4XCF32_XP(regout2, out, strideZ);
    }
  }
  if (M&1)
  {
    int _M = M & (~1);

    inX = (const xb_vecN_4xcf32 *)(x + _M*N*L);
    inY = (const xb_vecN_4xcf32 *)(y);
    out = (      xb_vecN_4xcf32 *)(z + _M*L);

    __Pragma("loop_count min=1");
    for (m = 0; m < (L >> (LOG2_BBE_SIMD_WIDTH-2)); m++)
    {
      regout0 = regout1 = BBE_ZERON_4XCF32();

      __Pragma("loop_count min=1");
      for (n = 0; n < N; n++)
      {
        BBE_LVN_4XCF32_XP(regY , inY, sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX0, inX, sz_cf32*L);
        BBE_MULMASN_4XCF32(regout0, regX0, regY, 0, 0x4);
        BBE_MULMASN_4XCF32(regout1, regX0, regY, 1, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout1);
      /* Save result */
      BBE_SVN_4XCF32_IP(regout0, out, 2*BBE_SIMD_WIDTH);

      strideX = -N*L*sz_cf32+2*BBE_SIMD_WIDTH;
      inX = (const xb_vecN_4xcf32 *)((intptr_t)inX + strideX);
      inY = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4-N*L, (uintptr_t)inY);
    }
  }

} /* cmatvmulnxmsf() */
#endif

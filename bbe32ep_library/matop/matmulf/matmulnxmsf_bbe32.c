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
DISCARD_FUN(void, matmulnxmsf,( float32_t * restrict z,
                          const float32_t * restrict x,
                          const float32_t * restrict y,
                          int N, int M, int L ))
#else

#include <string.h>
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

/* Streaming Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 8
*/
void matmulnxmsf ( float32_t * restrict z,
             const float32_t * restrict x,
             const float32_t * restrict y,
             int N, int M, int L )
{
  const xb_vecN_2xf32 * restrict inX;
  const xb_vecN_2xf32 * restrict inXmod;
  const xb_vecN_2xf32 * restrict inY;
        xb_vecN_2xf32 * restrict out = (xb_vecN_2xf32 *)z;

  xb_vecN_2xf32 regX0, regX1, regX2, regX3, regY;
  xb_vecN_2xf32 regout0, regout1, regout2, regout3;
  int m, n, p;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  if (M<=0 || L<=0) return;
  if (N <= 0)
  {
    memset(z, 0, M*M*L*sizeof(float32_t));
    return;
  }

  /* Compute matrices by 4 rows per iteration */
  for (m = 0; m < (M>>2); m++)
  {
    inY    = (const xb_vecN_2xf32 *)(y);
    inXmod = (const xb_vecN_2xf32 *)(x + m*4*N*L);
    out    = (      xb_vecN_2xf32 *)(z + m*4*M*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((float32_t *)inXmod + L));
    /* compute result by BBE_SIMD_WIDTH/2 elements for 4 rows and for L matrices */
    /* combine traverse of elements by M and L in one loop */
    __Pragma("loop_count min=1");
    for (p = 0; p < M*(L >> (LOG2_BBE_SIMD_WIDTH-1)); p++)
    {
      regout0 = BBE_ZERON_2XF32();
      regout1 = BBE_ZERON_2XF32();
      regout2 = BBE_ZERON_2XF32();
      regout3 = BBE_ZERON_2XF32();
      inX = inXmod;
      BBE_LVN_2XF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      /* compute by 4 values for each 4 rows *
       * in the innermost loop               */
      __Pragma("loop_count min=1");
      for (n = 0; n < N; n++)
      {
        /* load input data from x for each 4 rows */
        BBE_LVN_2XF32_XP(regX0, inX,    sz_f32*L*N);
        BBE_LVN_2XF32_XP(regX1, inX,    sz_f32*L*N);
        BBE_LVN_2XF32_XP(regX2, inX,    sz_f32*L*N);
        BBE_LVN_2XF32_XP(regX3, inX, -3*sz_f32*L*N+sz_f32*L);
        /* load input data from y for all 4 rows */
        BBE_LVN_2XF32_XP(regY, inY, sz_f32*L*M);
        /* perform multiplication */
        BBE_MULAN_2XF32(regout0, regX0, regY);
        BBE_MULAN_2XF32(regout1, regX1, regY);
        BBE_MULAN_2XF32(regout2, regX2, regY);
        BBE_MULAN_2XF32(regout3, regX3, regY);
      }
      BBE_SVN_2XF32_XP(regout0, out,    sz_f32*L*M);
      BBE_SVN_2XF32_XP(regout1, out,    sz_f32*L*M);
      BBE_SVN_2XF32_XP(regout2, out,    sz_f32*L*M);
      BBE_SVN_2XF32_XP(regout3, out, -3*sz_f32*L*M+2*BBE_SIMD_WIDTH);
      /* move pointer to the next matrices or next value in a row */
      inY = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH/2 - N*M*L, (uintptr_t)inY);
    }
  }

  /* compute last (M%4) rows for L matrices */
  if (M&2)
  {
    int _M = M & (~3);

    inY    = (const xb_vecN_2xf32 *)(y);
    inXmod = (const xb_vecN_2xf32 *)(x + _M*N*L);
    out    = (      xb_vecN_2xf32 *)(z + _M*M*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((float32_t *)inXmod + L));
    
    __Pragma("loop_count min=1");
    for (p = 0; p < M*(L >> (LOG2_BBE_SIMD_WIDTH-1)); p++)
    {
      regout0 = regout1 = BBE_ZERON_2XF32();
      inX = inXmod;
      BBE_LVN_2XF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      __Pragma("loop_count min=1");
      for (n = 0; n < N; n++)
      {
        BBE_LVN_2XF32_XP(regX0, inX,  sz_f32*L*N);
        BBE_LVN_2XF32_XP(regX1, inX, -sz_f32*L*N+sz_f32*L);
        BBE_LVN_2XF32_XP(regY , inY, sz_f32*L*M);
        BBE_MULAN_2XF32(regout0, regX0, regY);
        BBE_MULAN_2XF32(regout1, regX1, regY);
      }
      BBE_SVN_2XF32_XP(regout0, out,  sz_f32*L*M);
      BBE_SVN_2XF32_XP(regout1, out, -sz_f32*L*M+2*BBE_SIMD_WIDTH);
      inY = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH/2 - N*M*L, (uintptr_t)inY);
    }
  }
  if (M&1)
  {
    int _M = M & (~1);

    inY    = (const xb_vecN_2xf32 *)(y);
    inXmod = (const xb_vecN_2xf32 *)(x + _M*N*L);
    out    = (      xb_vecN_2xf32 *)(z + _M*M*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((float32_t *)inXmod + L));
    
    __Pragma("loop_count min=1");
    for (p = 0; p < M*(L >> (LOG2_BBE_SIMD_WIDTH-1)); p++)
    {
      regout0 = BBE_ZERON_2XF32();
      inX = inXmod;
      BBE_LVN_2XF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      __Pragma("loop_count min=1");
      for (n = 0; n < N; n++)
      {
        BBE_LVN_2XF32_XP(regY , inY, sz_f32*L*M);
        BBE_LVN_2XF32_XP(regX0, inX, sz_f32*L);
        BBE_MULAN_2XF32(regout0, regX0, regY);
      }
      BBE_SVN_2XF32_IP(regout0, out, 2*BBE_SIMD_WIDTH);
      inY = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH/2 - N*M*L, (uintptr_t)inY);
    }
  }

} /* matmulnxmsf() */
#endif

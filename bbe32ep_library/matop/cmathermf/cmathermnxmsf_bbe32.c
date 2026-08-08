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
DISCARD_FUN(void, cmathermnxmsf,( complex_float * restrict y, 
                            const complex_float * restrict x, 
                            int N, int M, int L ))
#else

#include <string.h>
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

/* Streaming Order, Floating-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     L must be a multiple of 4
*/
void cmathermnxmsf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int N, int M, int L )
{
#if 0
  const xb_vecN_4xcf32 * restrict inX;
  const xb_vecN_4xcf32 * restrict inXmod;
  const xb_vecN_4xcf32 * restrict inXc;
        xb_vecN_4xcf32 * restrict out;

  xb_vecN_4xcf32 regX0, regX1, regX2, regX3, regXc_, regXc;
  xb_vecN_4xcf32 regout0, regout1, regout2, regout3;
  int m, n, p;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  if (N<=0 || L<=0) return;
  if (M <= 0)
  {
    memset(y, 0, N*N*L*sizeof(complex_float));
    return;
  }

  /* Compute matrices by 4 rows per iteration */
  for (n = 0; n < (N>>2); n++)
  {
    inXc   = (const xb_vecN_4xcf32 *)(x);
    inXmod = (const xb_vecN_4xcf32 *)(x + n*4*L);
    out    = (      xb_vecN_4xcf32 *)(y + n*4*N*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((complex_float *)inXmod + L));
    /* compute result by BBE_SIMD_WIDTH/4 elements for 4 rows and for L matrices */
    /* combine traverse of elements by N and L in one loop */
    __Pragma("loop_count min=1");
    for (p = 0; p < N*(L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = BBE_ZERON_4XCF32();
      regout1 = BBE_ZERON_4XCF32();
      regout2 = BBE_ZERON_4XCF32();
      regout3 = BBE_ZERON_4XCF32();
      inX = inXmod;
      BBE_LVN_4XCF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      /* compute by 4 values for each 4 rows *
       * in the innermost loop               */
      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        /* load input data from x for each 4 rows */
        BBE_LVN_4XCF32_XP(regX0, inX,    sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX1, inX,    sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX2, inX,    sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX3, inX, -3*sz_cf32*L+sz_cf32*L*N);
        /* load input data from y for all 4 rows */
        regXc_ = BBE_LVN_4XCF32_X(inXc, 0);
        BBE_LVN_4XCF32_XP(regXc, inXc, sz_cf32*L*N);
        /* perform multiplication */
        BBE_MULJAN_4XCF32(regout0, regXc_, regX0);
        BBE_MULJAN_4XCF32(regout1, regXc_, regX1);
        BBE_MULJAN_4XCF32(regout2, regXc , regX2);
        BBE_MULJAN_4XCF32(regout3, regXc , regX3);
      }
      BBE_SVN_4XCF32_XP(regout0, out,    sz_cf32*L*N);
      BBE_SVN_4XCF32_XP(regout1, out,    sz_cf32*L*N);
      BBE_SVN_4XCF32_XP(regout2, out,    sz_cf32*L*N);
      BBE_SVN_4XCF32_XP(regout3, out, -3*sz_cf32*L*N+2*BBE_SIMD_WIDTH);
      /* move pointer to the next matrices or next value in a row */
      inXc = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4 - N*M*L, (uintptr_t)inXc);
    }
  }

  /* compute last (N%4) rows for L matrices */
  if (N&2)
  {
    int _N = N & (~3);

    inXc   = (const xb_vecN_4xcf32 *)(x);
    inXmod = (const xb_vecN_4xcf32 *)(x + _N*L);
    out    = (      xb_vecN_4xcf32 *)(y + _N*N*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((complex_float *)inXmod + L));
    
    __Pragma("loop_count min=1");
    for (p = 0; p < N*(L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = regout1 = regout2 = regout3 = BBE_ZERON_4XCF32();
      inX = inXmod;
      BBE_LVN_4XCF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        BBE_LVN_4XCF32_XP(regX0, inX,  sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX1, inX, -sz_cf32*L+sz_cf32*L*N);
        BBE_LVN_4XCF32_XP(regXc, inXc, sz_cf32*L*N);
        BBE_MULMASN_4XCF32(regout0, regX0, regXc, 0, 0x4);
        BBE_MULMASN_4XCF32(regout1, regX0, regXc, 2, 0xB);
        BBE_MULMASN_4XCF32(regout2, regX1, regXc, 0, 0x4);
        BBE_MULMASN_4XCF32(regout3, regX1, regXc, 2, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout1);
      regout2 = BBE_ADDN_4XCF32(regout2, regout3);
      BBE_SVN_4XCF32_XP(regout0, out,  sz_cf32*L*N);
      BBE_SVN_4XCF32_XP(regout2, out, -sz_cf32*L*N+2*BBE_SIMD_WIDTH);
      inXc = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4 - N*M*L, (uintptr_t)inXc);
    }
  }
  if (N&1)
  {
    int _N = N & (~1);

    inXc   = (const xb_vecN_4xcf32 *)(x);
    inXmod = (const xb_vecN_4xcf32 *)(x + _N*L);
    out    = (      xb_vecN_4xcf32 *)(y + _N*N*L);
    WUR_CBEGIN((uintptr_t)inXmod);
    WUR_CEND((uintptr_t)((complex_float *)inXmod + L));
    
    __Pragma("loop_count min=1");
    for (p = 0; p < N*(L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = regout1 = BBE_ZERON_4XCF32();
      inX = inXmod;
      BBE_LVN_4XCF32_IC(regX0, inXmod);/* make speculative load to move pointer using circular addressing */
      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        BBE_LVN_4XCF32_XP(regXc, inXc, sz_cf32*L*N);
        BBE_LVN_4XCF32_XP(regX0, inX , sz_cf32*L*N);
        BBE_MULMASN_4XCF32(regout0, regX0, regXc, 0, 0x4);
        BBE_MULMASN_4XCF32(regout1, regX0, regXc, 2, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout1);
      BBE_SVN_4XCF32_IP(regout0, out, 2*BBE_SIMD_WIDTH);
      inXc = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4 - N*M*L, (uintptr_t)inXc);
    }
  }
#else
  const xb_vecN_4xcf32 * restrict inX;
  const xb_vecN_4xcf32 * restrict inXmod;
  const xb_vecN_4xcf32 * restrict inXc;
        xb_vecN_4xcf32 * restrict out_row;
        xb_vecN_4xcf32 * restrict out_col;

  xb_vecN_4xcf32 regX0, regX1, regXc_, regXc;
  xb_vecN_4xcf32 regout0, regout1, regout2, regout3;
  int stridex, stridey;
  int k, modinc;
  int m, n, p;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  if (N<=0 || L<=0) return;
  if (M <= 0)
  {
    memset(y, 0, N*N*L*sizeof(complex_float));
    return;
  }

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  /* Main part:                                              */
  /* compute matrices by 2 rows and 2 columns per iteration. */
  for (n = 0; n < (N>>1); n++)
  {

    /*
     * First compute 2x2 block including 2 diagonal elements of matrix.
     */

    inXmod  = (const xb_vecN_4xcf32 *)(x + n*2*L);
    out_row = (      xb_vecN_4xcf32 *)(y + n*2*N*L + n*2*L);
    out_col = (      xb_vecN_4xcf32 *)((complex_float *)out_row + N*L);

    __Pragma("loop_count min=1");
    for (p = 0; p < (L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = BBE_ZERON_4XCF32();
      regout1 = BBE_ZERON_4XCF32();
      regout2 = BBE_ZERON_4XCF32();
      regout3 = BBE_ZERON_4XCF32();
      inX = inXmod;
      inXmod++;/* jump to next BBE_SIMD_WIDTH/4 matrices */

      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        /* load input data from X for 2x2 block */
        BBE_LVN_4XCF32_XP(regX0, inX,  sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX1, inX, -sz_cf32*L+sz_cf32*L*N);
        /* perform multiplication */
        BBE_MULJAN_4XCF32(regout0, regX0, regX0);
        BBE_MULJAN_4XCF32(regout2, regX0, regX1);
        BBE_MULJAN_4XCF32(regout3, regX1, regX1);
      }
      regout1 = BBE_CONJN_4XCF32(regout2);

      BBE_SVN_4XCF32_XP(regout0, out_row,  sz_cf32*L);
      BBE_SVN_4XCF32_XP(regout1, out_row, -sz_cf32*L+2*BBE_SIMD_WIDTH);
      BBE_SVN_4XCF32_XP(regout2, out_col,  sz_cf32*L);
      BBE_SVN_4XCF32_XP(regout3, out_col, -sz_cf32*L+2*BBE_SIMD_WIDTH);
    }

    /*
     * Second step: compute matrices by 2 rows and 2 columns per iteration,
     * columns are conjuncted with rows.
     */
    
    inXmod  = (const xb_vecN_4xcf32 *)(x + n*2*L);
    inXc    = (const xb_vecN_4xcf32 *)(x + n*2*L + 2*L);
    out_row = (      xb_vecN_4xcf32 *)((complex_float *)out_row + L);
    out_col = (      xb_vecN_4xcf32 *)((complex_float *)out_col + N*L - L);

    for (p = 0; p < (N-2-n*2)*(L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = BBE_ZERON_4XCF32();
      regout1 = BBE_ZERON_4XCF32();
      regout2 = BBE_ZERON_4XCF32();
      regout3 = BBE_ZERON_4XCF32();
      inX = inXmod;

      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        /* load input left-hand data from X */
        BBE_LVN_4XCF32_XP(regX0, inX,  sz_cf32*L);
        BBE_LVN_4XCF32_XP(regX1, inX, -sz_cf32*L+sz_cf32*L*N);
        /* load input right-hand data from X */
        regXc_ = BBE_LVN_4XCF32_X(inXc, 0);
        BBE_LVN_4XCF32_XP(regXc , inXc, sz_cf32*L*N);
        /* perform multiplication */
        BBE_MULMASN_4XCF32(regout0, regX0, regXc_, 0, 0x4);
        BBE_MULMASN_4XCF32(regout2, regX0, regXc_, 2, 0xB);
        BBE_MULMASN_4XCF32(regout1, regX1, regXc , 0, 0x4);
        BBE_MULMASN_4XCF32(regout3, regX1, regXc , 2, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout2);
      regout1 = BBE_ADDN_4XCF32(regout1, regout3);
      regout2 = BBE_CONJN_4XCF32(regout0);
      regout3 = BBE_CONJN_4XCF32(regout1);

      /* Save results */
      BBE_SVN_4XCF32_XP(regout0, out_row,  sz_cf32*L*N);
      BBE_SVN_4XCF32_XP(regout1, out_row, -sz_cf32*L*N+2*BBE_SIMD_WIDTH);
      BBE_SVN_4XCF32_XP(regout2, out_col,  sz_cf32*L);
      BBE_SVN_4XCF32_XP(regout3, out_col, -sz_cf32*L);

      /* move pointers to the next BBE_SIMD_WIDTH/4 matrices or next values of current matrices */
      k = BBE_ADDMOD16U(k, modinc);
      stridex = BBE_SIMD_WIDTH/4;
      XT_MOVEQZ(stridex, -L+BBE_SIMD_WIDTH/4, k);
      stridey = BBE_SIMD_WIDTH/4;
      XT_MOVEQZ(stridey, -L+BBE_SIMD_WIDTH/4+L*N, k);

      inXmod  = (const xb_vecN_4xcf32 *)XT_ADDX8(stridex, (uintptr_t)inXmod);
      inXc    = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4 - N*M*L, (uintptr_t)inXc);
      out_col = (      xb_vecN_4xcf32 *)XT_ADDX8(stridey, (uintptr_t)out_col);
    }
  }

  /* compute last (N*N & 1) value for L matrices */
  if (N&1)
  {
    int N_ = N & (~1);

    inXc    = (const xb_vecN_4xcf32 *)(x + N_*L);
    out_row = (      xb_vecN_4xcf32 *)(y + N_*L + N_*N*L);
    
    __Pragma("loop_count min=1");
    for (p = 0; p < (L >> (LOG2_BBE_SIMD_WIDTH-2)); p++)
    {
      regout0 = regout1 = BBE_ZERON_4XCF32();
      __Pragma("loop_count min=1");
      for (m = 0; m < M; m++)
      {
        BBE_LVN_4XCF32_XP(regXc, inXc, sz_cf32*L*N);
        BBE_MULMASN_4XCF32(regout0, regXc, regXc, 0, 0x4);
        BBE_MULMASN_4XCF32(regout1, regXc, regXc, 2, 0xB);
      }
      regout0 = BBE_ADDN_4XCF32(regout0, regout1);
      BBE_SVN_4XCF32_IP(regout0, out_row, 2*BBE_SIMD_WIDTH);

      inXc = (const xb_vecN_4xcf32 *)XT_ADDX8(BBE_SIMD_WIDTH/4 - N*M*L, (uintptr_t)inXc);
    }
  }
#endif
} /* cmathermnxmsf() */
#endif

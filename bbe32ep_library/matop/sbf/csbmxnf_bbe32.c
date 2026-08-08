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
    Streaming to Packed Conversion for Real and Complex Matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#define sz_cf32 (int)sizeof(complex_float)

/*-------------------------------------------------------------------------
Streaming to Packed Conversion for Real and Complex Matrices

Description: convert a sequence of L MxN matrices from streaming to packed
(block) order. Use rsb*() functions for real data, and csb*() functions - for
complex data.

Representation:
<r|c>sb<size>    16-bit fixed-point data
<r|c>sb<size>f   IEEE-754 Std single precision floating-point data

Storage size SY denotes the number of data elements required to store an
MxN matrix Y in block order. If matrix size M*N is less than the SIMD vector
size for appropriate data type, then the storage size equals M*N rounded up
to the next power of two, otherwise SY equals M*N rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for real fixed-point data 2*BBE_SIMD_WIDTH/sizeof(int16_t) == 16
  - for complex fixed-point data 2*BBE_SIMD_WIDTH/sizeof(complex_fract16) == 8
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4

Parameters:
Input:
x[M*N][L]  Input sequence of MxN matrices, streaming order
M          Number of rows in a matrix
N          Number of columns in a matrix
L          Number of matrices
Output:
y[L][Sy]   Output sequence of matrices, block order.

Restrictions:
x,y        Aligned on 32-byte boundary
x,y        Must not overlap
L          Must be a multiple of:
             16 for real fixed-point data 
              8 for complex fixed-point data and real floating-point data
              4 for complex floating-point data
-------------------------------------------------------------------------*/

/* Sx=M*N, Sy=<see the description> */
void csbmxnf ( complex_float * restrict y, const complex_float * restrict x, int M, int N, int L )
{
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
  int l, k, incmod;
  xb_vecNx16 X00, X01, X02, X03,
                 X10, X11, X12, X13;
  xb_vecNx16 Y00, Y01, Y02, Y03,
                 Y10, Y11, Y12, Y13;
  int numel;
  int stridex, stridex_nextL, stridex_nextNM;
  int stridey, stridey_nextL, stridey_nextNM;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);

  /* Process special cases */
  if (N<=0 || M<=0 || L<=0) return;
  if (N*M == 1)
  {
      px0 = (const xb_vecNx16 *)(x);
      py0 = (      xb_vecNx16 *)(y);
      __Pragma("loop_count min=1");
      for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
      {
          BBE_LVNX16_IP(X00, px0, 2*BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(X00, py0, 2*BBE_SIMD_WIDTH);
      }
      return;
  }
  else if (N*M == 2)
  {
      csb2x1f(y, x, L);
      return;
  }

  px0 = (const xb_vecNx16 *)(x);
  px1 = (const xb_vecNx16 *)(x+L);
  py0 = (      xb_vecNx16 *)(y);
  py1 = (      xb_vecNx16 *)(y+4);
  k = 0;
  incmod = (L<<16) | (BBE_SIMD_WIDTH/4);

  numel = (N*M + BBE_SIMD_WIDTH/4-1) & ~(BBE_SIMD_WIDTH/4-1);

  stridex_nextL  = -6*L*sz_cf32+2*BBE_SIMD_WIDTH;
  stridex_nextNM = L*sz_cf32+2*BBE_SIMD_WIDTH;
  stridey_nextL  = numel*sz_cf32;
  stridey_nextNM = -(L-1)*numel*sz_cf32+2*2*BBE_SIMD_WIDTH;

  /* Convert by 8 values */
  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*((N*M)/8); l++ )
  {
    /* k = (k+(BBE_SIMD_WIDTH/4))%L */
    k = BBE_ADDMOD16U(k, incmod);

    stridex = stridex_nextL;
    XT_MOVEQZ(stridex, stridex_nextNM, k);

    stridey = stridey_nextL;
    XT_MOVEQZ(stridey, stridey_nextNM, k);

    BBE_LVNX16_XP(X00, px0, 2*L*sz_cf32);
    BBE_LVNX16_XP(X01, px1, 2*L*sz_cf32);
    BBE_LVNX16_XP(X02, px0, 2*L*sz_cf32);
    BBE_LVNX16_XP(X03, px1, 2*L*sz_cf32);

    BBE_LVNX16_XP(X10, px0, 2*L*sz_cf32);
    BBE_LVNX16_XP(X11, px1, 2*L*sz_cf32);
    BBE_LVNX16_XP(X12, px0, stridex);
    BBE_LVNX16_XP(X13, px1, stridex);

    BBE_DSELNX16I(X02, X00, X02, X00, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X03, X01, X03, X01, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y01, Y00, X01, X00, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y03, Y02, X03, X02, BBE_DSELI_INTERLEAVE_4);

    BBE_DSELNX16I(X12, X10, X12, X10, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X13, X11, X13, X11, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y11, Y10, X11, X10, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y13, Y12, X13, X12, BBE_DSELI_INTERLEAVE_4);

    BBE_SVNX16_XP(Y00, py0, stridey_nextL);
    BBE_SVNX16_XP(Y01, py0, stridey_nextL);
    BBE_SVNX16_XP(Y02, py0, stridey_nextL);
    BBE_SVNX16_XP(Y03, py0, stridey);

    BBE_SVNX16_XP(Y10, py1, stridey_nextL);
    BBE_SVNX16_XP(Y11, py1, stridey_nextL);
    BBE_SVNX16_XP(Y12, py1, stridey_nextL);
    BBE_SVNX16_XP(Y13, py1, stridey);
  }
  /* Process last (M*N)%8 values for L matrices */
  if (N*M & 4)
  {
    int ix1, ix2, ix3;
    ix1 = 1*L*sz_cf32;
    ix2 = 2*L*sz_cf32;
    ix3 = 3*L*sz_cf32;

    __Pragma("loop_count min=1");
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
      X03 = BBE_LVNX16_X(px0, ix3);
      X02 = BBE_LVNX16_X(px0, ix2);
      X01 = BBE_LVNX16_X(px0, ix1);
      BBE_LVNX16_IP(X00, px0, 2*BBE_SIMD_WIDTH);

      BBE_DSELNX16I(X02, X00, X02, X00, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(X03, X01, X03, X01, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y01, Y00, X01, X00, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y03, Y02, X03, X02, BBE_DSELI_INTERLEAVE_4);

      BBE_SVNX16_XP(Y00, py0, stridey_nextL);
      BBE_SVNX16_XP(Y01, py0, stridey_nextL);
      BBE_SVNX16_XP(Y02, py0, stridey_nextL);
      BBE_SVNX16_XP(Y03, py0, stridey_nextL);
    }

    px0 = (const xb_vecNx16 *)((complex_float *)px0+3*L);
  }
  if (N*M & 3)
  {
    int NM_mod, ix1, ix2;
    NM_mod = (N*M) & 3;
    ix1 = XT_MIN(NM_mod-1, 1)*L*sz_cf32;
    ix2 = (NM_mod-1)*L*sz_cf32;

    py0 = (xb_vecNx16 *)(y + ((N*M) & ~3));
    X13 = BBE_ZERONX16();

    __Pragma("loop_count min=1");
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
      X02 = BBE_LVNX16_X(px0, ix2);
      X01 = BBE_LVNX16_X(px0, ix1);
      BBE_LVNX16_IP(X00, px0, 2*BBE_SIMD_WIDTH);

      BBE_DSELNX16I(X02, X00, X02, X00, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(X03, X01, X13, X01, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y01, Y00, X01, X00, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y03, Y02, X03, X02, BBE_DSELI_INTERLEAVE_4);

      BBE_SVNX16_XP(Y00, py0, stridey_nextL);
      BBE_SVNX16_XP(Y01, py0, stridey_nextL);
      BBE_SVNX16_XP(Y02, py0, stridey_nextL);
      BBE_SVNX16_XP(Y03, py0, stridey_nextL);
    }
  }
} /* csbmxnf() */

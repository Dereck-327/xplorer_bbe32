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

/* M=4, N=4, Sx=16, Sy=16 */
void csb4x4f ( complex_float * restrict y, const complex_float * restrict x, int L )
{
#if 0
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
  const xb_vecNx16 * restrict px2;
  const xb_vecNx16 * restrict px3;
        xb_vecNx16 * restrict py;
  xb_vecNx16 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  int l, k;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  px3 = (const xb_vecNx16 *)(x);
  py  = (      xb_vecNx16 *)(y);

  /* Convert by 4 values */
  for ( k=0; k<4; k++ )
  {
    px0 = (const xb_vecNx16 *)(px3);
    px1 = (const xb_vecNx16 *)((complex_float *)px0+L);
    px2 = (const xb_vecNx16 *)((complex_float *)px1+L);
    px3 = (const xb_vecNx16 *)((complex_float *)px2+L);
    py  = (      xb_vecNx16 *)(y);

    __Pragma("loop_count min=1");
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1, px1, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X2, px2, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X3, px3, 2*BBE_SIMD_WIDTH);

      BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

      BBE_SVNX16_IP(Y0, py, 4*2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y1, py, 4*2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y2, py, 4*2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y3, py, 4*2*BBE_SIMD_WIDTH);
    }
    y += 4;
  }
#else
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
  int l, k, incmod;
  xb_vecNx16 X00, X01, X02, X03,
                 X10, X11, X12, X13;
  xb_vecNx16 Y00, Y01, Y02, Y03,
                 Y10, Y11, Y12, Y13;
  int stridex, stridey;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  px0 = (const xb_vecNx16 *)(x);
  px1 = (const xb_vecNx16 *)(x+L);
  py0 = (      xb_vecNx16 *)(y);
  py1 = (      xb_vecNx16 *)(y+4);
  k = 0;
  incmod = (L<<16) | (BBE_SIMD_WIDTH/4);

  /* Convert by 8 values */
  __Pragma("loop_count min=2, factor=2");
  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*2; l++ )
  {
    k = BBE_ADDMOD16U(k, incmod);
    stridex = -6*L*sz_cf32+2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridex, L*sz_cf32+2*BBE_SIMD_WIDTH, k);
    stridey = 4*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(stridey, -L*16*sz_cf32+6*2*BBE_SIMD_WIDTH, k);

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

    BBE_SVNX16_IP(Y00, py0, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y01, py0, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y02, py0, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(Y03, py0, stridey);

    BBE_SVNX16_IP(Y10, py1, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y11, py1, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y12, py1, 4*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_XP(Y13, py1, stridey);
  }
#endif
} /* csb4x4f() */

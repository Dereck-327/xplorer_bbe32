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
#include "bs_common.h"
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

/* M=8, N=8, Sx=64, Sy=64 */
void rsb8x8 ( int16_t * restrict y, const int16_t * restrict x, int L )
{
  const xb_vecNx16* restrict w = (const xb_vecNx16*)x;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  int k, l, _8L, _7L, back, MN;
  int wstride, zstride;
  int wstridebig, zstridebig;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((2 * L) % BBE_SIMD_WIDTH == 0);
  if (L<=0) return;
  MN = 64;
  _8L = (MN << 4);
  _7L = (MN << 1) - _8L;
  back = 2 * BBE_SIMD_WIDTH - _8L + _7L;
  wstridebig = 2 * BBE_SIMD_WIDTH - (MN - 1)*(L << 1);
  zstridebig = back + 2 * BBE_SIMD_WIDTH*(MN - (MN >> 4));

  __Pragma("loop_count min=1");

  for (k = l = 0; l<(MN >> 4)*(L >> 4); l++)
  {
    k = BBE_ADDMOD16U(k, ((MN >> 4) << 16) | 1);
    wstride = (L << 1);
    zstride = back;
    XT_MOVEQZ(wstride, wstridebig, k);
    XT_MOVEQZ(zstride, zstridebig, k);

    BBE_LVNX16_XP(Y0, w, (L << 1));
    BBE_LVNX16_XP(Y1, w, (L << 1));
    BBE_LVNX16_XP(Y2, w, (L << 1));
    BBE_LVNX16_XP(Y3, w, (L << 1));
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y4, w, (L << 1));
    BBE_LVNX16_XP(Y5, w, (L << 1));
    BBE_LVNX16_XP(Y6, w, (L << 1));
    BBE_LVNX16_XP(Y7, w, (L << 1));
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y8, w, (L << 1));
    BBE_LVNX16_XP(Y9, w, (L << 1));
    BBE_LVNX16_XP(YA, w, (L << 1));
    BBE_LVNX16_XP(YB, w, (L << 1));
    BBE_DSELNX16I(X9, X8, Y9, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XB, XA, YB, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YA, Y8, XA, X8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YB, Y9, XB, X9, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(YC, w, (L << 1));
    BBE_LVNX16_XP(YD, w, (L << 1));
    BBE_LVNX16_XP(YE, w, (L << 1));
    BBE_LVNX16_XP(YF, w, wstride);
    BBE_DSELNX16I(XD, XC, YD, YC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XE, YF, YE, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YE, YC, XE, XC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YF, YD, XF, XD, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(X4, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X1, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X6, X2, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X3, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(XC, X8, YC, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XD, X9, YD, Y9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XE, XA, YE, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XB, YF, YB, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(Y8, Y0, X8, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y0, z, _8L);
    BBE_SVNX16_XP(Y8, z, _7L);
    BBE_DSELNX16I(Y9, Y1, X9, X1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y1, z, _8L);
    BBE_SVNX16_XP(Y9, z, _7L);
    BBE_DSELNX16I(YA, Y2, XA, X2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y2, z, _8L);
    BBE_SVNX16_XP(YA, z, _7L);
    BBE_DSELNX16I(YB, Y3, XB, X3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y3, z, _8L);
    BBE_SVNX16_XP(YB, z, _7L);
    BBE_DSELNX16I(YC, Y4, XC, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y4, z, _8L);
    BBE_SVNX16_XP(YC, z, _7L);
    BBE_DSELNX16I(YD, Y5, XD, X5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y5, z, _8L);
    BBE_SVNX16_XP(YD, z, _7L);
    BBE_DSELNX16I(YE, Y6, XE, X6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y6, z, _8L);
    BBE_SVNX16_XP(YE, z, _7L);
    BBE_DSELNX16I(YF, Y7, XF, X7, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y7, z, _8L);
    BBE_SVNX16_XP(YF, z, zstride);
  }
} /* rsb8x8() */

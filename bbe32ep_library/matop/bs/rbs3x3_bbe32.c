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
    Packed to Streaming Conversion for Real and Complex Matrices
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
Packed to Streaming Conversion for Real and Complex Matrices

Description: convert a sequence of L MxN matrices from packed (block) order
to streaming order. Use rbs*() functions for real data, and cbs*() functions -
for complex data.

Representation:
<r|c>bs<size>    16-bit fixed-point data
<r|c>bs<size>f   IEEE-754 Std single precision floating-point data

Storage size SX denotes the number of data elements required to store an
MxN matrix X in block order. If matrix size M*N is less than the SIMD vector
size for appropriate data type, then the storage size equals M*N rounded up
to the next power of two, otherwise SX equals M*N rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for real fixed-point data 2*BBE_SIMD_WIDTH/sizeof(int16_t) == 16
  - for complex fixed-point data 2*BBE_SIMD_WIDTH/sizeof(complex_fract16) == 8
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4

Parameters:
Input:
x[L][SX]   Input sequence of MxN matrices, block order
M          Number of rows in a matrix
N          Number of columns in a matrix
L          Number of matrices
Output:
y[M*N][L]  Output sequence of matrices, streaming order

Restrictions:
x,y      Aligned on 32-byte boundary
x,y      Must not overlap
L        Must be a multiple of:
           16 for real fixed-point data 
            8 for complex fixed-point data and real floating-point data
            4 for complex floating-point data
-------------------------------------------------------------------------*/

/* M=3, N=3, Sx=16, Sy=9 */
void rbs3x3 ( int16_t * restrict y, const int16_t * restrict x, int L )
{
  const xb_vecNx16* restrict X = (const xb_vecNx16*)x;
  xb_vecNx16* restrict Y = (xb_vecNx16*)y;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
  int k, back, _L;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  _L = (L << 1);
  back = 2 * BBE_SIMD_WIDTH - 8 * _L;
  if (L <= 0) return;
  for (k = 0; k<(L >> LOG2_BBE_SIMD_WIDTH); k++)
  {
    BBE_LVNX16_IP(Y0, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y2, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y3, X, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(Y4, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y5, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y6, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y7, X, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(Y8, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y9, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YA, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YB, X, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X9, X8, Y9, Y8, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XB, XA, YB, YA, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YA, Y8, XA, X8, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YB, Y9, XB, X9, BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(YC, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YD, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YE, X, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YF, X, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(XD, XC, YD, YC, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XF, XE, YF, YE, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YE, YC, XE, XC, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YF, YD, XF, XD, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(X4, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X5, X1, Y5, Y1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X6, X2, Y6, Y2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7, X3, Y7, Y3, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(XC, X8, YC, Y8, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XD, X9, YD, Y9, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XE, XA, YE, YA, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XF, XB, YF, YB, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(Y8, Y0, X8, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y9, Y1, X9, X1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YA, Y2, XA, X2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YB, Y3, XB, X3, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YC, Y4, XC, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YD, Y5, XD, X5, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YE, Y6, XE, X6, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YF, Y7, XF, X7, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(Y1, Y0, X8, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y0, Y, _L);
    BBE_SVNX16_XP(Y1, Y, _L);
    BBE_DSELNX16I(Y3, Y2, X9, X1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y2, Y, _L);
    BBE_SVNX16_XP(Y3, Y, _L);
    BBE_DSELNX16I(Y5, Y4, XA, X2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y4, Y, _L);
    BBE_SVNX16_XP(Y5, Y, _L);
    BBE_DSELNX16I(Y7, Y6, XB, X3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y6, Y, _L);
    BBE_SVNX16_XP(Y7, Y, _L);
    BBE_DSELNX16I(Y9, Y8, XC, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y8, Y, back);
  }
} /* rbs3x3() */

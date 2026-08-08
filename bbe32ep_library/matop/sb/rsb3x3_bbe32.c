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

#define BBE_MOVVW(hi,lo,w) { hi=BBE_MOVVWH(w);lo=BBE_MOVVWL(w);}
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

/* M=3, N=3, Sx=9, Sy=16 */
void rsb3x3 ( int16_t * restrict y, const int16_t * restrict x, int L )
{
  const xb_vecNx16* restrict X = (const xb_vecNx16*)x;
  xb_vecNx16* restrict Y = (xb_vecNx16*)y;

  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
  xb_vecNx40 W0, W1;

  int k, _2L;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  _2L = (L << 1);
  if (L <= 0) return;
  for (k = 0; k<(L >> LOG2_BBE_SIMD_WIDTH); k++)
  {
    BBE_LVNX16_XP(Y0, X, _2L);
    BBE_LVNX16_XP(Y1, X, _2L);
    BBE_DSELNX16I(X8, X0, Y1, Y0, BBE_DSELI_INTERLEAVE_1);
    W0 = BBE_MOVWV(X8, X0);
    BBE_LVNX16_XP(Y2, X, _2L);
    BBE_LVNX16_XP(Y3, X, _2L);
    BBE_DSELNX16I(X9, X1, Y3, Y2, BBE_DSELI_INTERLEAVE_1);
    W1 = BBE_MOVWV(X9, X1);
    BBE_LVNX16_XP(Y4, X, _2L);
    BBE_LVNX16_XP(Y5, X, _2L);
    BBE_DSELNX16I(XA, X2, Y5, Y4, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(Y6, X, _2L);
    BBE_LVNX16_XP(Y7, X, _2L);
    BBE_DSELNX16I(XB, X3, Y7, Y6, BBE_DSELI_INTERLEAVE_1);
    Y8 = BBE_LVNX16_I(X, 0);
    BBE_DSELNX16I(XC, X4, Y9, Y8, BBE_DSELI_INTERLEAVE_1);
    X = X + (2 * BBE_SIMD_WIDTH - 8 * _2L) / sizeof(*X);

    BBE_MOVVW(X8, X0, W0);
    BBE_DSELNX16I(Y4, Y0, X4, X0, BBE_DSELI_INTERLEAVE_2);
    BBE_MOVVW(X9, X1, W1);
    BBE_DSELNX16I(Y5, Y1, X5, X1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(Y6, Y2, X6, X2, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(Y7, Y3, X7, X3, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(YC, Y8, XC, X8, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(YD, Y9, XD, X9, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(YE, YA, XE, XA, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(YF, YB, XF, XB, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(X2, X0, Y2, Y0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X3, X1, Y3, Y1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X6, X4, Y6, Y4, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(X7, X5, Y7, Y5, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(XA, X8, YA, Y8, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(XB, X9, YB, Y9, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(XE, XC, YE, YC, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(XF, XD, YF, YD, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(Y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(Y2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y3, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(Y5, Y4, X5, X4, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(Y4, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y5, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(Y7, Y6, X7, X6, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(Y6, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y7, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(Y9, Y8, X9, X8, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(Y8, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Y9, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(YB, YA, XB, XA, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(YA, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(YB, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(YD, YC, XD, XC, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(YC, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(YD, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(YF, YE, XF, XE, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(YE, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(YF, Y, 2 * BBE_SIMD_WIDTH);
  }
} /* rsb3x3() */

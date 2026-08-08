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

/* M=4, N=4, Sx=16, Sy=16 */
void rsb4x4 ( int16_t * restrict y, const int16_t * restrict x, int L )
{
  const xb_vecNx16* restrict X = (const xb_vecNx16*)x;
  xb_vecNx16* restrict Y = (xb_vecNx16*)y;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;
  xb_vecNx16 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, ya, yb, yc, yd, ye, yf;
  xb_vecNx40 w0, w1;
  int k, _2L, back;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  _2L = (L << 1);
  back = 2 * BBE_SIMD_WIDTH - 15 * _2L;
  if (L <= 0) return;
  for (k = 0; k<(L >> LOG2_BBE_SIMD_WIDTH); k++)
  {
    BBE_LVNX16_XP(y0, X, _2L);
    BBE_LVNX16_XP(y1, X, _2L);
    BBE_DSELNX16I(x8, x0, y1, y0, BBE_DSELI_INTERLEAVE_1);
    w0 = BBE_MOVWV(x8, x0);
    BBE_LVNX16_XP(y2, X, _2L);
    BBE_LVNX16_XP(y3, X, _2L);
    BBE_DSELNX16I(x9, x1, y3, y2, BBE_DSELI_INTERLEAVE_1);
    w1 = BBE_MOVWV(x9, x1);
    BBE_LVNX16_XP(y4, X, _2L);
    BBE_LVNX16_XP(y5, X, _2L);
    BBE_DSELNX16I(xa, x2, y5, y4, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(y6, X, _2L);
    BBE_LVNX16_XP(y7, X, _2L);
    BBE_DSELNX16I(xb, x3, y7, y6, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(y8, X, _2L);
    BBE_LVNX16_XP(y9, X, _2L);
    BBE_DSELNX16I(xc, x4, y9, y8, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(ya, X, _2L);
    BBE_LVNX16_XP(yb, X, _2L);
    BBE_DSELNX16I(xd, x5, yb, ya, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(yc, X, _2L);
    BBE_LVNX16_XP(yd, X, _2L);
    BBE_DSELNX16I(xe, x6, yd, yc, BBE_DSELI_INTERLEAVE_1);
    BBE_LVNX16_XP(ye, X, _2L);
    BBE_LVNX16_XP(yf, X, back);
    BBE_DSELNX16I(xf, x7, yf, ye, BBE_DSELI_INTERLEAVE_1);

    BBE_MOVVW(x8, x0, w0);
    BBE_DSELNX16I(y4, y0, x4, x0, BBE_DSELI_INTERLEAVE_2);
    BBE_MOVVW(x9, x1, w1);
    BBE_DSELNX16I(y5, y1, x5, x1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y6, y2, x6, x2, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(y7, y3, x7, x3, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(yc, y8, xc, x8, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(yd, y9, xd, x9, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(ye, ya, xe, xa, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(yf, yb, xf, xb, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(x2, x0, y2, y0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x3, x1, y3, y1, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x6, x4, y6, y4, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(x7, x5, y7, y5, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(xa, x8, ya, y8, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(xb, x9, yb, y9, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(xe, xc, ye, yc, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELNX16I(xf, xd, yf, yd, BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(y1, y0, x1, x0, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y1, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(y3, y2, x3, x2, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(y2, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y3, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(y5, y4, x5, x4, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(y4, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y5, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(y7, y6, x7, x6, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(y6, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y7, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(y9, y8, x9, x8, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(y8, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(y9, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(yb, ya, xb, xa, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(ya, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(yb, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(yd, yc, xd, xc, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(yc, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(yd, Y, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(yf, ye, xf, xe, BBE_DSELI_INTERLEAVE_2);
    BBE_SVNX16_IP(ye, Y, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(yf, Y, 2 * BBE_SIMD_WIDTH);
  }
} /* rsb4x4() */

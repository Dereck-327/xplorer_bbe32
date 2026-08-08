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

/* M=8, N=1, Sx=8, Sy=8 */
void rbs8x1 ( int16_t * restrict y, const int16_t * restrict x, int L )
{
  const xb_vecNx16* restrict w = (const xb_vecNx16*)x;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  int l;
  xb_vecNx16* restrict z1 = z + (L >> 4);
  xb_vecNx16* restrict z2 = z + (L >> 3);
  xb_vecNx16* restrict z3 = z1 + (L >> 3);
  xb_vecNx16* restrict z4 = z2 + (L >> 3);
  xb_vecNx16* restrict z5 = z3 + (L >> 3);
  xb_vecNx16* restrict z6 = z4 + (L >> 3);
  xb_vecNx16* restrict z7 = z5 + (L >> 3);
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  if (L <= 0) return;
  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    BBE_LVNX16_IP(Y0, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y2, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y3, w, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(Y4, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y5, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y6, w, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y7, w, 2 * BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(X1, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X4, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_SVNX16_IP(X0, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1, z1, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2, z2, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X3, z3, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X4, z4, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X5, z5, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X6, z6, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X7, z7, 2 * BBE_SIMD_WIDTH);
  }
} /* rbs8x1() */

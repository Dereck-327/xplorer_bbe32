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

#define sz_cf32 (int)sizeof(complex_float)

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
void cbs8x1f ( complex_float * restrict y, const complex_float * restrict x, int L )
{
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
  int l;
  xb_vecNx16 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  int ix0, ix1, ix2;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);

  px0 = (const xb_vecNx16 *)(x);
  px1 = (const xb_vecNx16 *)(x+4);
  py0 = (      xb_vecNx16 *)(y+3*L);
  py1 = (      xb_vecNx16 *)(y+7*L);
  ix0 = -3*L*sz_cf32;
  ix1 = -2*L*sz_cf32;
  ix2 = -1*L*sz_cf32;

  /* Convert by 8 values */
  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
  {

    BBE_LVNX16_IP(X0, px0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, px0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, px0, 2*2*BBE_SIMD_WIDTH);

    BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

    BBE_SVNX16_X (Y0, py0, ix0);
    BBE_SVNX16_X (Y1, py0, ix1);
    BBE_SVNX16_X (Y2, py0, ix2);
    BBE_SVNX16_IP(Y3, py0, 2*BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(X0, px1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, px1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, px1, 2*2*BBE_SIMD_WIDTH);

    BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

    BBE_SVNX16_X (Y0, py1, ix0);
    BBE_SVNX16_X (Y1, py1, ix1);
    BBE_SVNX16_X (Y2, py1, ix2);
    BBE_SVNX16_IP(Y3, py1, 2*BBE_SIMD_WIDTH);
  }
} /* cbs8x1f() */

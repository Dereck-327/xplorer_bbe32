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

/* M=4, N=1, Sx=4, Sy=4 */
void csb4x1f ( complex_float * restrict y, const complex_float * restrict x, int L )
{
#if 1
    csb2x2f(y, x, L);
#else
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict px2;
  const xb_vecN_4xcf32 * restrict px3;
        xb_vecN_4xcf32 * restrict py;
  xb_vecN_4xcf32 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  int l;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)((complex_float *)px0+L);
  px2 = (const xb_vecN_4xcf32 *)((complex_float *)px1+L);
  px3 = (const xb_vecN_4xcf32 *)((complex_float *)px2+L);
  py  = (      xb_vecN_4xcf32 *)(y);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
  {
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px3, 2*BBE_SIMD_WIDTH);

    BBE_DSELN_4XCF32I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y2, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y3, py, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* csb4x1f() */

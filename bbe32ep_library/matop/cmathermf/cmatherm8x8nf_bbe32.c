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
DISCARD_FUN(void, cmatherm8x8nf,( complex_float * restrict y, 
                            const complex_float * restrict x, 
                            int L ))
#else
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

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     None
*/
void cmatherm8x8nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  const complex_float  * restrict px_row;
  const xb_vecN_4xcf32 * restrict px_col0;
  const xb_vecN_4xcf32 * restrict px_col1;
        xb_vecN_4xcf32 * restrict py;
  int l;
  int i, rowstride, colstride;

  xb_vecN_4xcf32 X00, X01, X02, X03, X04, X05, X06, X07,
                 X10, X11, X12, X13, X14, X15, X16, X17;
  xb_vecN_4xcf32 Xc0, Xc1, Xc2, Xc3, Xc4, Xc5, Xc6, Xc7;
  xb_vecN_4xcf32 Y0, Y1;
  xb_vecN_4xcf32 y00, y01, y10, y11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  px_row  = (const complex_float  *)(x);
  px_col0 = (const xb_vecN_4xcf32 *)(x);
  px_col1 = (const xb_vecN_4xcf32 *)(x+4);
  py      = (      xb_vecN_4xcf32 *)(y);
  i = 0;

  /* Compute matrices by 1 row per iteration */
  __Pragma("loop_count factor=4");
  for (l=0; l<L*8; l++)
  {
    /* i=(i+1)%8; */
    i = BBE_ADDMOD16U(i, 0x080001);
    /* rowstride= (i==0)? 0 : -14*2*BBE_SIMD_WIDTH; */
    rowstride = -14*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(rowstride, i, i);
    /* colstride= (i==0)? 0 : -16*2*BBE_SIMD_WIDTH; */
    colstride = -16*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(colstride, i, i);

    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X00, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X01, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X02, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X03, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X04, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X05, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X06, px_col0, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X07, px_col0, 2*2*BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(X10, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X11, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X12, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X13, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X14, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X15, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X16, px_col1, 2*2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X17, px_col1, 2*2*BBE_SIMD_WIDTH);

    px_col0 = (const xb_vecN_4xcf32 *)((intptr_t)px_col0 + colstride);
    px_col1 = (const xb_vecN_4xcf32 *)((intptr_t)px_col1 + colstride);

    /* Load 1 column of matrix X */
    BBE_LSN_4XCF32_IP(Xc0, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc1, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc2, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc3, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc4, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc5, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc6, px_row, 2*2*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_IP(Xc7, px_row, sizeof(complex_float));
    px_row = (const complex_float *)((intptr_t)px_row + rowstride);
    Xc0 = BBE_REPN_4XCF32(Xc0, 0);
    Xc1 = BBE_REPN_4XCF32(Xc1, 0);
    Xc2 = BBE_REPN_4XCF32(Xc2, 0);
    Xc3 = BBE_REPN_4XCF32(Xc3, 0);
    Xc4 = BBE_REPN_4XCF32(Xc4, 0);
    Xc5 = BBE_REPN_4XCF32(Xc5, 0);
    Xc6 = BBE_REPN_4XCF32(Xc6, 0);
    Xc7 = BBE_REPN_4XCF32(Xc7, 0);

    /* Compute 1 row of Hermitian matrix Y=X'*X */
    y00 = BBE_MULMN_4XCF32( Xc0, X00, 0, 0x4);
    y01 = BBE_MULMN_4XCF32( Xc0, X00, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc1, X01, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc1, X01, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc2, X02, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc2, X02, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc3, X03, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc3, X03, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc4, X04, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc4, X04, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc5, X05, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc5, X05, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc6, X06, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc6, X06, 2, 0xB);
    BBE_MULMASN_4XCF32(y00, Xc7, X07, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, Xc7, X07, 2, 0xB);

    y10 = BBE_MULMN_4XCF32( Xc0, X10, 0, 0x4);
    y11 = BBE_MULMN_4XCF32( Xc0, X10, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc1, X11, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc1, X11, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc2, X12, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc2, X12, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc3, X13, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc3, X13, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc4, X14, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc4, X14, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc5, X15, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc5, X15, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc6, X16, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc6, X16, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, Xc7, X17, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, Xc7, X17, 2, 0xB);

    Y0 = BBE_ADDN_4XCF32(y00, y01);
    Y1 = BBE_ADDN_4XCF32(y10, y11);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
  }
} /* cmatherm8x8nf() */
#endif

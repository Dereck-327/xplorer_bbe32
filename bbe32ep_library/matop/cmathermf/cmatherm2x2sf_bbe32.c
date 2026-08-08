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
DISCARD_FUN(void, cmatherm2x2sf,( complex_float * restrict y, 
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

/* Streaming Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm2x2sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
        xb_vecN_4xcf32 * restrict py0;
        xb_vecN_4xcf32 * restrict py1;
  int l;

  xb_vecN_4xcf32 x00, x01, x10, x11;
  xb_vecN_4xcf32 y00, y01, y10, y11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);

  px0 = (const xb_vecN_4xcf32 *)(x+L*2);
  px1 = (const xb_vecN_4xcf32 *)(x+L*3);
  py0 = (      xb_vecN_4xcf32 *)(y+L*2);
  py1 = (      xb_vecN_4xcf32 *)(y+L*3);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
  {
    /* Load input matrices X and Y */
    x00 = BBE_LVN_4XCF32_X(px0, -2*L*(int)sizeof(complex_float));
    x01 = BBE_LVN_4XCF32_X(px1, -2*L*(int)sizeof(complex_float));
    BBE_LVN_4XCF32_IP(x10, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(x11, px1, 2*BBE_SIMD_WIDTH);

    /* Multiply matrices X' and X */
    y00 = BBE_MULJN_4XCF32(x00, x00);
    BBE_MULJAN_4XCF32(y00, x10, x10);
    y01 = BBE_MULJN_4XCF32(x01, x00);
    BBE_MULJAN_4XCF32(y01, x11, x10);
    y10 = BBE_CONJN_4XCF32(y01);
    y11 = BBE_MULJN_4XCF32(x01, x01);
    BBE_MULJAN_4XCF32(y11, x11, x11);

    /* Save results */
    BBE_SVN_4XCF32_X (y00, py0, -2*L*(int)sizeof(complex_float));
    BBE_SVN_4XCF32_X (y01, py1, -2*L*(int)sizeof(complex_float));
    BBE_SVN_4XCF32_IP(y10, py0, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(y11, py1, 2*BBE_SIMD_WIDTH);
  }
} /* cmatherm2x2sf() */
#endif

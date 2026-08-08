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
DISCARD_FUN(void, cmatherm2x2nf,( complex_float * restrict y, 
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

/* Block Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     None
*/
void cmatherm2x2nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 X, Y;
  xb_vecN_4xcf32 x00, x01, x10, x11, y0, y1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_4xcf32 *)x;
  py = (      xb_vecN_4xcf32 *)y;

  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X, px, 2*BBE_SIMD_WIDTH);

    BBE_DSELN_4XCF32I(x01, x00, X, X, BBE_DSELI_INTERLEAVE_4);
    x10 = BBE_SHFLN_4XCF32I(X, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    x11 = BBE_SHFLN_4XCF32I(X, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

    /* Multiply matrices X' and X */
    y0 = BBE_MULMN_4XCF32( x00, x10, 0, 0x4);
    y1 = BBE_MULMN_4XCF32( x00, x10, 2, 0xB);
    BBE_MULMASN_4XCF32(y0, x01, x11, 0, 0x4);
    BBE_MULMASN_4XCF32(y1, x01, x11, 2, 0xB);
    Y = BBE_ADDN_4XCF32(y0, y1);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
  }
} /* cmatherm2x2nf() */
#endif

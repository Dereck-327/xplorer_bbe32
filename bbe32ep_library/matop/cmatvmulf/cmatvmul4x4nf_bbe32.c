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
    Complex Matrix-Matrix/Matrix-Vector Multiply
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
DISCARD_FUN(void, cmatvmul4x4nf,( complex_float * restrict z, 
                            const complex_float * restrict x, 
                            const complex_float * restrict y, 
                            int L ))
#else
/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     None
*/
void cmatvmul4x4nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3;
  xb_vecN_4xcf32 Y;
  xb_vecN_4xcf32 Z, Z0, Z1, Z2, Z3;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x+8);
  py  = (const xb_vecN_4xcf32 *)(y);
  pz  = (      xb_vecN_4xcf32 *)(z);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_XP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(X1, px0, 6*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(X2, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(X3, px1, 6*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(Y , py, 2*BBE_SIMD_WIDTH);

    /* Multiply input matrices X and Y */
    Z0 = BBE_MULN_4XCF32(X0, Y);
    Z1 = BBE_MULN_4XCF32(X1, Y);
    Z2 = BBE_MULN_4XCF32(X2, Y);
    Z3 = BBE_MULN_4XCF32(X3, Y);

    BBE_DSELN_4XCF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Z3, Z1, Z3, Z1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Z1, Z0, Z1, Z0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Z3, Z2, Z3, Z2, BBE_DSELI_INTERLEAVE_4);

    Z0 = BBE_ADDN_4XCF32(Z0, Z1);
    Z2 = BBE_ADDN_4XCF32(Z2, Z3);
    Z  = BBE_ADDN_4XCF32(Z0, Z2);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
} /* cmatvmul4x4nf() */
#endif

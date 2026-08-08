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
DISCARD_FUN(void, cmatvmul2x2nf,( complex_float * restrict z, 
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

/* Block Order, Floating-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be even
*/
void cmatvmul2x2nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L )
{
#if 1
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, Y, Z;
  xb_vecN_4xcf32 x0, x1, y0, y1, z0, z1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_4xcf32 *)x;
  py = (const xb_vecN_4xcf32 *)y;
  pz = (      xb_vecN_4xcf32 *)z;

  for (l=0; l<(L>>1); l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y, py, 2*BBE_SIMD_WIDTH);

    x0 = BBE_SELN_4XCF32I(X1, X0, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    x1 = BBE_SELN_4XCF32I(X1, X0, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    y0 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_DUPLICATE_4_EVEN);
    y1 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_DUPLICATE_4_ODD);

    /* Multiply input matrices X and Y */
    z0 = BBE_MULMN_4XCF32( x0, y0, 0, 0x4);
    z1 = BBE_MULMN_4XCF32( x0, y0, 1, 0xB);
    BBE_MULMASN_4XCF32(z0, x1, y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z1, x1, y1, 1, 0xB);
    Z = BBE_ADDN_4XCF32(z0, z1);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, Y, Z;
  xb_vecN_4xcf32 x0, x1, y0, y1, z0, z1, t0, t1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_4xcf32 *)x;
  py = (const xb_vecN_4xcf32 *)y;
  pz = (      xb_vecN_4xcf32 *)z;

  for (l=0; l<(L>>1); l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y, py, 2*BBE_SIMD_WIDTH);

    x0 = X0;
    x1 = X1;
    y0 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    y1 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

    /* Multiply input matrices X and Y */
    t0 = BBE_MULN_4XCF32(x0, y0);
    t1 = BBE_MULN_4XCF32(x1, y1);
    z0 = BBE_SELN_4XCF32I(t1, t0, BBE_SELI_EXTRACT_4_OF_8_OFF_0);
    z1 = BBE_SELN_4XCF32I(t1, t0, BBE_SELI_EXTRACT_4_OF_8_OFF_4);
    Z = BBE_ADDN_4XCF32(z0, z1);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatvmul2x2nf() */
#endif

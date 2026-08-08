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
DISCARD_FUN(void, cmatmul4x4nf,( complex_float * restrict z, 
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

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatmul4x4nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L )
{
#if 0
  const complex_float  * restrict px0;
  const complex_float  * restrict px1;
  const complex_float  * restrict px2;
  const complex_float  * restrict px3;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  xb_vecN_4xcf32 Z0, Z1, Z2, Z3;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px0 = (const complex_float  *)(x+0);
  px1 = (const complex_float  *)(x+4);
  px2 = (const complex_float  *)(x+8);
  px3 = (const complex_float  *)(x+12);
  py  = (const xb_vecN_4xcf32 *)(y);
  pz  = (      xb_vecN_4xcf32 *)(z);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LSN_4XCF32_IP(X0, px0, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X1, px1, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X2, px2, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X3, px3, sizeof(complex_float));
    BBE_LVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    X0 = BBE_REPN_4XCF32(X0, 0);
    X1 = BBE_REPN_4XCF32(X1, 0);
    X2 = BBE_REPN_4XCF32(X2, 0);
    X3 = BBE_REPN_4XCF32(X3, 0);
    Z0 = BBE_MULN_4XCF32(X0, Y0);
    Z1 = BBE_MULN_4XCF32(X1, Y0);
    Z2 = BBE_MULN_4XCF32(X2, Y0);
    Z3 = BBE_MULN_4XCF32(X3, Y0);

    BBE_LSN_4XCF32_IP(X0, px0, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X1, px1, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X2, px2, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X3, px3, sizeof(complex_float));
    BBE_LVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
    X0 = BBE_REPN_4XCF32(X0, 0);
    X1 = BBE_REPN_4XCF32(X1, 0);
    X2 = BBE_REPN_4XCF32(X2, 0);
    X3 = BBE_REPN_4XCF32(X3, 0);
    BBE_MULAN_4XCF32(Z0, X0, Y1);
    BBE_MULAN_4XCF32(Z1, X1, Y1);
    BBE_MULAN_4XCF32(Z2, X2, Y1);
    BBE_MULAN_4XCF32(Z3, X3, Y1);

    BBE_LSN_4XCF32_IP(X0, px0, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X1, px1, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X2, px2, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X3, px3, sizeof(complex_float));
    BBE_LVN_4XCF32_IP(Y2, py, 2*BBE_SIMD_WIDTH);
    X0 = BBE_REPN_4XCF32(X0, 0);
    X1 = BBE_REPN_4XCF32(X1, 0);
    X2 = BBE_REPN_4XCF32(X2, 0);
    X3 = BBE_REPN_4XCF32(X3, 0);
    BBE_MULAN_4XCF32(Z0, X0, Y2);
    BBE_MULAN_4XCF32(Z1, X1, Y2);
    BBE_MULAN_4XCF32(Z2, X2, Y2);
    BBE_MULAN_4XCF32(Z3, X3, Y2);

    BBE_LSN_4XCF32_IP(X0, px0, 13*sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X1, px1, 13*sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X2, px2, 13*sizeof(complex_float));
    BBE_LSN_4XCF32_IP(X3, px3, 13*sizeof(complex_float));
    BBE_LVN_4XCF32_IP(Y3, py, 2*BBE_SIMD_WIDTH);
    X0 = BBE_REPN_4XCF32(X0, 0);
    X1 = BBE_REPN_4XCF32(X1, 0);
    X2 = BBE_REPN_4XCF32(X2, 0);
    X3 = BBE_REPN_4XCF32(X3, 0);
    BBE_MULAN_4XCF32(Z0, X0, Y3);
    BBE_MULAN_4XCF32(Z1, X1, Y3);
    BBE_MULAN_4XCF32(Z2, X2, Y3);
    BBE_MULAN_4XCF32(Z3, X3, Y3);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z1, pz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z2, pz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z3, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;
  int i, ystride;

  xb_vecN_4xcf32 X0, X1, Y0, Y1, Y2, Y3;
  xb_vecN_4xcf32 x00, x01, x02, x03, x10, x11, x12, x13;
  xb_vecN_4xcf32 Z0, Z1;
  xb_vecN_4xcf32 z00, z01, z10, z11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x);
  py  = (const xb_vecN_4xcf32 *)(y);
  pz  = (      xb_vecN_4xcf32 *)(z);
  i = 0;

  __Pragma("loop_count factor=2");
  for (l=0; l<L*2; l++)
  {
    i = XT_XOR(i, 1);
    /* ystride= (i==0)? 0 : -4*2*BBE_SIMD_WIDTH; */
    ystride = -4 * 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(ystride, i, i);

    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y2, py, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y3, py, 2*BBE_SIMD_WIDTH);
    py = (const xb_vecN_4xcf32 *)((intptr_t)py + ystride);

    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 2*BBE_SIMD_WIDTH);
    x00 = BBE_REPN_4XCF32(X0, 0);
    x10 = BBE_REPN_4XCF32(X1, 0);
    x01 = BBE_REPN_4XCF32(X0, 1);
    x11 = BBE_REPN_4XCF32(X1, 1);
    BBE_LVN_4XCF32_IP(X0, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px1, 2*BBE_SIMD_WIDTH);
    x02 = BBE_REPN_4XCF32(X0, 2);
    x12 = BBE_REPN_4XCF32(X1, 2);
    x03 = BBE_REPN_4XCF32(X0, 3);
    x13 = BBE_REPN_4XCF32(X1, 3);

    z00 = BBE_MULMN_4XCF32( x00, Y0, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x00, Y0, 1, 0xB);
    z10 = BBE_MULMN_4XCF32( x10, Y0, 0, 0x4);
    z11 = BBE_MULMN_4XCF32( x10, Y0, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, Y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, Y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x11, Y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x11, Y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, Y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, Y2, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x12, Y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x12, Y2, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x03, Y3, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x03, Y3, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x13, Y3, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x13, Y3, 1, 0xB);
    Z0 = BBE_ADDN_4XCF32(z00, z01);
    Z1 = BBE_ADDN_4XCF32(z10, z11);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z1, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatmul4x4nf() */
#endif

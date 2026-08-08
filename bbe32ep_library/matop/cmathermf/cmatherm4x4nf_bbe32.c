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
DISCARD_FUN(void, cmatherm4x4nf,( complex_float * restrict y, 
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

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     None
*/
void cmatherm4x4nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
#if 0
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3, Y0, Y1;
  xb_vecN_4xcf32 t0, t1;
  xb_vecN_4xcf32 x00, x01, x02, x03;
  xb_vecN_4xcf32 x10, x11, x12, x13;
  xb_vecN_4xcf32 y00, y01, y10, y11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  if (L <= 0) return;

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x);
  py  = (      xb_vecN_4xcf32 *)(y);

  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px0, 2*BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(t0, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(t1, px1, 6*BBE_SIMD_WIDTH);

    x00 = BBE_REPN_4XCF32(t0, 0);
    x10 = BBE_REPN_4XCF32(t0, 1);
    x01 = BBE_REPN_4XCF32(t1, 0);
    x11 = BBE_REPN_4XCF32(t1, 1);
    x02 = BBE_REPN_4XCF32(X2, 0);
    x12 = BBE_REPN_4XCF32(X2, 1);
    x03 = BBE_REPN_4XCF32(X3, 0);
    x13 = BBE_REPN_4XCF32(X3, 1);

    /* Multiply matrices X' and X */
    y00 = BBE_MULMN_4XCF32( x00, X0, 0, 0x4);
    y01 = BBE_MULMN_4XCF32( x00, X0, 2, 0xB);
    y10 = BBE_MULMN_4XCF32( x10, X0, 0, 0x4);
    y11 = BBE_MULMN_4XCF32( x10, X0, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x01, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x01, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x11, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x11, X1, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x02, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x02, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x12, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x12, X2, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x03, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x03, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x13, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x13, X3, 2, 0xB);

    /* Save results */
    Y0 = BBE_ADDN_4XCF32(y00, y01);
    Y1 = BBE_ADDN_4XCF32(y10, y11);
    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 6*BBE_SIMD_WIDTH);
  }

  /*
   * Compute 2-nd part of matrices
   */

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x);
  py  = (      xb_vecN_4xcf32 *)(y+8);

  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px0, 2*BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(t0, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(t1, px1, 6*BBE_SIMD_WIDTH);

    x00 = BBE_REPN_4XCF32(t0, 2);
    x10 = BBE_REPN_4XCF32(t0, 3);
    x01 = BBE_REPN_4XCF32(t1, 2);
    x11 = BBE_REPN_4XCF32(t1, 3);
    x02 = BBE_REPN_4XCF32(X2, 2);
    x12 = BBE_REPN_4XCF32(X2, 3);
    x03 = BBE_REPN_4XCF32(X3, 2);
    x13 = BBE_REPN_4XCF32(X3, 3);

    /* Multiply matrices X' and X */
    y00 = BBE_MULMN_4XCF32( x00, X0, 0, 0x4);
    y01 = BBE_MULMN_4XCF32( x00, X0, 2, 0xB);
    y10 = BBE_MULMN_4XCF32( x10, X0, 0, 0x4);
    y11 = BBE_MULMN_4XCF32( x10, X0, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x01, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x01, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x11, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x11, X1, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x02, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x02, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x12, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x12, X2, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x03, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x03, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x13, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x13, X3, 2, 0xB);

    /* Save results */
    Y0 = BBE_ADDN_4XCF32(y00, y01);
    Y1 = BBE_ADDN_4XCF32(y10, y11);
    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 6*BBE_SIMD_WIDTH);
  }
#elif 0
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3,
                 Y0, Y1, Y2, Y3;
  xb_vecN_4xcf32 t0, t1, t2, t3;
  xb_vecN_4xcf32 x00, x01, x02, x03;
  xb_vecN_4xcf32 x10, x11, x12, x13;
  xb_vecN_4xcf32 x20, x21, x22, x23;
  xb_vecN_4xcf32 x30, x31, x32, x33;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)(x);
  px1 = (const xb_vecN_4xcf32 *)(x);
  py  = (      xb_vecN_4xcf32 *)(y);

  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px0, 2*BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(t0, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(t1, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(t2, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(t3, px1, 2*BBE_SIMD_WIDTH);

    x00 = BBE_REPN_4XCF32(t0, 0);
    x10 = BBE_REPN_4XCF32(X0, 1);
    x20 = BBE_REPN_4XCF32(t0, 2);
    x30 = BBE_REPN_4XCF32(X0, 3);
    x01 = BBE_REPN_4XCF32(t1, 0);
    x11 = BBE_REPN_4XCF32(X1, 1);
    x21 = BBE_REPN_4XCF32(t1, 2);
    x31 = BBE_REPN_4XCF32(X1, 3);
    x02 = BBE_REPN_4XCF32(t2, 0);
    x12 = BBE_REPN_4XCF32(X2, 1);
    x22 = BBE_REPN_4XCF32(t2, 2);
    x32 = BBE_REPN_4XCF32(X2, 3);
    x03 = BBE_REPN_4XCF32(t3, 0);
    x13 = BBE_REPN_4XCF32(X3, 1);
    x23 = BBE_REPN_4XCF32(t3, 2);
    x33 = BBE_REPN_4XCF32(X3, 3);

    /* Multiply matrices X' and X */
    Y0 = BBE_MULMN_4XCF32( x00, X0, 0, 0x4);
    Y1 = BBE_MULMN_4XCF32( x10, X0, 0, 0x4);
    Y2 = BBE_MULMN_4XCF32( x20, X0, 0, 0x4);
    Y3 = BBE_MULMN_4XCF32( x30, X0, 0, 0x4);
    BBE_MULMASN_4XCF32(Y0, x00, X0, 2, 0xB);
    BBE_MULMASN_4XCF32(Y1, x10, X0, 2, 0xB);
    BBE_MULMASN_4XCF32(Y2, x20, X0, 2, 0xB);
    BBE_MULMASN_4XCF32(Y3, x30, X0, 2, 0xB);

    BBE_MULMASN_4XCF32(Y0, x01, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(Y1, x11, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(Y2, x21, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(Y3, x31, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(Y0, x01, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(Y1, x11, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(Y2, x21, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(Y3, x31, X1, 2, 0xB);

    BBE_MULMASN_4XCF32(Y0, x02, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(Y1, x12, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(Y2, x22, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(Y3, x32, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(Y0, x02, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(Y1, x12, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(Y2, x22, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(Y3, x32, X2, 2, 0xB);

    BBE_MULMASN_4XCF32(Y0, x03, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(Y1, x13, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(Y2, x23, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(Y3, x33, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(Y0, x03, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(Y1, x13, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(Y2, x23, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(Y3, x33, X3, 2, 0xB);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y2, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y3, py, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict ptmp;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, X3, Y0, Y1;
  xb_vecNx16     X0_,X1_,X2_,X3_;
  xb_vecN_4xcf32 x00, x01, x02, x03;
  xb_vecN_4xcf32 x10, x11, x12, x13;
  xb_vecN_4xcf32 y00, y01, y10, y11;
  xb_vecNx16 vTmp0, vTmp1, vTmp2, vTmp3;
  vselN sel0, sel1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_4xcf32 *)x;
  px1 = (const xb_vecN_4xcf32 *)x;
  py  = (      xb_vecN_4xcf32 *)y;

  /* Initialize select registers */
  vTmp1 = BBE_SEQNX16();
  vTmp0 = BBE_SHFLNX16I(vTmp1, BBE_SHFLI_REP_0X4);
  vTmp1 = BBE_SHFLNX16I(vTmp1, BBE_SHFLI_REP_1X4);
  sel0 = BBE_MOVVSELNX16(vTmp0, 0);
  sel1 = BBE_MOVVSELNX16(vTmp1, 0);

  /* Compute matrices by 2 rows per iteration */
  __Pragma("loop_count factor=2");
  for (l=0; l<L*2; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X2, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X3, px0, 2*BBE_SIMD_WIDTH);
    ptmp = px0;
    px0 = px1;
    px1 = ptmp;

    /* Replicate 2 columns of matrix X */
    X0_ = BBE_MOVNX16_FROMN_4XCF32(X0);
    X1_ = BBE_MOVNX16_FROMN_4XCF32(X1);
    X2_ = BBE_MOVNX16_FROMN_4XCF32(X2);
    X3_ = BBE_MOVNX16_FROMN_4XCF32(X3);
    vTmp0 = BBE_SHFLNX16(X0_, sel0);
    vTmp1 = BBE_SHFLNX16(X1_, sel0);
    vTmp2 = BBE_SHFLNX16(X2_, sel0);
    BBE_SHFLUNX16(vTmp3, X3_, sel0, 8);
    x00 = BBE_MOVN_4XCF32_FROMNX16(vTmp0);
    x01 = BBE_MOVN_4XCF32_FROMNX16(vTmp1);
    x02 = BBE_MOVN_4XCF32_FROMNX16(vTmp2);
    x03 = BBE_MOVN_4XCF32_FROMNX16(vTmp3);
    vTmp0 = BBE_SHFLNX16(X0_, sel1);
    vTmp1 = BBE_SHFLNX16(X1_, sel1);
    vTmp2 = BBE_SHFLNX16(X2_, sel1);
    BBE_SHFLUNX16(vTmp3, X3_, sel1, 8);
    x10 = BBE_MOVN_4XCF32_FROMNX16(vTmp0);
    x11 = BBE_MOVN_4XCF32_FROMNX16(vTmp1);
    x12 = BBE_MOVN_4XCF32_FROMNX16(vTmp2);
    x13 = BBE_MOVN_4XCF32_FROMNX16(vTmp3);

    /* Compute 2 rows of Hermitian matrix Y=X'*X */
    y00 = BBE_MULMN_4XCF32( x00, X0, 0, 0x4);
    y01 = BBE_MULMN_4XCF32( x00, X0, 2, 0xB);
    y10 = BBE_MULMN_4XCF32( x10, X0, 0, 0x4);
    y11 = BBE_MULMN_4XCF32( x10, X0, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x01, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x01, X1, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x11, X1, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x11, X1, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x02, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x02, X2, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x12, X2, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x12, X2, 2, 0xB);

    BBE_MULMASN_4XCF32(y00, x03, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y01, x03, X3, 2, 0xB);
    BBE_MULMASN_4XCF32(y10, x13, X3, 0, 0x4);
    BBE_MULMASN_4XCF32(y11, x13, X3, 2, 0xB);

    Y0 = BBE_ADDN_4XCF32(y00, y01);
    Y1 = BBE_ADDN_4XCF32(y10, y11);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y0, py, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Y1, py, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatherm4x4nf() */
#endif

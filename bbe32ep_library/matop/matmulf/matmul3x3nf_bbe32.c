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
    Real Matrix-Matrix/Matrix-Vector Multiply
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
DISCARD_FUN(void, matmul3x3nf,( float32_t * restrict z, 
                          const float32_t * restrict x, 
                          const float32_t * restrict y, 
                          int L ))
#else

#ifndef BBE_SELN_2XF32
#define BBE_SELN_2XF32(a, b, c) (BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(a), BBE_MOVNX16_FROMN_2XF32(b), (c) )))
#endif
/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/                                                                                                                                                      
void matmul3x3nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
#if 0
  static const uint16_t ALIGN(32) aSel[BBE_SIMD_WIDTH]=
  {	0x4,  0x5,  0x4,  0x5,  0x4,  0x5,
    0xa,  0xb,  0xa,  0xb,  0xa,  0xb,
    0x10, 0x11, 0x10, 0x11 };

  const xb_vecN_2xf32 * restrict vpx;
  const xb_vecN_2xf32 * restrict vpy;
        xb_vecN_2xf32 * restrict vpz;
  const float32_t * restrict spx;
  const float32_t * restrict spy;
        float32_t * restrict spz;
  int l;

  xb_vecN_2xf32 X0, X1, Y0, Y1, Z0, Z1;
  xb_vecN_2xf32 x00, x01, x02, y00, y01, y02;
  xb_vecN_2xf32 x10, x11, x12, y10, y11, y12;
  xb_vecNx16 vTmp;
  vselN x02_sel;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  vpx = (const xb_vecN_2xf32 *)x;
  vpy = (const xb_vecN_2xf32 *)y;
  vpz = (      xb_vecN_2xf32 *)z;
  spx = (const float32_t *)(vpx+1);
  spy = (const float32_t *)(vpy+1);
  spz = (      float32_t *)(vpz+1);
  
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,0);
  x02_sel = BBE_MOVVSELNX16(vTmp,0);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, vpx, 4*BBE_SIMD_WIDTH);
    BBE_LSN_2XF32_XP(X1, spx, 4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y0, vpy, 4*BBE_SIMD_WIDTH);
    BBE_LSN_2XF32_XP(Y1, spy, 4*BBE_SIMD_WIDTH);

    x00 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_1);
    x01 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_2);
    x02 = BBE_SELN_2XF32(X1, X0, x02_sel);
    x10 = BBE_MOVN_2XF32_FROMF32(BBE_SELSN_2XF32(X0, 6));
    x11 = BBE_MOVN_2XF32_FROMF32(BBE_SELSN_2XF32(X0, 7));
    x12 = X1;

    y00 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_1);
    y01 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_2);
    y02 = BBE_SELN_2XF32I(Y1, Y0, BBE_SELI_MMC3X3X3X3_OFFSET_M2_0_STEP_3);
    y10 = BBE_MOVN_2XF32_FROMF32(BBE_SELSN_2XF32(Y0, 2));
    y11 = BBE_MOVN_2XF32_FROMF32(BBE_SELSN_2XF32(Y0, 5));
    y12 = Y1;

    /* Multiply input matrices X and Y */
    Z0 = BBE_MULN_2XF32(x00, y00);
    BBE_MULAN_2XF32(Z0, x01, y01);
    BBE_MULAN_2XF32(Z0, x02, y02);
    
    Z1 = BBE_MULN_2XF32(x10, y10);
    BBE_MULAN_2XF32(Z1, x11, y11);
    BBE_MULAN_2XF32(Z1, x12, y12);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, vpz, 4*BBE_SIMD_WIDTH);
    BBE_SSN_2XF32_XP(Z1, spz, 4*BBE_SIMD_WIDTH);
  }
#else
  static const uint16_t ALIGN(32) aSel[BBE_SIMD_WIDTH]=
  {	0x4,  0x5,  0x4,  0x5,  0x4,  0x5,
    0xa,  0xb,  0xa,  0xb,  0xa,  0xb,
    0x10, 0x11, 0x10, 0x11 };

  const xb_vecN_2xf32 * restrict vpx;
  const xb_vecN_2xf32 * restrict vpy;
        xb_vecN_2xf32 * restrict vpz;
        xb_vecN_2xf32 * restrict spz;
  const float32_t * restrict spx;
  const float32_t * restrict spy;
  int l;

  xb_vecN_2xf32 X0, X1, Y0, Y1, Z0, Z1;
  xb_vecN_2xf32 x00, x01, x02, y00, y01, y02;
  xb_vecN_2xf32 x10, x11, x12, y10, y11, y12;
  xb_vecNx16 vTmp;
  vselN_2 x02_sel;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  vpx = (const xb_vecN_2xf32 *)(x);
  vpy = (const xb_vecN_2xf32 *)(y);
  vpz = (      xb_vecN_2xf32 *)(z);
  spz = (      xb_vecN_2xf32 *)(z+8);
  spx = x+6;
  spy = y+2;
  
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,0);
  x02_sel = BBE_MOVVSELN_2NX16(vTmp,0);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, vpx, 4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y0, vpy, 4*BBE_SIMD_WIDTH);

    BBE_LSN_2XF32_IP(x10, spx, sizeof(float32_t));
    BBE_LSN_2XF32_IP(x11, spx, sizeof(float32_t));
    BBE_LSN_2XF32_XP(x12, spx, 4*BBE_SIMD_WIDTH-2*sizeof(float32_t));
    
    BBE_LSN_2XF32_IP(y10, spy, 3*sizeof(float32_t));
    BBE_LSN_2XF32_IP(y11, spy, 3*sizeof(float32_t));
    BBE_LSN_2XF32_XP(y12, spy, 4*BBE_SIMD_WIDTH-6*sizeof(float32_t));

    X1 = x12;
    Y1 = y12;

    x00 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_1);
    x01 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_2);
    x02 = BBE_SELN_2XF32(X1, X0, x02_sel);

    y00 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_1);
    y01 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_2);
    y02 = BBE_SELN_2XF32I(Y1, Y0, BBE_SELI_MMC3X3X3X3_OFFSET_M2_0_STEP_3);

    /* Multiply input matrices X and Y */
    Z0 = BBE_MULN_2XF32(x02, y02);
    BBE_MULAN_2XF32(Z0, x01, y01);
    BBE_MULAN_2XF32(Z0, x00, y00);
    
    Z1 = BBE_MULN_2XF32(x10, y10);
    BBE_MULAN_2XF32(Z1, x11, y11);
    BBE_MULAN_2XF32(Z1, x12, y12);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, vpz, 4*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(Z1, spz, 4*BBE_SIMD_WIDTH);
  }
#endif
} /* matmul3x3nf() */
#endif

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
DISCARD_FUN(void, matmul4x4nf,( float32_t * restrict z, 
                          const float32_t * restrict x, 
                          const float32_t * restrict y, 
                          int L ))
#else

#undef BBE_SHFLN_2XF32
#ifndef BBE_SHFLN_2XF32
#define BBE_SHFLN_2XF32(a, b) (BBE_MOVN_2XF32_FROMNX16(BBE_SHFLNX16(BBE_MOVNX16_FROMN_2XF32(a), (b) )))
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

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void matmul4x4nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L )
{
#if 0
  int l;

  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  xb_vecN_2xf32 X0, Y0, X1, Y1, Z0, Z1;
  xb_vecN_2xf32 sel_x00, sel_x01, sel_x023;
  xb_vecN_2xf32 sel_x10, sel_x11, sel_x123;
  xb_vecN_2xf32 sel_y00, sel_y01, sel_y02, sel_y03;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y1, py, 2 * BBE_SIMD_WIDTH);

    sel_x00 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_1);
    sel_x01 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_2);
    sel_x023= BBE_SHFLN_2XF32I(X0, BBE_SHFLI_DUPLICATE_4_ODD);
    sel_x10 = BBE_SHFLN_2XF32I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_1);
    sel_x11 = BBE_SHFLN_2XF32I(X1, BBE_SHFLI_MMC4X4X4X4_M1_STEP_2);
    sel_x123= BBE_SHFLN_2XF32I(X1, BBE_SHFLI_DUPLICATE_4_ODD);
    sel_y00 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
    sel_y01 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
    sel_y02 = BBE_SHFLN_2XF32I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_3);
    sel_y03 = BBE_SHFLN_2XF32I(Y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_4);

    /* Compute matrix Z (first 2 rows) */
    Z0 = BBE_MULN_2XF32(sel_x00, sel_y00);
    BBE_MULAN_2XF32(Z0, sel_x01, sel_y01);
    BBE_MULMASN_2XF32(Z0, sel_x023, sel_y02, 0, 0x4);
    BBE_MULMASN_2XF32(Z0, sel_x023, sel_y03, 0, 0xE);
    /* Compute matrix Z (second 2 rows) */
    Z1 = BBE_MULN_2XF32(sel_x10, sel_y00);
    BBE_MULAN_2XF32(Z1, sel_x11, sel_y01);
    BBE_MULMASN_2XF32(Z1, sel_x123, sel_y02, 0, 0x4);
    BBE_MULMASN_2XF32(Z1, sel_x123, sel_y03, 0, 0xE);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(Z1, pz, 2 * BBE_SIMD_WIDTH);
  }
#elif 0
  #define BBE_SHFLN_2XF32(a, b) (BBE_MOVN_2XF32_FROMNX16(BBE_SHFLNX16(BBE_MOVNX16_FROMN_2XF32(a), (b) )))
  static const uint16_t ALIGN(32) tblSel[2*BBE_SIMD_WIDTH]=
  {
    0x00, 0x01, 0x04, 0x05, 0x00, 0x01, 0x04, 0x05, 0x0A, 0x0B, 0x0E, 0x0F, 0x0A, 0x0B, 0x0E, 0x0F,
    0x02, 0x03, 0x06, 0x07, 0x02, 0x03, 0x06, 0x07, 0x08, 0x09, 0x0C, 0x0D, 0x08, 0x09, 0x0C, 0x0D
  };
  int l;

  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  xb_vecN_2xf32 X0, Y0, X1, Y1, Z0, Z1;
  xb_vecN_2xf32 sel_x00, sel_x01, sel_x10, sel_x11;
  xb_vecN_2xf32 sel_y00, sel_y01, sel_y10, sel_y11;
  xb_vecN_2xf32 z00, z01, z10, z11;
  xb_vecNx16 vTmp;
  vselN sel0, sel1;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;

  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,0*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel0 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,1*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel1 = BBE_MOVVSELNX16(vTmp,0);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y1, py, 2 * BBE_SIMD_WIDTH);

    sel_x00 = BBE_SHFLN_2XF32(X0, sel0);
    sel_x01 = BBE_SHFLN_2XF32(X0, sel1);
    sel_x10 = BBE_SHFLN_2XF32(X1, sel0);
    sel_x11 = BBE_SHFLN_2XF32(X1, sel1);
    sel_y00 = Y0;
    sel_y01 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_SWAP_8);
    sel_y10 = Y1;
    sel_y11 = BBE_SHFLN_2XF32I(Y1, BBE_SHFLI_SWAP_8);

    /* Compute matrix Z (first 2 rows) */
    z00 = BBE_MULMN_2XF32( sel_x00, sel_y00, 0, 0x4);
    z01 = BBE_MULMN_2XF32( sel_x00, sel_y10, 0, 0xE);
    BBE_MULMASN_2XF32(z00, sel_x01, sel_y01, 0, 0x4);
    BBE_MULMASN_2XF32(z01, sel_x01, sel_y11, 0, 0xE);
    Z0 = BBE_ADDN_2XF32(z00, z01);

    /* Compute matrix Z (second 2 rows) */
    z10 = BBE_MULMN_2XF32( sel_x10, sel_y00, 0, 0x4);
    z11 = BBE_MULMN_2XF32( sel_x10, sel_y10, 0, 0xE);
    BBE_MULMASN_2XF32(z10, sel_x11, sel_y01, 0, 0x4);
    BBE_MULMASN_2XF32(z11, sel_x11, sel_y11, 0, 0xE);
    Z1 = BBE_ADDN_2XF32(z10, z11);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(Z1, pz, 2 * BBE_SIMD_WIDTH);
  }
#else
  static const uint16_t ALIGN(32) tblSel[2*BBE_SIMD_WIDTH]=
  {
    0x00, 0x01, 0x04, 0x05, 0x00, 0x01, 0x04, 0x05, 0x0A, 0x0B, 0x0E, 0x0F, 0x0A, 0x0B, 0x0E, 0x0F,
    0x02, 0x03, 0x06, 0x07, 0x02, 0x03, 0x06, 0x07, 0x08, 0x09, 0x0C, 0x0D, 0x08, 0x09, 0x0C, 0x0D
  };
  int l;

  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict py0;
  const xb_vecN_2xf32 * restrict py1;
        xb_vecN_2xf32 * restrict pz0;
        xb_vecN_2xf32 * restrict pz1;
  xb_vecN_2xf32 X0, Y0, X1, Y1, Z0, Z1;
  xb_vecN_2xf32 sel_x00, sel_x01, sel_x10, sel_x11;
  xb_vecN_2xf32 sel_y00, sel_y01, sel_y10, sel_y11;
  xb_vecNx16 vTmp;
  vselN sel0, sel1;

  /* check restrictions */
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px0 = (const xb_vecN_2xf32 *)(x);
  px1 = (const xb_vecN_2xf32 *)(x+8);
  py0 = (const xb_vecN_2xf32 *)(y);
  py1 = (const xb_vecN_2xf32 *)(y+8);
  pz0 = (      xb_vecN_2xf32 *)(z);
  pz1 = (      xb_vecN_2xf32 *)(z+8);

  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,0*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel0 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,1*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel1 = BBE_MOVVSELNX16(vTmp,0);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, px0, 4 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X1, px1, 4 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y0, py0, 4 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y1, py1, 4 * BBE_SIMD_WIDTH);

    sel_x00 = BBE_SHFLN_2XF32(X0, sel0);
    sel_x01 = BBE_SHFLN_2XF32(X0, sel1);
    sel_x10 = BBE_SHFLN_2XF32(X1, sel0);
    sel_x11 = BBE_SHFLN_2XF32(X1, sel1);
    sel_y00 = Y0;
    sel_y01 = BBE_SHFLN_2XF32I(Y0, BBE_SHFLI_SWAP_8);
    sel_y10 = Y1;
    sel_y11 = BBE_SHFLN_2XF32I(Y1, BBE_SHFLI_SWAP_8);

    /* Compute matrix Z (first 2 rows) */
    Z0 = BBE_MULMN_2XF32( sel_x00, sel_y00, 0, 0x4);
    BBE_MULMASN_2XF32(Z0, sel_x00, sel_y10, 0, 0xE);
    BBE_MULMASN_2XF32(Z0, sel_x01, sel_y01, 0, 0x4);
    BBE_MULMASN_2XF32(Z0, sel_x01, sel_y11, 0, 0xE);

    /* Compute matrix Z (second 2 rows) */
    Z1 = BBE_MULMN_2XF32( sel_x10, sel_y00, 0, 0x4);
    BBE_MULMASN_2XF32(Z1, sel_x10, sel_y10, 0, 0xE);
    BBE_MULMASN_2XF32(Z1, sel_x11, sel_y01, 0, 0x4);
    BBE_MULMASN_2XF32(Z1, sel_x11, sel_y11, 0, 0xE);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, pz0, 4 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(Z1, pz1, 4 * BBE_SIMD_WIDTH);
  }
#endif
} /* matmul4x4nf() */
#endif

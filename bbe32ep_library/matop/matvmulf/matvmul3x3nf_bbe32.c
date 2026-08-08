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
DISCARD_FUN(void, matvmul3x3nf,( float32_t * restrict z, 
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

/* Block Order, Floating-Point, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 2
*/
void matvmul3x3nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L )
{
  static const uint16_t ALIGN(32) tblSel[3*BBE_SIMD_WIDTH]=
  {
    0x00, 0x01, 0x06, 0x07, 0x0C, 0x0D, 0x0C, 0x0D, 0x10, 0x11, 0x16, 0x17, 0x1C, 0x1D, 0x1C, 0x1D,
    0x02, 0x03, 0x08, 0x09, 0x0E, 0x0F, 0x0E, 0x0F, 0x12, 0x13, 0x18, 0x19, 0x1E, 0x1F, 0x1E, 0x1F,
    0x06, 0x07, 0x0C, 0x0D, 0x00, 0x01, 0x00, 0x01, 0x16, 0x17, 0x1C, 0x1D, 0x10, 0x11, 0x10, 0x11
  };
  const xb_vecN_2xf32 * restrict px0;
  const float32_t     * restrict px1;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l;

  xb_vecN_2xf32 X00, X01, X10, X11, Y, Z;
  xb_vecN_2xf32 x0, x1, x2, y01, y2;
  xb_vecNx16 vTmp;
  vselN_2 sel0, sel1, sel2;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&1) == 0);

  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,0*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel0 = BBE_MOVVSELN_2NX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,1*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel1 = BBE_MOVVSELN_2NX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,2*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel2 = BBE_MOVVSELN_2NX16(vTmp,0);

  px0 = (const xb_vecN_2xf32 *)x;
  px1 = (const float32_t     *)x+8;
  py  = (const xb_vecN_2xf32 *)y;
  pz  = (      xb_vecN_2xf32 *)z;

  for (l = 0; l<(L>>1); l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X00, px0, 4*BBE_SIMD_WIDTH);
    BBE_LSN_2XF32_XP(X01, px1, 4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X10, px0, 4*BBE_SIMD_WIDTH);
    BBE_LSN_2XF32_XP(X11, px1, 4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y , py, 2*BBE_SIMD_WIDTH);

    x0 = BBE_SELN_2XF32(X10, X00, sel0);
    x1 = BBE_SELN_2XF32(X10, X00, sel1);
    X00 = BBE_SELN_2XF32I(X00, X01, BBE_SELI_PACK_2);
    X10 = BBE_SELN_2XF32I(X10, X11, BBE_SELI_PACK_2);
    x2 = BBE_SELN_2XF32(X10, X00, sel2);
    y01 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_DUPLICATE_4_EVEN);
    y2  = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC4X4X4X4_M1_STEP_3);

    /* Multiply input matrices X and Y */
    Z = BBE_MULMN_2XF32( x0, y01, 0, 0x8);
    BBE_MULMASN_2XF32(Z, x1, y01, 0, 0xD);
    BBE_MULAN_2XF32(Z, x2, y2);

    /* Save results */
    BBE_SVN_2XF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
} /* matvmul3x3nf() */
#endif

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
DISCARD_FUN(void, cmatvmul3x3nf,( complex_float * restrict z, 
                            const complex_float * restrict x, 
                            const complex_float * restrict y, 
                            int L ))
#else

#ifndef BBE_SELN_4XCF32
#define BBE_SELN_4XCF32(a, b, c) (BBE_MOVN_4XCF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_4XCF32(a), BBE_MOVNX16_FROMN_4XCF32(b), (c) )))
#endif

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

/* Block Order, Floating-Point, 3x3*3x1->3x1, Sx=12, Sy=4, Sz=4
   Restrictions:
     None
*/
void cmatvmul3x3nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L )
{
#if 0
  static const uint16_t ALIGN(32) tblSel[3*BBE_SIMD_WIDTH]=
  {
    0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B, 0x18, 0x19, 0x1A, 0x1B,
    0x04, 0x05, 0x06, 0x07, 0x10, 0x11, 0x12, 0x13, 0x1C, 0x1D, 0x1E, 0x1F, 0x1C, 0x1D, 0x1E, 0x1F,
    0x08, 0x09, 0x0A, 0x0B, 0x14, 0x15, 0x16, 0x17, 0x08, 0x09, 0x0A, 0x0B, 0x14, 0x15, 0x16, 0x17
  };
  const xb_vecN_4xcf32 * restrict px0;
  const complex_float  * restrict px1;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, Y, Z;
  xb_vecN_4xcf32 x0, x1, x2, y0, y1, y2, z0, z1;
  xb_vecNx16 vTmp;
  vselN sel0, sel1, sel2;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,0*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel0 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,1*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel1 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,2*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel2 = BBE_MOVVSELNX16(vTmp,0);

  px0 = (const xb_vecN_4xcf32 *)x;
  px1 = (const complex_float  *)x+8;
  py  = (const xb_vecN_4xcf32 *)y;
  pz  = (      xb_vecN_4xcf32 *)z;

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, px0, 4*BBE_SIMD_WIDTH);
    BBE_LSN_4XCF32_XP(X2, px1, 6*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y , py, 2*BBE_SIMD_WIDTH);

    x0 = BBE_SELN_4XCF32(X1, X0, sel0);
    x1 = BBE_SELN_4XCF32(X1, X0, sel1);
    x2 = BBE_SELN_4XCF32(X1, X0, sel2);
    x2 = BBE_SELN_4XCF32I(X2, x2, BBE_SELI_PACK_8);
    y0 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_REP_0X4);
    y1 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_REP_1X4);
    y2 = BBE_SHFLN_4XCF32I(Y, BBE_SHFLI_REP_2X4);

    /* Multiply input matrices X and Y */
    z0 = BBE_MULMN_4XCF32( x0, y0, 0, 0x4);
    z1 = BBE_MULMN_4XCF32( x0, y0, 1, 0xB);
    BBE_MULMASN_4XCF32(z0, x1, y1, 0, 0x4);
    BBE_MULMASN_4XCF32(z1, x1, y1, 1, 0xB);
    BBE_MULMASN_4XCF32(z0, x2, y2, 0, 0x4);
    BBE_MULMASN_4XCF32(z1, x2, y2, 1, 0xB);
    Z = BBE_ADDN_4XCF32(z0, z1);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_4xcf32 * restrict px;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, Xskip;
  xb_vecN_4xcf32 Y;
  xb_vecN_4xcf32 Z, Z0, Z1, Z2, Zt;
  valign alx;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_4xcf32 *)x;
  py = (const xb_vecN_4xcf32 *)y;
  pz = (      xb_vecN_4xcf32 *)z;
  alx = BBE_LAN_4XCF32_PP(px);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LAVN_4XCF32_XP(X0,    alx, px, 3*sizeof(complex_float));
    BBE_LAVN_4XCF32_XP(X1,    alx, px, 3*sizeof(complex_float));
    BBE_LAVN_4XCF32_XP(X2,    alx, px, 3*sizeof(complex_float));
    BBE_LAVN_4XCF32_XP(Xskip, alx, px, 3*sizeof(complex_float));
    BBE_LVN_4XCF32_IP(Y, py, 2*BBE_SIMD_WIDTH);

    /* Multiply input matrices X and Y */
    Z0 = BBE_MULN_4XCF32(X0, Y);
    Z1 = BBE_MULN_4XCF32(X1, Y);
    Z2 = BBE_MULN_4XCF32(X2, Y);

    BBE_DSELN_4XCF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Zt, Z1, Zt, Z1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Z1, Z0, Z1, Z0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(Zt, Z2, Zt, Z2, BBE_DSELI_INTERLEAVE_4);

    Z = BBE_ADDN_4XCF32(Z0, Z1);
    Z = BBE_ADDN_4XCF32(Z, Z2);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* cmatvmul3x3nf() */
#endif

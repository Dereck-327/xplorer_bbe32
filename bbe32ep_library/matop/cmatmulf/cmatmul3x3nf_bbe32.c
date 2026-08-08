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
DISCARD_FUN(void, cmatmul3x3nf,( complex_float * restrict z, 
                           const complex_float * restrict x, 
                           const complex_float * restrict y,
                           int L ))
#else

#ifndef BBE_SELN_4XCF32
#define BBE_SELN_4XCF32(a, b, c) (BBE_MOVN_4XCF32_FROMNX16(BBE_SELNX16( BBE_MOVNX16_FROMN_4XCF32(a), BBE_MOVNX16_FROMN_4XCF32(b), (c) )))
#endif
#ifndef BBE_SHFLN_4XCF32
#define BBE_SHFLN_4XCF32(a, b) (BBE_MOVN_4XCF32_FROMNX16(BBE_SHFLNX16( BBE_MOVNX16_FROMN_4XCF32(a), (b) )))
#endif
#ifndef BBE_SELUN_4XCF32
#define BBE_SELUN_4XCF32(a, b, c, d, e) \
{ \
  xb_vecNx16 v_tmp; \
  v_tmp = BBE_MOVNX16_FROMN_4XCF32(a); \
  BBE_SELUNX16( v_tmp, BBE_MOVNX16_FROMN_4XCF32(b), BBE_MOVNX16_FROMN_4XCF32(c), (d), (e) ); \
  a = BBE_MOVN_4XCF32_FROMNX16(v_tmp); \
}
#endif
#ifndef BBE_SHFLUN_4XCF32
#define BBE_SHFLUN_4XCF32(a, b, c, d) \
{ \
  xb_vecNx16 v_tmp; \
  v_tmp = BBE_MOVNX16_FROMN_4XCF32(a); \
  BBE_SHFLUNX16( v_tmp, BBE_MOVNX16_FROMN_4XCF32(b), (c), (d) ); \
  a = BBE_MOVN_4XCF32_FROMNX16(v_tmp); \
}
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

/* Block Order, Floating-Point, 3x3*3x3->3x3, Sx=12, Sy=12, Sz=12
   Restrictions:
     None
*/
void cmatmul3x3nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y,
              int L )
{
  static const uint16_t ALIGN(32) aSel[4*BBE_SIMD_WIDTH]=
  {
      0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F,
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x00, 0x01, 0x02, 0x03,
      0x0C, 0x0D, 0x0E, 0x0F, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B, 0x18, 0x19, 0x1A, 0x1B,
      0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07
  };

  const xb_vecN_4xcf32 * restrict vpx;
  const xb_vecN_4xcf32 * restrict vpy;
        xb_vecN_4xcf32 * restrict vpz;
        xb_vecN_4xcf32 * restrict spz;
  const complex_float * restrict spx;
  const complex_float * restrict spy;
  int l;

  xb_vecN_4xcf32 X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2;
  xb_vecN_4xcf32 x00, x01, x02, x10, x11, x12, x20, x21, x22;
  xb_vecN_4xcf32 y00, y01, y02, y10, y11, y12, y20, y21, y22;
  xb_vecN_4xcf32 z00, z01, z10, z11, z20, z21;
  xb_vecNx16 vTmp;
  vselN sel_x0, sel_x1, sel_y0, sel_y1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  vpx = (const xb_vecN_4xcf32 *)(x);
  vpy = (const xb_vecN_4xcf32 *)(y);
  vpz = (      xb_vecN_4xcf32 *)(z);
  spz = (      xb_vecN_4xcf32 *)(z+8);
  spx = x+6;
  spy = y+2;
  
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,0*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel_x0 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,1*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel_y0 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,2*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel_x1 = BBE_MOVVSELNX16(vTmp,0);
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)aSel,3*BBE_SIMD_WIDTH*sizeof(int16_t));
  sel_y1 = BBE_MOVVSELNX16(vTmp,0);

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_4XCF32_IP(X0, vpx, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(X1, vpx, 4*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y0, vpy, 2*BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(Y1, vpy, 4*BBE_SIMD_WIDTH);

    BBE_LSN_4XCF32_IP(x20, spx, sizeof(complex_float));
    BBE_LSN_4XCF32_IP(x21, spx, sizeof(complex_float));
    BBE_LSN_4XCF32_XP(x22, spx, 6*BBE_SIMD_WIDTH-2*sizeof(complex_float));
    
    BBE_LSN_4XCF32_IP(y20, spy, 3*sizeof(complex_float));
    BBE_LSN_4XCF32_IP(y21, spy, 3*sizeof(complex_float));
    BBE_LSN_4XCF32_XP(y22, spy, 6*BBE_SIMD_WIDTH-6*sizeof(complex_float));

    X2 = x22;
    Y2 = y22;

    BBE_SHFLUN_4XCF32(x00,     X0, sel_x0, 4);
    BBE_SELUN_4XCF32 (x01, X1, X0, sel_x0, 4);
    BBE_SELUN_4XCF32 (x02, X1, X0, sel_x0, -8);

    BBE_SELUN_4XCF32 (x10, X1, X0, sel_x1, -12);
    BBE_SHFLUN_4XCF32(x11,     X1, sel_x1, 4);
    BBE_SELUN_4XCF32 (x12, X2, X1, sel_x1, 8);

    BBE_SHFLUN_4XCF32(y00,     Y0, sel_y0, 12);
    BBE_SELUN_4XCF32 (y01, Y1, Y0, sel_y0, -4);
    BBE_SELUN_4XCF32 (y02, Y2, Y1, sel_y0, -8);

    BBE_SHFLUN_4XCF32(y10,     Y0, sel_y1, 12);
    BBE_SELUN_4XCF32 (y11, Y1, Y0, sel_y1, -4);
    BBE_SELUN_4XCF32 (y12, Y2, Y1, sel_y1, -8);

    /* Multiply input matrices X and Y */
    z00 = BBE_MULMN_4XCF32( x00, y00, 0, 0x4);
    z01 = BBE_MULMN_4XCF32( x00, y00, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x01, y01, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x01, y01, 1, 0xB);
    BBE_MULMASN_4XCF32(z00, x02, y02, 0, 0x4);
    BBE_MULMASN_4XCF32(z01, x02, y02, 1, 0xB);
    Z0 = BBE_ADDN_4XCF32(z00, z01);

    z10 = BBE_MULMN_4XCF32( x10, y10, 0, 0x4);
    z11 = BBE_MULMN_4XCF32( x10, y10, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x11, y11, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x11, y11, 1, 0xB);
    BBE_MULMASN_4XCF32(z10, x12, y12, 0, 0x4);
    BBE_MULMASN_4XCF32(z11, x12, y12, 1, 0xB);
    Z1 = BBE_ADDN_4XCF32(z10, z11);

    z20 = BBE_MULMN_4XCF32( x20, y20, 0, 0x4);
    z21 = BBE_MULMN_4XCF32( x20, y20, 1, 0xB);
    BBE_MULMASN_4XCF32(z20, x21, y21, 0, 0x4);
    BBE_MULMASN_4XCF32(z21, x21, y21, 1, 0xB);
    BBE_MULMASN_4XCF32(z20, x22, y22, 0, 0x4);
    BBE_MULMASN_4XCF32(z21, x22, y22, 1, 0xB);
    Z2 = BBE_ADDN_4XCF32(z20, z21);

    /* Save results */
    BBE_SVN_4XCF32_IP(Z0, vpz, 2*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z1, vpz, 4*BBE_SIMD_WIDTH);
    BBE_SVN_4XCF32_IP(Z2, spz, 6*BBE_SIMD_WIDTH);
  }
} /* cmatmul3x3nf() */
#endif

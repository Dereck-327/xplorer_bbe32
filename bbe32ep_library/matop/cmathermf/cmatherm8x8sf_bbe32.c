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
DISCARD_FUN(void, cmatherm8x8sf,( complex_float * restrict y, 
                            const complex_float * restrict x, 
                            int L ))
#else

#define sz_cf32 ((int)sizeof(complex_float))
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

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm8x8sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
#if 0
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict pxr;
        xb_vecN_4xcf32 * restrict py0;
        xb_vecN_4xcf32 * restrict py1;
  int stridex;
  int k, modinc;
  int l, n;

  xb_vecN_4xcf32 x00, x10, x20, x30, x40, x50, x60, x70,
                 x01, x11, x21, x31, x41, x51, x61, x71;
  xb_vecN_4xcf32 xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7;
  xb_vecN_4xcf32 y00, y01, y10, y11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  /*
   * Compute matrices by 2 rows per iteration
   */
  for ( n=0; n<8; n+=2 )
  {
    px0 = (const xb_vecN_4xcf32 *)(x+n*L);
    px1 = (const xb_vecN_4xcf32 *)(x+n*L+L);
    pxr = (const xb_vecN_4xcf32 *)(x);
    py0 = (      xb_vecN_4xcf32 *)(y+n*8*L);
    py1 = (      xb_vecN_4xcf32 *)(y+n*8*L+8*L);

    __Pragma("loop_count min=4, factor=4");
    for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*8; l++)
    {
      k = BBE_ADDMOD16U(k, modinc);
      stridex = -56*L*sz_cf32+2*BBE_SIMD_WIDTH;
      XT_MOVEQZ(stridex, -56*L*sz_cf32+2*BBE_SIMD_WIDTH-L*sz_cf32, k);

      /* Load 2 columns of input matrix X */
      BBE_LVN_4XCF32_XP(x00, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x10, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x20, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x30, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x40, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x50, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x60, px0, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x70, px0, stridex);
      BBE_LVN_4XCF32_XP(x01, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x11, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x21, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x31, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x41, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x51, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x61, px1, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(x71, px1, stridex);

      /* Load 1 column of input matrix X */
      BBE_LVN_4XCF32_XP(xr0, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr1, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr2, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr3, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr4, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr5, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr6, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr7, pxr, -56*L*sz_cf32+2*BBE_SIMD_WIDTH);

      /* Compute 2 values for 2 rows of matrix Y */
      y00 = BBE_MULMN_4XCF32( x00, xr0, 0, 0x4);
      y01 = BBE_MULMN_4XCF32( x00, xr0, 2, 0xB);
      y10 = BBE_MULMN_4XCF32( x01, xr0, 0, 0x4);
      y11 = BBE_MULMN_4XCF32( x01, xr0, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x10, xr1, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x10, xr1, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x11, xr1, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x11, xr1, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x20, xr2, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x20, xr2, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x21, xr2, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x21, xr2, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x30, xr3, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x30, xr3, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x31, xr3, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x31, xr3, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x40, xr4, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x40, xr4, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x41, xr4, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x41, xr4, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x50, xr5, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x50, xr5, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x51, xr5, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x51, xr5, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x60, xr6, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x60, xr6, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x61, xr6, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x61, xr6, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, x70, xr7, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, x70, xr7, 2, 0xB);
      BBE_MULMASN_4XCF32(y10, x71, xr7, 0, 0x4);
      BBE_MULMASN_4XCF32(y11, x71, xr7, 2, 0xB);
      /* Make dummy operations for better loop scheduling */
      y00 = BBE_ADDN_4XCF32(y00, 0.0f);
      y10 = BBE_ADDN_4XCF32(y10, 0.0f);

      y00 = BBE_ADDN_4XCF32(y00, y01);
      y10 = BBE_ADDN_4XCF32(y10, y11);

      /* Save results (by one values for 2 rows) */
      BBE_SVN_4XCF32_IP(y00, py0, 2*BBE_SIMD_WIDTH);
      BBE_SVN_4XCF32_IP(y10, py1, 2*BBE_SIMD_WIDTH);
    }
  }
#else
  const xb_vecN_4xcf32 * restrict pxl;
  const xb_vecN_4xcf32 * restrict pxr;
        xb_vecN_4xcf32 * restrict pyrow;
        xb_vecN_4xcf32 * restrict pycol;
  int stridex, stridey;
  int k, modinc;
  int l, n;

  xb_vecN_4xcf32 xl0, xl1, xl2, xl3, xl4, xl5, xl6, xl7;
  xb_vecN_4xcf32 xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7;
  xb_vecN_4xcf32 y00, y01, y02, y03, y00conj;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  /*
   * Compute matrices by 2 rows per iteration
   */
  for ( n=0; n<8; n+=1 )
  {
    pxl = pxr = (const xb_vecN_4xcf32 *)(x+n*L);
    pyrow = pycol = (  xb_vecN_4xcf32 *)(y+n*9*L);

    __Pragma("loop_count min=1");
    for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*(8-n); l++)
    {
      k = BBE_ADDMOD16U(k, modinc);
      stridex = -56*L*sz_cf32+2*BBE_SIMD_WIDTH;
      XT_MOVEQZ(stridex, -56*L*sz_cf32+2*BBE_SIMD_WIDTH-L*sz_cf32, k);
      stridey = 2*BBE_SIMD_WIDTH;
      XT_MOVEQZ(stridey, 7*L*sz_cf32+2*BBE_SIMD_WIDTH, k);

      /* Load 1 column of input matrix X */
      BBE_LVN_4XCF32_XP(xl0, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl1, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl2, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl3, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl4, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl5, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl6, pxl, 8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xl7, pxl, stridex);

      /* Load 1 column of input matrix X */
      BBE_LVN_4XCF32_XP(xr0, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr1, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr2, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr3, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr4, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr5, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr6, pxr,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(xr7, pxr, -56*L*sz_cf32+2*BBE_SIMD_WIDTH);

      /* Compute 2 values for 2 rows of matrix Y */
      y00 = BBE_MULMN_4XCF32( xr0, xl0, 0, 0x4);
      y01 = BBE_MULMN_4XCF32( xr0, xl0, 2, 0xB);
      y02 = BBE_MULMN_4XCF32( xr1, xl1, 0, 0x4);
      y03 = BBE_MULMN_4XCF32( xr1, xl1, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, xr2, xl2, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, xr2, xl2, 2, 0xB);
      BBE_MULMASN_4XCF32(y02, xr3, xl3, 0, 0x4);
      BBE_MULMASN_4XCF32(y03, xr3, xl3, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, xr4, xl4, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, xr4, xl4, 2, 0xB);
      BBE_MULMASN_4XCF32(y02, xr5, xl5, 0, 0x4);
      BBE_MULMASN_4XCF32(y03, xr5, xl5, 2, 0xB);
      BBE_MULMASN_4XCF32(y00, xr6, xl6, 0, 0x4);
      BBE_MULMASN_4XCF32(y01, xr6, xl6, 2, 0xB);
      BBE_MULMASN_4XCF32(y02, xr7, xl7, 0, 0x4);
      BBE_MULMASN_4XCF32(y03, xr7, xl7, 2, 0xB);
      /* Make dummy operations for better loop scheduling */
      y00 = BBE_ADDN_4XCF32(y00, 0.0f);
      y02 = BBE_ADDN_4XCF32(y02, 0.0f);

      y00 = BBE_ADDN_4XCF32(y00, y01);
      y02 = BBE_ADDN_4XCF32(y02, y03);
      y00conj = BBE_ADDN_4XCF32(y00, y02);
      y00 = BBE_CONJN_4XCF32(y00conj);

      /* Save results (by one values for 1 row and 1 column) */
      BBE_SVN_4XCF32_IP(y00    , pyrow, 2*BBE_SIMD_WIDTH);
      BBE_SVN_4XCF32_XP(y00conj, pycol, stridey);
    }
  }
#endif
} /* cmatherm8x8sf() */
#endif

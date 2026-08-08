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
DISCARD_FUN(void, cmatmul8x8sf, ( complex_float * restrict z, 
                            const complex_float * restrict x, 
                            const complex_float * restrict y, 
                            int L ))
#else

#define sz_cf32 ((int)sizeof(complex_float))
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

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul8x8sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz0;
        xb_vecN_4xcf32 * restrict pz1;
  int stridex;
  int k, modinc;
  int l, n;

  xb_vecN_4xcf32 x00, x01, x02, x03, x04, x05, x06, x07,
                 x10, x11, x12, x13, x14, x15, x16, x17;
  xb_vecN_4xcf32 y00, y10, y20, y30, y40, y50, y60, y70;
  xb_vecN_4xcf32 z00, z01, z10, z11;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;

  k = 0;
  modinc = (L<<16) | (BBE_SIMD_WIDTH/4);

  /*
   * Compute matrices by 2 rows per iteration
   */
  for ( n=0; n<8; n+=2 )
  {
    px0 = (const xb_vecN_4xcf32 *)(x+n*8*L);
    px1 = (const xb_vecN_4xcf32 *)(x+n*8*L+8*L);
    py  = (const xb_vecN_4xcf32 *)(y);
    pz0 = (      xb_vecN_4xcf32 *)(z+n*8*L);
    pz1 = (      xb_vecN_4xcf32 *)(z+n*8*L+8*L);

    __Pragma("loop_count min=4, factor=4");
    for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2))*8; l++)
    {
      k = BBE_ADDMOD16U(k, modinc);
      stridex = -7*L*sz_cf32+2*BBE_SIMD_WIDTH;
      XT_MOVEQZ(stridex, -7*L*sz_cf32+2*BBE_SIMD_WIDTH-L*sz_cf32, k);

      /* Load 2 rows of input matrix X */
      BBE_LVN_4XCF32_XP(x00, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x01, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x02, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x03, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x04, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x05, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x06, px0, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x07, px0, stridex);
      BBE_LVN_4XCF32_XP(x10, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x11, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x12, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x13, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x14, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x15, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x16, px1, L*sz_cf32);
      BBE_LVN_4XCF32_XP(x17, px1, stridex);

      /* Load 1 column of input matrix Y */
      BBE_LVN_4XCF32_XP(y00, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y10, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y20, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y30, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y40, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y50, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y60, py,   8*L*sz_cf32);
      BBE_LVN_4XCF32_XP(y70, py, -56*L*sz_cf32+2*BBE_SIMD_WIDTH);

      /* Compute 2 values for 2 rows of matrix Z */
      z00 = BBE_MULMN_4XCF32( x00, y00, 0, 0x4);
      z01 = BBE_MULMN_4XCF32( x00, y00, 1, 0xB);
      z10 = BBE_MULMN_4XCF32( x10, y00, 0, 0x4);
      z11 = BBE_MULMN_4XCF32( x10, y00, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x01, y10, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x01, y10, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x11, y10, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x11, y10, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x02, y20, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x02, y20, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x12, y20, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x12, y20, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x03, y30, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x03, y30, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x13, y30, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x13, y30, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x04, y40, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x04, y40, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x14, y40, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x14, y40, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x05, y50, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x05, y50, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x15, y50, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x15, y50, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x06, y60, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x06, y60, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x16, y60, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x16, y60, 1, 0xB);
      BBE_MULMASN_4XCF32(z00, x07, y70, 0, 0x4);
      BBE_MULMASN_4XCF32(z01, x07, y70, 1, 0xB);
      BBE_MULMASN_4XCF32(z10, x17, y70, 0, 0x4);
      BBE_MULMASN_4XCF32(z11, x17, y70, 1, 0xB);
      /* Make dummy operations for better loop scheduling */
      z00 = BBE_ADDN_4XCF32(z00, 0.0f);
      z10 = BBE_ADDN_4XCF32(z10, 0.0f);

      z00 = BBE_ADDN_4XCF32(z00, z01);
      z10 = BBE_ADDN_4XCF32(z10, z11);

      /* Save results (by one values for 2 rows) */
      BBE_SVN_4XCF32_IP(z00, pz0, 2*BBE_SIMD_WIDTH);
      BBE_SVN_4XCF32_IP(z10, pz1, 2*BBE_SIMD_WIDTH);
    }
  }
} /* cmatmul8x8sf() */
#endif

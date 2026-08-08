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
DISCARD_FUN(void, matvmul4x4nf,( float32_t * restrict z, 
                           const float32_t * restrict x, 
                           const float32_t * restrict y, 
                           int L ))
#else
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

/* Block Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 2
*/
void matvmul4x4nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L )
{
#if 1
  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l;

  xb_vecN_2xf32 X00, X01, X10, X11, Y, Z;
  xb_vecN_2xf32 x0, x1, x2, x3, y01, y23, z0, z1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&1) == 0);

  px0 = (const xb_vecN_2xf32 *)x;
  px1 = (const xb_vecN_2xf32 *)x + 2;
  py  = (const xb_vecN_2xf32 *)y;
  pz  = (      xb_vecN_2xf32 *)z;

  for (l = 0; l<(L>>1); l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_XP(X00, px0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X01, px0, 6*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X10, px1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X11, px1, 6*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(Y  , py , 2*BBE_SIMD_WIDTH);
    
    BBE_DSELN_2XF32I(X01, X00, X01, X00, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X11, X10, X11, X10, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(x2 , x0 , X10, X00, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(x3 , x1 , X11, X01, BBE_DSELI_DEINTERLEAVE_2);
    y01 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_DUPLICATE_4_EVEN);
    y23 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_DUPLICATE_4_ODD);

    /* Multiply input matrices X and Y */
    z0 = BBE_MULMN_2XF32( x0, y01, 0, 0x8);
    BBE_MULMASN_2XF32(z0, x2, y23, 0, 0x8);
    z1 = BBE_MULMN_2XF32( x1, y01, 0, 0xD);
    BBE_MULMASN_2XF32(z1, x3, y23, 0, 0xD);

    /* Save results */
    Z = BBE_ADDN_2XF32(z0, z1);
    BBE_SVN_2XF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l;

  xb_vecN_2xf32 X00, X01, X10, X11, Y, Z;
  xb_vecN_2xf32 x0, x1, x2, x3, y01, y23, z0, z1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&1) == 0);

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;

  for (l = 0; l<(L>>1); l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X00, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X01, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X10, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X11, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y  , py, 2*BBE_SIMD_WIDTH);
    
    BBE_DSELN_2XF32I(X01, X00, X01, X00, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X11, X10, X11, X10, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(x2 , x0 , X10, X00, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(x3 , x1 , X11, X01, BBE_DSELI_DEINTERLEAVE_2);
    y01 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_DUPLICATE_4_EVEN);
    y23 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_DUPLICATE_4_ODD);
    BBE_DSELN_2XF32I(y23, y01, Y  , Y  , BBE_DSELI_DEINTERLEAVE_2);

    /* Multiply input matrices X and Y */
    z0 = BBE_MULMN_2XF32( x0, y01, 0, 0x8);
    BBE_MULMASN_2XF32(z0, x1, y01, 0, 0xD);
    z1 = BBE_MULMN_2XF32( x2, y23, 0, 0x8);
    BBE_MULMASN_2XF32(z1, x3, y23, 0, 0xD);

    /* Save results */
    Z = BBE_ADDN_2XF32(z0, z1);
    BBE_SVN_2XF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* matvmul4x4nf() */
#endif

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
DISCARD_FUN(void, matvmul8x8nf,( float32_t * restrict z, 
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

/* Block Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
      None
*/
void matvmul8x8nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L )
{
#if 0
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l;

  xb_vecN_2xf32 X0, X1, X2, X3, X4, X5, X6, X7, Y, Z;
  xb_vecN_2xf32 y01, y23, y45, y67, z0, z1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;

  for (l=0; l<L; l++)
  {
    /* Load input matrices X and Y */
    BBE_LVN_2XF32_IP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X2, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X3, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X4, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X5, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X6, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(X7, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(Y , py, 2*BBE_SIMD_WIDTH);
    
    BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X5, X4, X5, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X7, X6, X7, X6, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X6, X4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X7, X5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X4, X0, X4, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X5, X1, X5, X1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X6, X2, X6, X2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(X7, X3, X7, X3, BBE_DSELI_DEINTERLEAVE_2);
    y01 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_0X4);
    y23 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_1X4);
    y45 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_2X4);
    y67 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_3X4);

    /* Multiply input matrices X and Y */
    z0 = BBE_MULMN_2XF32( X0, y01, 0, 0x8);
    BBE_MULMASN_2XF32(z0, X1, y01, 0, 0xD);
    BBE_MULMASN_2XF32(z0, X2, y23, 0, 0x8);
    BBE_MULMASN_2XF32(z0, X3, y23, 0, 0xD);
    z1 = BBE_MULMN_2XF32( X4, y45, 0, 0x8);
    BBE_MULMASN_2XF32(z1, X5, y45, 0, 0xD);
    BBE_MULMASN_2XF32(z1, X6, y67, 0, 0x8);
    BBE_MULMASN_2XF32(z1, X7, y67, 0, 0xD);

    /* Save results */
    Z = BBE_ADDN_2XF32(z0, z1);
    BBE_SVN_2XF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l;

  xb_vecN_2xf32 X0, X1, X2, X3, X4, X5, X6, X7, Y;
  xb_vecN_2xf32 Z0, Z2, Z4, Z6;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;

  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_2XF32_XP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X2, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X3, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X4, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X5, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X6, px, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(X7, px, 2*BBE_SIMD_WIDTH);
    /* Load vector Y */
    BBE_LVN_2XF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
    Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

    BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_2);

    Z0 = BBE_MULMN_2XF32( X0, Y, 0, 0x8);
    BBE_MULMASN_2XF32(Z0, X1, Y, 0, 0xD);
    Z2 = BBE_MULMN_2XF32( X2, Y, 0, 0x8);
    BBE_MULMASN_2XF32(Z2, X3, Y, 0, 0xD);
    Z4 = BBE_MULMN_2XF32( X4, Y, 0, 0x8);
    BBE_MULMASN_2XF32(Z4, X5, Y, 0, 0xD);
    Z6 = BBE_MULMN_2XF32( X6, Y, 0, 0x8);
    BBE_MULMASN_2XF32(Z6, X7, Y, 0, 0xD);

    BBE_DSELN_2XF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_2XF32I(Z6, Z4, Z6, Z4, BBE_DSELI_INTERLEAVE_4);
    Z0 = BBE_ADDN_2XF32(Z0, Z2);
    Z4 = BBE_ADDN_2XF32(Z4, Z6);

    Z2 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_LO_HALVES);
    Z6 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_HI_HALVES);
    Z0 = BBE_ADDN_2XF32(Z2, Z6);

    /* Save results */
    BBE_SVN_2XF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);
  }
#endif
} /* matvmul8x8nf() */
#endif

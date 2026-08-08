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
DISCARD_FUN(void, matvmul3x3sf,( float32_t * restrict z, 
                           const float32_t * restrict x, 
                           const float32_t * restrict y, 
                           int L ))
#else

#define sz_f32 ((int)sizeof(float32_t))
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

/* Streaming Order, Floating-Point, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 8
*/
void matvmul3x3sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L )
{
#if 1
  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict px2;
  const xb_vecN_2xf32 * restrict py0;
  const xb_vecN_2xf32 * restrict py1;
  const xb_vecN_2xf32 * restrict py2;
        xb_vecN_2xf32 * restrict pz0;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x10, x11, x12, x20, x21, x22;
  xb_vecN_2xf32 y0, y1, y2;
  xb_vecN_2xf32 z0, z1, z2;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);

  px0 = (const xb_vecN_2xf32 *)(x+L*0);
  px1 = (const xb_vecN_2xf32 *)(x+L*1);
  px2 = (const xb_vecN_2xf32 *)(x+L*2);
  py0 = (const xb_vecN_2xf32 *)(y+L*0);
  py1 = (const xb_vecN_2xf32 *)(y+L*1);
  py2 = (const xb_vecN_2xf32 *)(y+L*2);
  pz0 = (      xb_vecN_2xf32 *)(z);

  for (l = 0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
  {
    /* Load input matrix X */
    BBE_LVN_2XF32_XP(x00, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px0, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x21, px1, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x22, px2, -6*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_2XF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);

    /* Compute matrix Z */
    z0 = BBE_MULN_2XF32(x00, y0);
    z1 = BBE_MULN_2XF32(x10, y0);
    z2 = BBE_MULN_2XF32(x20, y0);
    BBE_MULAN_2XF32(z0, x01, y1);
    BBE_MULAN_2XF32(z1, x11, y1);
    BBE_MULAN_2XF32(z2, x21, y1);
    BBE_MULAN_2XF32(z0, x02, y2);
    BBE_MULAN_2XF32(z1, x12, y2);
    BBE_MULAN_2XF32(z2, x22, y2);

    /* Save results */
    BBE_SVN_2XF32_XP(z0, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z1, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z2, pz0, -2*L*sz_f32+2*BBE_SIMD_WIDTH);
  }
#else
  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict px2;
  const xb_vecN_2xf32 * restrict px0_;
  const xb_vecN_2xf32 * restrict px1_;
  const xb_vecN_2xf32 * restrict px2_;
  const xb_vecN_2xf32 * restrict py0;
  const xb_vecN_2xf32 * restrict py1;
  const xb_vecN_2xf32 * restrict py2;
        xb_vecN_2xf32 * restrict pz0;
  int l;

  xb_vecN_2xf32 x00, x01, x02, x10, x11, x12, x20, x21, x22;
  xb_vecN_2xf32 y0, y1, y2;
  xb_vecN_2xf32 z0, z1, z2;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT((L&(BBE_SIMD_WIDTH/2-1)) == 0);
  if (L <= 0) return;

  px0 = (const xb_vecN_2xf32 *)(x+L*0);
  px1 = (const xb_vecN_2xf32 *)(x+L*1);
  px2 = (const xb_vecN_2xf32 *)(x+L*2);
  px0_= (const xb_vecN_2xf32 *)(x+L*0+(BBE_SIMD_WIDTH/2));
  px1_= (const xb_vecN_2xf32 *)(x+L*1+(BBE_SIMD_WIDTH/2));
  px2_= (const xb_vecN_2xf32 *)(x+L*2+(BBE_SIMD_WIDTH/2));
  py0 = (const xb_vecN_2xf32 *)(y+L*0);
  py1 = (const xb_vecN_2xf32 *)(y+L*1);
  py2 = (const xb_vecN_2xf32 *)(y+L*2);
  pz0 = (      xb_vecN_2xf32 *)(z);

  for (l = 0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++)
  {
    /* Load input matrix X */
    BBE_LVN_2XF32_XP(x00, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px0, -6*L*sz_f32+4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x21, px1, -6*L*sz_f32+4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x22, px2, -6*L*sz_f32+4*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_2XF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);

    /* Compute matrix Z */
    z0 = BBE_MULN_2XF32(x00, y0);
    z1 = BBE_MULN_2XF32(x10, y0);
    z2 = BBE_MULN_2XF32(x20, y0);
    BBE_MULAN_2XF32(z0, x01, y1);
    BBE_MULAN_2XF32(z1, x11, y1);
    BBE_MULAN_2XF32(z2, x21, y1);
    BBE_MULAN_2XF32(z0, x02, y2);
    BBE_MULAN_2XF32(z1, x12, y2);
    BBE_MULAN_2XF32(z2, x22, y2);

    /* Save results */
    BBE_SVN_2XF32_XP(z0, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z1, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z2, pz0, -2*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Load input matrix X */
    BBE_LVN_2XF32_XP(x00, px0_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px1_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px2_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px0_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px1_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px2_,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px0_, -6*L*sz_f32+4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x21, px1_, -6*L*sz_f32+4*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x22, px2_, -6*L*sz_f32+4*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_2XF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);

    /* Compute matrix Z */
    z0 = BBE_MULN_2XF32(x00, y0);
    z1 = BBE_MULN_2XF32(x10, y0);
    z2 = BBE_MULN_2XF32(x20, y0);
    BBE_MULAN_2XF32(z0, x01, y1);
    BBE_MULAN_2XF32(z1, x11, y1);
    BBE_MULAN_2XF32(z2, x21, y1);
    BBE_MULAN_2XF32(z0, x02, y2);
    BBE_MULAN_2XF32(z1, x12, y2);
    BBE_MULAN_2XF32(z2, x22, y2);

    /* Save results */
    BBE_SVN_2XF32_XP(z0, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z1, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z2, pz0, -2*L*sz_f32+2*BBE_SIMD_WIDTH);
  }
  if (L&(BBE_SIMD_WIDTH/2))
  {
    px0 = (const xb_vecN_2xf32 *)(x+L*0+L-(BBE_SIMD_WIDTH/2));
    px1 = (const xb_vecN_2xf32 *)(x+L*1+L-(BBE_SIMD_WIDTH/2));
    px2 = (const xb_vecN_2xf32 *)(x+L*2+L-(BBE_SIMD_WIDTH/2));
    py0 = (const xb_vecN_2xf32 *)(y+L*0+L-(BBE_SIMD_WIDTH/2));
    py1 = (const xb_vecN_2xf32 *)(y+L*1+L-(BBE_SIMD_WIDTH/2));
    py2 = (const xb_vecN_2xf32 *)(y+L*2+L-(BBE_SIMD_WIDTH/2));
    pz0 = (      xb_vecN_2xf32 *)(z    +L-(BBE_SIMD_WIDTH/2));

    /* Load input matrix X */
    BBE_LVN_2XF32_XP(x00, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x01, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x02, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x10, px0,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x11, px1,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x12, px2,  3*L*sz_f32);
    BBE_LVN_2XF32_XP(x20, px0, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x21, px1, -6*L*sz_f32+2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(x22, px2, -6*L*sz_f32+2*BBE_SIMD_WIDTH);

    /* Load input matrix Y */
    BBE_LVN_2XF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);

    /* Compute matrix Z */
    z0 = BBE_MULN_2XF32(x00, y0);
    z1 = BBE_MULN_2XF32(x10, y0);
    z2 = BBE_MULN_2XF32(x20, y0);
    BBE_MULAN_2XF32(z0, x01, y1);
    BBE_MULAN_2XF32(z1, x11, y1);
    BBE_MULAN_2XF32(z2, x21, y1);
    BBE_MULAN_2XF32(z0, x02, y2);
    BBE_MULAN_2XF32(z1, x12, y2);
    BBE_MULAN_2XF32(z2, x22, y2);

    /* Save results */
    BBE_SVN_2XF32_XP(z0, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z1, pz0,    L*sz_f32);
    BBE_SVN_2XF32_XP(z2, pz0, -2*L*sz_f32+2*BBE_SIMD_WIDTH);
  }
#endif
} /* matvmul3x3sf() */
#endif

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
    Real Matrix by Complex Matrix/Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"


/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

/*-------------------------------------------------------------------------
Real Matrix by Complex Matrix/Vector Multiply 

Description: These functions perform pairwise multiplication of left-hand
real matrices by right-hand complex matrices or vectors. Both the block order
and streaming order are allowed for input/output matrix sequences.

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand real matrices
y[L*Sy]     Sequence of right-hand complex matrices or vectors
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

/* Streaming Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul8x8s ( complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q )
{
  int l;

  int L8 = (L >> LOG2_BBE_SIMD_WIDTH);

  int _l = 0;
  int16_t _k = 0;

  const int __k = (L8 << 16) + 1;

  xb_vecNx16 xx00, xx01, xx02, xx03, zero;
  xb_vecNx16 x00l, x00h, x01l, x01h, x02l, x02h, x03l, x03h;
  xb_vecNx16 y0l, y0h, y1l, y1h;

  xb_vecNx16 zl, zh;

  xb_vecNx40 Zl, Zh;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  if (L <= 0) return;
  zero = 0;

  __Pragma("ymemory( py )");
  __Pragma("loop_count min=8");
  for (l = 0; l<8 * L8; l++)
  {
    _k = BBE_ADDMOD16U(_k, __k);
    XT_MOVEQZ(_l, L, _k);

    BBE_LVNX16_XP(xx00, px, 2 * L);
    BBE_LVNX16_XP(xx01, px, 2 * L);
    BBE_LVNX16_XP(xx02, px, 2 * L);
    BBE_LVNX16_XP(xx03, px, 2 * L);

    BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x02h, x02l, zero, xx02, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x03h, x03l, zero, xx03, BBE_DSELI_INTERLEAVE_1);

    y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y0l, py, 2 * 2 * L);
    y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1l, py, 2 * 2 * L);

    Zl = BBE_MULNX16C(x00l, y0l); Zh = BBE_MULNX16C(x00h, y0h);
    BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

    y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y0l, py, 2 * 2 * L);
    y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1l, py, 2 * 2 * L);

    BBE_MULANX16C(Zl, x02l, y0l); BBE_MULANX16C(Zh, x02h, y0h);
    BBE_MULANX16C(Zl, x03l, y1l); BBE_MULANX16C(Zh, x03h, y1h);

    BBE_LVNX16_XP(xx00, px, 2 * L);
    BBE_LVNX16_XP(xx01, px, 2 * L);
    BBE_LVNX16_XP(xx02, px, 2 * L);
    BBE_LVNX16_XP(xx03, px, -7 * 2 * L + 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x02h, x02l, zero, xx02, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(x03h, x03l, zero, xx03, BBE_DSELI_INTERLEAVE_1);

    y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y0l, py, 2 * 2 * L);
    y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1l, py, 2 * 2 * L);

    BBE_MULANX16C(Zl, x00l, y0l); BBE_MULANX16C(Zh, x00h, y0h);
    BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

    y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y0l, py, 2 * 2 * L);
    y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1l, py, -7 * 2 * 2 * L + 2 * 2 * BBE_SIMD_WIDTH);

    BBE_MULANX16C(Zl, x02l, y0l); BBE_MULANX16C(Zh, x02h, y0h);
    BBE_MULANX16C(Zl, x03l, y1l); BBE_MULANX16C(Zh, x03h, y1h);

    /* Save results */

    zl = BBE_PACKVNX40(Zl, q);   zh = BBE_PACKVNX40(Zh, q);

    BBE_SVNX16_IP(zl, pz, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(zh, pz, 2 * BBE_SIMD_WIDTH);

    px = (xb_vecNx16 *)XT_ADDX2(7 * _l, (int32_t)px);
    py = (xb_vecNx16 *)XT_SUB((int32_t)py, 2 * 2 * _l);

    _l = 0;
  }
} /* rcmatvmul8x8s() */

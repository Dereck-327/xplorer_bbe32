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

/* Streaming Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul4x4s ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  int l, i;

  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 x0l, x1l, x2l, x3l;
  xb_vecNx16 x0h, x1h, x2h, x3h;

  xb_vecNx16 y0l, y1l, y2l, y3l;
  xb_vecNx16 y0h, y1h, y2h, y3h;

  xb_vecNx16 zero;

  xb_vecNx16 zl, zh;

  xb_vecNx40 Acc_l, Acc_h;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 *          px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *          py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

  NASSERT(L%BBE_SIMD_WIDTH == 0);

  zero = 0;
  if (L <= 0) return;
  for (i = 0; i<4; i++)
  {
    __Pragma("ymemory( py )");
    __Pragma("loop_count min=16");
    for (l = 0; l<L; l += BBE_SIMD_WIDTH)
    {
      /*Load input matrix X*/
      BBE_LVNX16_XP(x0, px, 2 * L);
      BBE_LVNX16_XP(x1, px, 2 * L);
      BBE_LVNX16_XP(x2, px, 2 * L);
      BBE_LVNX16_XP(x3, px, -3 * 2 * L + 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(x0h, x0l, zero, x0, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x1h, x1l, zero, x1, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x2h, x2l, zero, x2, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x3h, x3l, zero, x3, BBE_DSELI_INTERLEAVE_1);

      /*Load zero columns of input matrix Y*/
      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y0l, py, 4 * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y1l, py, 4 * 2 * 2 * L);
      y2h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y2l, py, 4 * 2 * 2 * L);
      y3h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y3l, py, -11 * 2 * 2 * L);

      /* Multiply */
      Acc_l = BBE_MULRNX16C(x0l, y0l, q); Acc_h = BBE_MULRNX16C(x0h, y0h, q);
      BBE_MULANX16C(Acc_l, x1l, y1l);     BBE_MULANX16C(Acc_h, x1h, y1h);
      BBE_MULANX16C(Acc_l, x2l, y2l);     BBE_MULANX16C(Acc_h, x2h, y2h);
      BBE_MULANX16C(Acc_l, x3l, y3l);     BBE_MULANX16C(Acc_h, x3h, y3h);

      /* Pack and save rezult */
      zl = BBE_PACKVNX40(Acc_l, q); zh = BBE_PACKVNX40(Acc_h, q);
      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(zl, pz, 2 * 2 * L);

      /* Load first columns of input matrix Y */
      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y0l, py, 4 * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y1l, py, 4 * 2 * 2 * L);
      y2h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y2l, py, 4 * 2 * 2 * L);
      y3h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y3l, py, -11 * 2 * 2 * L);

      /* Multiply */
      Acc_l = BBE_MULRNX16C(x0l, y0l, q); Acc_h = BBE_MULRNX16C(x0h, y0h, q);
      BBE_MULANX16C(Acc_l, x1l, y1l);     BBE_MULANX16C(Acc_h, x1h, y1h);
      BBE_MULANX16C(Acc_l, x2l, y2l);     BBE_MULANX16C(Acc_h, x2h, y2h);
      BBE_MULANX16C(Acc_l, x3l, y3l);     BBE_MULANX16C(Acc_h, x3h, y3h);

      /* Pack and save rezult */
      zl = BBE_PACKVNX40(Acc_l, q); zh = BBE_PACKVNX40(Acc_h, q);
      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(zl, pz, 2 * 2 * L);

      /* Load second columns of input matrix Y */
      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y0l, py, 4 * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y1l, py, 4 * 2 * 2 * L);
      y2h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y2l, py, 4 * 2 * 2 * L);
      y3h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y3l, py, -11 * 2 * 2 * L);

      /* Multiply */
      Acc_l = BBE_MULRNX16C(x0l, y0l, q); Acc_h = BBE_MULRNX16C(x0h, y0h, q);
      BBE_MULANX16C(Acc_l, x1l, y1l);     BBE_MULANX16C(Acc_h, x1h, y1h);
      BBE_MULANX16C(Acc_l, x2l, y2l);     BBE_MULANX16C(Acc_h, x2h, y2h);
      BBE_MULANX16C(Acc_l, x3l, y3l);     BBE_MULANX16C(Acc_h, x3h, y3h);

      /* Pack and save rezult */
      zl = BBE_PACKVNX40(Acc_l, q); zh = BBE_PACKVNX40(Acc_h, q);
      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(zl, pz, 2 * 2 * L);

      /* Load third columns of input matrix Y */
      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y0l, py, 4 * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y1l, py, 4 * 2 * 2 * L);
      y2h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y2l, py, 4 * 2 * 2 * L);
      y3h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(y3l, py, -15 * 2 * 2 * L + 2 * 2 * BBE_SIMD_WIDTH);

      /* Multiply */
      Acc_l = BBE_MULRNX16C(x0l, y0l, q); Acc_h = BBE_MULRNX16C(x0h, y0h, q);
      BBE_MULANX16C(Acc_l, x1l, y1l);     BBE_MULANX16C(Acc_h, x1h, y1h);
      BBE_MULANX16C(Acc_l, x2l, y2l);     BBE_MULANX16C(Acc_h, x2h, y2h);
      BBE_MULANX16C(Acc_l, x3l, y3l);     BBE_MULANX16C(Acc_h, x3h, y3h);

      /* Pack and save rezult */
      zl = BBE_PACKVNX40(Acc_l, q); zh = BBE_PACKVNX40(Acc_h, q);
      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH); BBE_SVNX16_XP(zl, pz, -3 * 2 * 2 * L + 2 * 2 * BBE_SIMD_WIDTH);
    }

    px = (const xb_vecNx16 *)XT_ADDX2(3 * L, (int32_t)px);
    py = (const xb_vecNx16 *)XT_SUB((int32_t)py, 2 * 2 * L);
    pz = (xb_vecNx16 *)XT_ADDX4(3 * L, (int32_t)pz);
  }
} /* rcmatmul4x4s() */

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
#include <string.h>

void rcmatmulnxms_0 (complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int M, int L, int Q);

void rcmatmulnxms_1 (complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int M, int L, int Q);
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

/* Streaming Order, MxN*NxM ->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmulnxms ( complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q )
{
  NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT(Q>=0 && Q<=16);
  NASSERT(L%BBE_SIMD_WIDTH==0);
  
  if ((L <= 0) || (M <= 0)) return;
  if (N <= 0)
  {
    memset(z, 0, 2 * M*M*L*sizeof(int16_t));
    return;
  }
  if (N&1)
  {
    rcmatmulnxms_1 (z, 
                    x, 
                    y, 
                    N, M, L, Q);
  }
  else
  {
    rcmatmulnxms_0 (z, 
                    x, 
                    y, 
                    N, M, L, Q);
  }
} /* rcmatmulnxms() */

void rcmatmulnxms_0(complex_fract16 * restrict z,
              const int16_t * restrict x,
              const complex_fract16 * restrict y,
              int N, int M, int L, int Q)
{
  int l, j, k, p;

  xb_vecNx16 xx00, xx01, zero;
  xb_vecNx16 x00l, x00h, x01l, x01h;
  xb_vecNx16 y0l, y0h, y1l, y1h;

  xb_vecNx16 zl, zh;

  xb_vecNx40 Zl, Zh;

  vsaN q = BBE_MOVVSA32(Q);

  const xb_vecNx16 * restrict px;
  const xb_vecNx16 * restrict py;
  xb_vecNx16 * restrict pz;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L%BBE_SIMD_WIDTH == 0);

  zero = BBE_ZERONX16();

  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    pz = (xb_vecNx16 *)(z + BBE_SIMD_WIDTH*l);
    py = (const xb_vecNx16 *)(y + BBE_SIMD_WIDTH*l);
    px = (const xb_vecNx16 *)(x + BBE_SIMD_WIDTH*l);
    __Pragma("loop_count min=1");
    for (k = 0; k<M; k++)
    {
      BBE_LVNX16_XP(xx00, px, 2 * L);
      BBE_LVNX16_XP(xx01, px, 2 * L);

      BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

      Zl = BBE_MULRNX16C(x00l, y0l, q); Zh = BBE_MULRNX16C(x00h, y0h, q);
      BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

      for (j = 0; j<M - 1; j++)
      {
        __Pragma("ymemory( py )");
        for (p = 0; p<N - 2; p += 2)
        {
          BBE_LVNX16_XP(xx00, px, 2 * L);
          BBE_LVNX16_XP(xx01, px, 2 * L);

          BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
          BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

          y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
          y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

          BBE_MULANX16C(Zl, x00l, y0l); BBE_MULANX16C(Zh, x00h, y0h);
          BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
        }

        px = (const xb_vecNx16 *)XT_ADDX2(-L*N, (int32_t)px);
        py = (const xb_vecNx16 *)XT_ADDX4(-L*N*M + L, (int32_t)py);

        zl = BBE_PACKVNX40(Zl, q);    zh = BBE_PACKVNX40(Zh, q);

        BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_XP(zl, pz, 2 * 2 * L);

        BBE_LVNX16_XP(xx00, px, 2 * L);
        BBE_LVNX16_XP(xx01, px, 2 * L);

        BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
        BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

        y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
        y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

        Zl = BBE_MULRNX16C(x00l, y0l, q); Zh = BBE_MULRNX16C(x00h, y0h, q);
        BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
      }

      __Pragma("ymemory( py )");
      for (p = 0; p<N - 2; p += 2)
      {
        BBE_LVNX16_XP(xx00, px, 2 * L);
        BBE_LVNX16_XP(xx01, px, 2 * L);

        BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
        BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

        y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
        y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

        BBE_MULANX16C(Zl, x00l, y0l); BBE_MULANX16C(Zh, x00h, y0h);
        BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
      }

      py = (const xb_vecNx16 *)XT_ADDX4(-L*N*M + L, (int32_t)py);

      zl = BBE_PACKVNX40(Zl, q);    zh = BBE_PACKVNX40(Zh, q);

      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_XP(zl, pz, 2 * 2 * L);

      py = (const xb_vecNx16 *)XT_SUB((int32_t)py, M * 2 * 2 * L);
    }
  }
} /* rcmatmulnxms_0() */

void rcmatmulnxms_1(complex_fract16 * restrict z,
              const int16_t * restrict x,
              const complex_fract16 * restrict y,
              int N, int M, int L, int Q)
{
  int l, j, k, p;

  xb_vecNx16 xx00, xx01, zero;
  xb_vecNx16 x00l, x00h, x01l, x01h;
  xb_vecNx16 y0l, y0h, y1l, y1h;

  xb_vecNx16 zl, zh;

  xb_vecNx40 Zl, Zh;

  vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 * restrict px;
  xb_vecNx16 * restrict py;
  xb_vecNx16 * restrict pz;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L%BBE_SIMD_WIDTH == 0);

  zero = BBE_ZERONX16();

  pz = (xb_vecNx16 *)z;

  if (N == 1)
  {
    for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
    {
      pz = (xb_vecNx16 *)(z + BBE_SIMD_WIDTH*l);

      for (k = 0; k<M; k++)
      {
        py = (xb_vecNx16 *)(y + BBE_SIMD_WIDTH*l);
        px = (xb_vecNx16 *)(x + L*k + BBE_SIMD_WIDTH*l);

        __Pragma("ymemory( py )");
        for (j = 0; j<M; j++)
        {
          Zl = BBE_ZERONX40(); Zh = BBE_ZERONX40();

          BBE_LVNX16_IP(xx01, px, 0);

          BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

          y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_XP(y1l, py, 2 * 2 * L);

          BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

          zl = BBE_PACKVNX40(Zl, q);    zh = BBE_PACKVNX40(Zh, q);

          BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_XP(zl, pz, 2 * 2 * L);
        }
      }
    }

    return;
  }

  for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
  {
    pz = (xb_vecNx16 *)(z + BBE_SIMD_WIDTH*l);

    for (k = 0; k<M; k++)
    {
      py = (xb_vecNx16 *)(y + BBE_SIMD_WIDTH*l);
      px = (xb_vecNx16 *)(x + N*L*k + BBE_SIMD_WIDTH*l);

      BBE_LVNX16_XP(xx00, px, 2 * L);
      BBE_LVNX16_XP(xx01, px, 2 * L);

      BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

      y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

      Zl = BBE_MULRNX16C(x00l, y0l, q); Zh = BBE_MULRNX16C(x00h, y0h, q);
      BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

      for (j = 0; j<M - 1; j++)
      {
        __Pragma("ymemory( py )");
        for (p = 0; p<(N >> 1) - 1; p++)
        {
          BBE_LVNX16_XP(xx00, px, 2 * L);
          BBE_LVNX16_XP(xx01, px, 2 * L);

          BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
          BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

          y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
          y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

          BBE_MULANX16C(Zl, x00l, y0l); BBE_MULANX16C(Zh, x00h, y0h);
          BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
        }

        BBE_LVNX16_XP(xx01, px, -(N - 1) * 2 * L);

        BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

        y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1l, py, -(N - 1)*M * 2 * 2 * L + 2 * 2 * L);

        BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

        zl = BBE_PACKVNX40(Zl, q);    zh = BBE_PACKVNX40(Zh, q);

        BBE_LVNX16_XP(xx00, px, 2 * L);
        BBE_LVNX16_XP(xx01, px, 2 * L);

        BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
        BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

        y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
        y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

        Zl = BBE_MULRNX16C(x00l, y0l, q); Zh = BBE_MULRNX16C(x00h, y0h, q);
        BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

        BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_XP(zl, pz, 2 * 2 * L);
      }

      __Pragma("ymemory( py )");
      for (p = 0; p<(N >> 1) - 1; p++)
      {
        BBE_LVNX16_XP(xx00, px, 2 * L);
        BBE_LVNX16_XP(xx01, px, 2 * L);

        BBE_DSELNX16I(x00h, x00l, zero, xx00, BBE_DSELI_INTERLEAVE_1);
        BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

        y0h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y0l, py, M * 2 * 2 * L);
        y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1l, py, M * 2 * 2 * L);

        BBE_MULANX16C(Zl, x00l, y0l); BBE_MULANX16C(Zh, x00h, y0h);
        BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);
      }

      BBE_LVNX16_XP(xx01, px, -(N - 1) * 2 * L);

      BBE_DSELNX16I(x01h, x01l, zero, xx01, BBE_DSELI_INTERLEAVE_1);

      y1h = BBE_LVNX16_I(py, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y1l, py, -(N - 1)*M * 2 * 2 * L + 2 * 2 * L);

      BBE_MULANX16C(Zl, x01l, y1l); BBE_MULANX16C(Zh, x01h, y1h);

      zl = BBE_PACKVNX40(Zl, q);    zh = BBE_PACKVNX40(Zh, q);

      BBE_SVNX16_I(zh, pz, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_XP(zl, pz, 2 * 2 * L);
    }
  }
} /* rcmatmulnxms_1() */


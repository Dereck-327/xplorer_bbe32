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

/* Streaming Order, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm8x8s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q )
{
  const int N = 8;
  const int M = N;
  const int L8 = L >> 3;
  vsaN q = BBE_MOVVSA32(Q);
  int i, j, k;

  xb_vecNx16 x00, x01;
  xb_vecNx16 x10, x11;
  xb_vecNx16 y00, y10;

  xb_vecNx40 zout_00, zout_01, zout_10, zout_11;

  const xb_vecNx16 * px00;
  const xb_vecNx16 * restrict px01;
  const xb_vecNx16 * restrict px02;
  const xb_vecNx16 * restrict py01;
  const xb_vecNx16 * restrict py02;
  xb_vecNx16 * restrict pz00;
  xb_vecNx16 * restrict pz01;

  px00 = (const xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;

  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);

  if (L <= 0) return;

  xb_vecNx16 Zero = 0;

  px02 = px00;
  pz00 = ((xb_vecNx16*)y);
  __Pragma("ymemory( px02 )");
  __Pragma("loop_count min=4, factor=4");
  for (i = k = 0; k < L8 * 4; k++)
  {
    int xstride, zstride;
    i = BBE_ADDMOD16U(i, (L8 << 16) | 1);
    xstride = (-N + 1) * 4 * N*L + 2 * BBE_SIMD_WIDTH;
    zstride = 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(xstride, (-N + 1) * 4 * N*L + 2 * BBE_SIMD_WIDTH + 4 * L, i);
    XT_MOVEQZ(zstride, 2 * BBE_SIMD_WIDTH + 17 * 4 * L, i);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    zout_00 = BBE_MAGIRNX16C(x11, x00, q);
    zout_10 = BBE_MULRNX16J(x11, x00, q);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, xstride);
    BBE_MAGIANX16C(zout_00, x11, x00);
    BBE_MULANX16J(zout_10, x11, x00);

    x00 = BBE_PACKVNX40(zout_00, q);
    x11 = BBE_SELNX16I(Zero, x00, BBE_SELI_INTERLEAVE_1_ODD);
    x00 = BBE_SELNX16I(Zero, x00, BBE_SELI_INTERLEAVE_1_EVEN);
    x10 = BBE_PACKVNX40(zout_10, q);
    x01 = BBE_CONJSNX16C(x10);

    BBE_SVNX16_X(x01, pz00, 4 * N*L);
    BBE_SVNX16_X(x10, pz00, 4 * L);
    BBE_SVNX16_X(x11, pz00, 4 * L*(N + 1));
    BBE_SVNX16_XP(x00, pz00, zstride);
  }

  px00 = (const xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;
  for (i = 0; i < N; i += 2) 
  {
    px01 = px00 + i*L8;
    pz00 = ((xb_vecNx16*)y) + i*(1 + N)*L8;
    px02 = px01;

    py01 = (xb_vecNx16 *)XT_ADDX2(4 * L, (int32_t)px01);

    for (j = i + 2; j < N; j += 2) 
    {
      pz00 = ((xb_vecNx16*)y) + (i + j*N)*L8;
      pz01 = ((xb_vecNx16*)y) + (j + i*N)*L8;

      px02 = px01;
      py02 = py01;

      __Pragma("ymemory( py02 )");
      __Pragma("ymemory( px02 )");
      __Pragma("loop_count min=1");
      for (k = 0; k < L8; k++)
      {
        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        zout_00 = BBE_MULRNX16J(y00, x00, q);
        zout_01 = BBE_MULRNX16J(y00, x01, q);
        zout_10 = BBE_MULRNX16J(y10, x00, q);
        zout_11 = BBE_MULRNX16J(y10, x01, q);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x01 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L - 4 * N*M*L + 2 * BBE_SIMD_WIDTH);
        y10 = BBE_LVNX16_X(py02, 4 * L);
        BBE_LVNX16_XP(y00, py02, 4 * N*L - 4 * N*M*L + 2 * BBE_SIMD_WIDTH);

        BBE_MULANX16J(zout_00, y00, x00);
        BBE_MULANX16J(zout_01, y00, x01);
        BBE_MULANX16J(zout_10, y10, x00);
        BBE_MULANX16J(zout_11, y10, x01);

        x00 = BBE_PACKVNX40(zout_10, q);
        x01 = BBE_CONJSNX16C(x00);

        x10 = BBE_PACKVNX40(zout_11, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x01, pz00, 4 * N*L);
        BBE_SVNX16_X(x00, pz01, 4 * L);

        BBE_SVNX16_X(x11, pz00, 4 * (N + 1)*L);
        BBE_SVNX16_X(x10, pz01, 4 * (N + 1)*L);

        x00 = BBE_PACKVNX40(zout_00, q);
        x01 = BBE_CONJSNX16C(x00);

        x10 = BBE_PACKVNX40(zout_01, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x11, pz00, 4 * L);
        BBE_SVNX16_X(x10, pz01, 4 * N*L);

        BBE_SVNX16_IP(x01, pz00, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x00, pz01, 2 * BBE_SIMD_WIDTH);
      }
      py01 = (xb_vecNx16 *)XT_ADDX2(4 * L, (int32_t)py01);
    }
  }
} /* cmatherm8x8s() */

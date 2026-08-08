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

/* Streaming Order, NxM * MxN->NxN, Sx=MxN,Sy=NxN
   Restrictions:
   Q=0..15, L is a multiple of 8
*/
void cmathermnxms_1(complex_fract16 * restrict y,
                    const complex_fract16 * restrict x,
                    int M, int N, int L, int Q);

/* Streaming Order, NxM * MxN->NxN, Sx=MxN,Sy=NxN
   Restrictions:
   Q=0..15, L is a multiple of 8
*/
void cmathermnxms_2(complex_fract16 * restrict y,
                    const complex_fract16 * restrict x, 
                    int M, int N, int L, int Q);

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

/* Streaming Order, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     L must be a multiple of 8
*/
void cmathermnxms(complex_fract16 * restrict y,
              const complex_fract16 * restrict x, 
              int N, int M, int L, int Q )
{
  /* check restrictions */
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);
  if (L<=0 || N<=0 || M<=0) return;
  if (N & 1)
  {
    cmathermnxms_1(y, x, M, N, L, Q);
  }
  else
  {
    cmathermnxms_2(y, x, M, N, L, Q);
  }
} /* cmathermnxms() */

/* Streaming Order, NxM * MxN->NxN, Sx=MxN,Sy=NxN
Restrictions:
L is a multiple of 16
*/
void cmathermnxms_1(complex_fract16 * restrict y,
  const complex_fract16 * restrict x,
  int M, int N, int L, int Q)
{
  const int L8 = L >> 3;
  vsaN q = BBE_MOVVSA32(Q);
  int i, j, k, r;

  xb_vecNx16 x00;
  xb_vecNx16 x01;
  xb_vecNx16 x10;
  xb_vecNx16 x11;
  xb_vecNx16 y00;
  xb_vecNx16 y01;

  xb_vecNx40 z_out00, z_out01;
  xb_vecNx40 z_out10, z_out11;

  xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict px01;
  xb_vecNx16 * restrict px02;
  xb_vecNx16 * restrict py01;
  xb_vecNx16 * restrict py02;
  xb_vecNx16 * restrict pz00;
  xb_vecNx16 * restrict pz01;

  px00 = (xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;

  /* check restrictions */
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);

  for (i = 0; i < N - 1; i += 2) 
  {
    px01 = px00 + i*L8;

    pz00 = ((xb_vecNx16*)y) + i*(1 + N)*L8;
    __Pragma("loop_count min=1");
    for (k = 0; k < L8; k++)
    {
      px02 = px01 + k;

      z_out00 = 0;

      __Pragma("ymemory( px02 )");
      __Pragma("loop_count min=1");
      for (r = 0; r<M; r++)
      {
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        BBE_MULANX16J(z_out00, x00, x00);
      }

      x00 = BBE_PACKVNX40(z_out00, q);

      BBE_SVNX16_IP(x00, pz00, 2 * BBE_SIMD_WIDTH);
    }

    for (j = i + 1; j < N; j += 2)
    {
      py01 = px00 + j*L8;

      pz00 = ((xb_vecNx16*)y) + (i + j*N)*L8;
      pz01 = ((xb_vecNx16*)y) + (j + i*N)*L8;

      for (k = 0; k < L8; k++)
      {
        px02 = px01 + k;
        py02 = py01 + k;

        z_out00 = 0;
        z_out01 = 0;
        z_out10 = 0;
        z_out11 = 0;

        __Pragma("ymemory( px02 )");
        __Pragma("ymemory( py02 )");
        __Pragma("loop_count min=1");
        for (r = 0; r<M; r++)
        {
          x01 = BBE_LVNX16_X(px02, 4 * L);
          y01 = BBE_LVNX16_X(py02, 4 * L);
          BBE_LVNX16_XP(x00, px02, 4 * N*L);
          BBE_LVNX16_XP(y00, py02, 4 * N*L);

          BBE_MULANX16J(z_out00, y00, x00);
          BBE_MULANX16J(z_out01, y00, x01);
          BBE_MULANX16J(z_out10, y01, x00);
          BBE_MULANX16J(z_out11, y01, x01);
        }

        x10 = BBE_PACKVNX40(z_out11, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x11, pz00, 4 * (N + 1)*L);
        BBE_SVNX16_X(x10, pz01, 4 * (N + 1)*L);

        x00 = BBE_PACKVNX40(z_out10, q);
        x01 = BBE_CONJSNX16C(x00);

        BBE_SVNX16_X(x01, pz00, 4 * N*L);
        BBE_SVNX16_X(x00, pz01, 4 * L);

        x10 = BBE_PACKVNX40(z_out01, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x11, pz00, 4 * L);
        BBE_SVNX16_X(x10, pz01, 4 * N*L);

        x00 = BBE_PACKVNX40(z_out00, q);
        x01 = BBE_CONJSNX16C(x00);

        BBE_SVNX16_IP(x01, pz00, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x00, pz01, 2 * BBE_SIMD_WIDTH);
      }
    }
  }

  px01 = px00 + (N - 1)*L8;

  pz00 = ((xb_vecNx16*)y) + (N - 1 + (N - 1)*N)*L8;

  for (k = 0; k < L8; k++)
  {
    px02 = px01 + k;

    z_out00 = 0;

    __Pragma("ymemory( px02 )");
    __Pragma("loop_count min=1");
    for (r = 0; r<M; r++)
    {
      BBE_LVNX16_XP(x00, px02, 4 * N*L);
      BBE_MULANX16J(z_out00, x00, x00);
    }

    x00 = BBE_PACKVNX40(z_out00, q);
    x01 = BBE_PACKVNX40(z_out10, q);

    BBE_SVNX16_IP(x00, pz00, 2 * BBE_SIMD_WIDTH);
  }
}

/* Streaming Order, NxM * MxN->NxN, Sx=MxN,Sy=NxN
Restrictions:
L is a multiple of 16
*/
void cmathermnxms_2(complex_fract16 * restrict y,
              const complex_fract16 * restrict x,
              int M, int N, int L, int Q)
{
  const int L8 = L >> 3;
  vsaN q = BBE_MOVVSA32(Q);
  int i, j, k, r;

  xb_vecNx16 x00;
  xb_vecNx16 x01;
  xb_vecNx16 x10;
  xb_vecNx16 x11;
  xb_vecNx16 y00;
  xb_vecNx16 y01;

  xb_vecNx16 Zero = 0;

  xb_vecNx40 z_out00, z_out01;
  xb_vecNx40 z_out10, z_out11;

  xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict px01;
  xb_vecNx16 * restrict px02;
  xb_vecNx16 * restrict py01;
  xb_vecNx16 * restrict py02;
  xb_vecNx16 * restrict pz00;
  xb_vecNx16 * restrict pz01;

  px00 = (xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;

  /* check restrictions */
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);

  for (i = 0; i < N; i += 2)
  {
    px01 = (xb_vecNx16*)XT_ADD(i * 4 * L, (int32_t)px00);

    pz00 = ((xb_vecNx16*)y) + i*(1 + N)*L8;

    __Pragma("loop_count min=1");
    for (k = 0; k < L8; k++)
    {
      px02 = px01 + k;

      z_out10 = 0;
      z_out00 = 0;

      __Pragma("ymemory( px02 )");
      __Pragma("loop_count min=1");
      for (r = 0; r<M; r++)
      {
        y00 = BBE_LVNX16_X(px02, 4 * L);
        BBE_LVNX16_XP(x00, px02, 4 * N*L);
        BBE_MULANX16J(z_out10, y00, x00);
        BBE_MAGIANX16C(z_out00, y00, x00);
      }

      x10 = BBE_PACKVNX40(z_out10, q);
      x01 = BBE_CONJSNX16C(x10);
      y00 = BBE_PACKVNX40(z_out00, q);
      x00 = BBE_SELNX16I(Zero, y00, BBE_SELI_INTERLEAVE_1_EVEN);
      x11 = BBE_SELNX16I(Zero, y00, BBE_SELI_INTERLEAVE_1_ODD);

      BBE_SVNX16_X(x11, pz00, 4 * (N + 1)*L);
      BBE_SVNX16_X(x01, pz00, 4 * N*L);
      BBE_SVNX16_X(x10, pz00, 4 * L);
      BBE_SVNX16_IP(x00, pz00, 2 * BBE_SIMD_WIDTH);
    }

    for (j = i + 2; j < N; j += 2) 
    {
      py01 = (xb_vecNx16*)XT_ADD(j * 4 * L, (int32_t)px00);

      pz00 = ((xb_vecNx16*)y) + (i + j*N)*L8;
      pz01 = ((xb_vecNx16*)y) + (j + i*N)*L8;

      __Pragma("loop_count min=1");
      for (k = 0; k < L8; k++)
      {
        px02 = px01 + k;
        py02 = py01 + k;

        z_out00 = 0;
        z_out01 = 0;
        z_out10 = 0;
        z_out11 = 0;

        __Pragma("ymemory( px02 )");
        __Pragma("ymemory( py02 )");
        __Pragma("loop_count min=1");
        for (r = 0; r<M; r++)
        {
          x01 = BBE_LVNX16_X(px02, 4 * L);
          y01 = BBE_LVNX16_X(py02, 4 * L);
          BBE_LVNX16_XP(x00, px02, 4 * N*L);
          BBE_LVNX16_XP(y00, py02, 4 * N*L);

          BBE_MULANX16J(z_out00, y00, x00);
          BBE_MULANX16J(z_out01, y00, x01);
          BBE_MULANX16J(z_out10, y01, x00);
          BBE_MULANX16J(z_out11, y01, x01);
        }

        x00 = BBE_PACKVNX40(z_out10, q);
        x01 = BBE_CONJSNX16C(x00);

        x10 = BBE_PACKVNX40(z_out11, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x11, pz00, 4 * (N + 1)*L);
        BBE_SVNX16_X(x10, pz01, 4 * (N + 1)*L);

        BBE_SVNX16_X(x01, pz00, 4 * N*L);
        BBE_SVNX16_X(x00, pz01, 4 * L);

        x00 = BBE_PACKVNX40(z_out00, q);
        x01 = BBE_CONJSNX16C(x00);

        x10 = BBE_PACKVNX40(z_out01, q);
        x11 = BBE_CONJSNX16C(x10);

        BBE_SVNX16_X(x11, pz00, 4 * L);
        BBE_SVNX16_X(x10, pz01, 4 * N*L);

        BBE_SVNX16_IP(x01, pz00, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x00, pz01, 2 * BBE_SIMD_WIDTH);
      }
    }
  }
} /* cmathermnxms_2() */

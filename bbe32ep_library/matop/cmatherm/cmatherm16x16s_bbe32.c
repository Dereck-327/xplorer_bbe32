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

#define ITER(z1,z0,px0,py,N,L)     \
{                                  \
    xb_vecNx16 x0,y0,x1;           \
    x1=BBE_LVNX16_X(  px0,   4*L); \
    BBE_LVNX16_XP(x0, px0, 4*N*L); \
    BBE_LVNX16_XP(y0, py , 4*N*L); \
    BBE_MULANX16J(z0, y0 ,x0);     \
    BBE_MULANX16J(z1, y0 ,x1);     \
}

#define ITER1(z1,z0,px0,px1,py,N,L)\
{                                  \
    xb_vecNx16 x0,y0,x1;           \
    BBE_LVNX16_XP(x1, px1, 4*N*L); \
    BBE_LVNX16_XP(x0, px0, 4*N*L); \
    BBE_LVNX16_XP(y0, py , 4*N*L); \
    BBE_MULANX16J(z0, y0 ,x0);     \
    BBE_MULANX16J(z1, y0 ,x1);     \
}

#define ITER0(z1,z0,px,N,L)             \
{                                       \
    xb_vecNx16 x00,x11;                 \
    x11 = BBE_LVNX16_X(px, 4*L);        \
    BBE_LVNX16_XP(x00, px, 4*N*L);      \
    BBE_MAGIANX16C(z0, x11 ,x00);       \
    BBE_MULANX16J (z1, x11 ,x00);       \
}
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

/* Streaming Order, 16x16*16x16->16x16, Sx=256, Sy=256
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm16x16s ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L, int Q )
{
  int temp[3];
  const int N = 16;
  const int L8 = L >> 3;
  vsaN q = BBE_MOVVSA32(Q);
  int i, j, k;

  xb_vecNx16 x00, x01;
  xb_vecNx16 x10, x11;
  xb_vecNx16 y00;
  xb_vecNx16 Zero = 0;
  xb_vecNx40 zout_00, zout_01, zout_10;

  const xb_vecNx16 * px00;
  const xb_vecNx16 * restrict px01;
  const xb_vecNx16 * restrict px02;
  const xb_vecNx16 * restrict py02;
  xb_vecNx16 * restrict pz00;
  xb_vecNx16 * restrict pz01;

  /* cached constants */
  temp[0] = 2 * BBE_SIMD_WIDTH;
  temp[1] = 2 * BBE_SIMD_WIDTH + 4 * (N - 1)*L;

  px00 = (const xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;

  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);

  if (L <= 0) return;

  px02 = px00;
  pz00 = ((xb_vecNx16*)y);
  __Pragma("loop_count min=8, factor=8");
  for (i = k = 0; k < L8 * 8; k++)
  {
    int xstride, zstride;
    i = BBE_ADDMOD16U(i, (L8 << 16) | 1);
    xstride = (-N + 1) * 4 * N*L + 2 * BBE_SIMD_WIDTH;
    zstride = 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(xstride, (-N + 1) * 4 * N*L + 2 * BBE_SIMD_WIDTH + 4 * L, i);
    XT_MOVEQZ(zstride, 2 * BBE_SIMD_WIDTH + (N * 2 + 1) * 4 * L, i);

    x11 = BBE_LVNX16_X(px02, 4 * L);
    BBE_LVNX16_XP(x00, px02, 4 * N*L);
    zout_00 = BBE_MAGIRNX16C(x11, x00, q);
    zout_10 = BBE_MULRNX16J(x11, x00, q);

    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)
    ITER0(zout_10, zout_00, px02, N, L)

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

  
  temp[0] = 2 * BBE_SIMD_WIDTH;
  temp[1] = 2 * BBE_SIMD_WIDTH + 4 * (N - 1)*L;
  temp[2] = -15 * 4 * N*L + 2 * BBE_SIMD_WIDTH;

  for (i = 0; i < N; i += 2) 
  {
    px02 = px00;
    px01 = py02 = (const xb_vecNx16*)XT_ADD(4 * L, (uintptr_t)px02);
    pz00 = (xb_vecNx16*)XT_ADD(4 * N*L, (uintptr_t)y);
    pz01 = (xb_vecNx16*)XT_ADD(4 * L, (uintptr_t)y);
    __Pragma("ymemory( px02 )");
    __Pragma("loop_count min=2");
    for (j = k = 0; k < L8*(N - i - 1); k++)
    {
      int xstride, xstride0, zstride, zstride1;
      xstride0 = XT_L32I(temp, 8);
      zstride = XT_L32I(temp, 0);
      zstride1 = XT_L32I(temp, 4);
      xstride = XT_SUB(xstride0, 4 * L);
      //xstride=-15*4*N*L + 2*BBE_SIMD_WIDTH-4*L;
      j = BBE_ADDMOD16U(j, (L8 << 16) | 1);
      XT_MOVNEZ(xstride, xstride0, j);
      XT_MOVEQZ(zstride, zstride1, j);

      BBE_LVNX16_XP(x01, px01, 4 * N*L);
      BBE_LVNX16_XP(x00, px02, 4 * N*L);
      BBE_LVNX16_XP(y00, py02, 4 * N*L);
      zout_01 = BBE_MULNX16J(y00, x01);
      zout_00 = BBE_MULNX16J(y00, x00);

      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)
      ITER1(zout_01, zout_00, px02, px01, py02, N, L)

      BBE_LVNX16_XP(x01, px01, xstride);
      BBE_LVNX16_XP(x00, px02, xstride);
      BBE_LVNX16_XP(y00, py02, xstride0);
      BBE_MULANX16J(zout_00, y00, x00);
      BBE_MULANX16J(zout_01, y00, x01);

      x00 = BBE_PACKVNX40(zout_00, q);
      x01 = BBE_CONJSNX16C(x00);
      x10 = BBE_PACKVNX40(zout_01, q);
      x11 = BBE_CONJSNX16C(x10);

      BBE_SVNX16_X(x11, pz00, 4 * L);
      BBE_SVNX16_X(x10, pz01, 4 * N*L);
      BBE_SVNX16_XP(x01, pz00, zstride);
      BBE_SVNX16_IP(x00, pz01, 2 * BBE_SIMD_WIDTH);
    }  
    px00 = (const xb_vecNx16*)XT_ADDX2(4 * L, (uintptr_t)px00);
    y = (complex_fract16*)XT_ADDX2(4 * N*L, (uintptr_t)y);
    y = (complex_fract16*)XT_ADDX2(4 * L, (uintptr_t)y);
  }           
} /* cmatherm16x16s() */

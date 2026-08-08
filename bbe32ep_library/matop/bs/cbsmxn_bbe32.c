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
    Packed to Streaming Conversion for Real and Complex Matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"
#include <string.h>
static void _cbs2x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
    (void)MN;
    cbs2x1(y, x, L);
}
static void _cbs4x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
    (void)MN;
    cbs4x1(y, x, L);
}
static void _cbs8x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
    (void)MN;
    cbs8x1(y, x, L);
}

/*-------------------------------------------------------------------------
Packed to Streaming Conversion for Real and Complex Matrices

Description: convert a sequence of L MxN matrices from packed (block) order
to streaming order. Use rbs*() functions for real data, and cbs*() functions -
for complex data.

Representation:
<r|c>bs<size>    16-bit fixed-point data
<r|c>bs<size>f   IEEE-754 Std single precision floating-point data

Storage size SX denotes the number of data elements required to store an
MxN matrix X in block order. If matrix size M*N is less than the SIMD vector
size for appropriate data type, then the storage size equals M*N rounded up
to the next power of two, otherwise SX equals M*N rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for real fixed-point data 2*BBE_SIMD_WIDTH/sizeof(int16_t) == 16
  - for complex fixed-point data 2*BBE_SIMD_WIDTH/sizeof(complex_fract16) == 8
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4

Parameters:
Input:
x[L][SX]   Input sequence of MxN matrices, block order
M          Number of rows in a matrix
N          Number of columns in a matrix
L          Number of matrices
Output:
y[M*N][L]  Output sequence of matrices, streaming order

Restrictions:
x,y      Aligned on 32-byte boundary
x,y      Must not overlap
L        Must be a multiple of:
           16 for real fixed-point data 
            8 for complex fixed-point data and real floating-point data
            4 for complex floating-point data
-------------------------------------------------------------------------*/
void cbsmxn_odd(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
  const xb_vecNx16* restrict X = (const xb_vecNx16*)x;
  xb_vecNx16* restrict Y = (xb_vecNx16*)y;
  xb_vecNx16* y1; const xb_vecNx16* x1;
  int xstride, ystride;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  xb_vecNx40 W0;
  vboolN b;
  int k, l, Sx, off1, off2, off3, off4, off5, off6;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
  NASSERT((MN & 7) != 0);
  if (L <= 0) return;
  Sx = ((MN + 7) >> 3) << 5;

  x1 = X;
  y1 = Y;
  __Pragma("loop_count min=1");
  for (k = l = 0; l<(MN >> 3)*(L >> 3); l++)
  {
    k = BBE_ADDMOD16U(k, ((MN >> 3) << 16) | 1);
    xstride = 2 * BBE_SIMD_WIDTH - 7 * Sx;
    ystride = (L << 2);
    XT_MOVEQZ(xstride, 2 * BBE_SIMD_WIDTH*(1 + (Sx >> 2) - (MN >> 3)) - 7 * Sx, k);
    XT_MOVEQZ(ystride, (L << 2)*(1 - 8 * (MN >> 3)) + 2 * BBE_SIMD_WIDTH, k);
    BBE_LVNX16_XP(X0, X, Sx);
    BBE_LVNX16_XP(X1, X, Sx);
    BBE_LVNX16_XP(X2, X, Sx);
    BBE_LVNX16_XP(X3, X, Sx);
    DEINTLV1(Y1, Y0, X1, X0);
    DEINTLV1(Y3, Y2, X3, X2);
    DEINTLV2(X1, X0, Y2, Y0);
    DEINTLV2(X3, X2, Y3, Y1);

    BBE_LVNX16_XP(X4, X, Sx);
    BBE_LVNX16_XP(X5, X, Sx);
    BBE_LVNX16_XP(X6, X, Sx);
    BBE_LVNX16_XP(X7, X, xstride);
    DEINTLV1(Y5, Y4, X5, X4);
    DEINTLV1(Y7, Y6, X7, X6);
    DEINTLV2(X5, X4, Y6, Y4);
    DEINTLV2(X7, X6, Y7, Y5);

    DEINTLV3(Y1, Y0, X4, X0);
    DEINTLV3(Y3, Y2, X5, X1);
    DEINTLV3(Y5, Y4, X6, X2);
    DEINTLV3(Y7, Y6, X7, X3);
    BBE_SVNX16_XP(Y0, Y, (L << 2));
    BBE_SVNX16_XP(Y1, Y, (L << 2));
    BBE_SVNX16_XP(Y2, Y, (L << 2));
    BBE_SVNX16_XP(Y3, Y, (L << 2));
    BBE_SVNX16_XP(Y4, Y, (L << 2));
    BBE_SVNX16_XP(Y5, Y, (L << 2));
    BBE_SVNX16_XP(Y6, Y, (L << 2));
    BBE_SVNX16_XP(Y7, Y, ystride);
  }

  /* last MN&7 samples */

  X0 = BBE_SEQNX16();
  X1 = BBE_MOVVA16((MN & 7));
  b = BBE_LTNX16(X0, X1);
  X2 = BBE_MOVVA16((int16_t)(L << 2));
  X3 = 0;
  X0 = BBE_MOVNX16T(X0, X3, b);
  W0 = BBE_MULNX16(X0, X2);
  X0 = BBE_MOVVWL(W0);
  off1 = BBE_EXTRNX16C(X0, 1);
  off2 = BBE_EXTRNX16C(X0, 2);
  off3 = BBE_EXTRNX16C(X0, 3);
  off4 = BBE_EXTRNX16C(X0, 4);
  off5 = BBE_EXTRNX16C(X0, 5);
  off6 = BBE_EXTRNX16C(X0, 6);

  X = (const xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH*(MN >> 3), (uintptr_t)x1);
  Y = (xb_vecNx16*)XT_ADD(8 * (L << 2)*(MN >> 3), (uintptr_t)y1);
  __Pragma("loop_count min=1");
  for (k = 0; k<(L >> 3); k++)
  {
    BBE_LVNX16_XP(X0, X, Sx);
    BBE_LVNX16_XP(X1, X, Sx);
    BBE_LVNX16_XP(X2, X, Sx);
    BBE_LVNX16_XP(X3, X, Sx);
    DEINTLV1(Y1, Y0, X1, X0);
    DEINTLV1(Y3, Y2, X3, X2);
    DEINTLV2(X1, X0, Y2, Y0);
    DEINTLV2(X3, X2, Y3, Y1);

    BBE_LVNX16_XP(X4, X, Sx);
    BBE_LVNX16_XP(X5, X, Sx);
    BBE_LVNX16_XP(X6, X, Sx);
    BBE_LVNX16_XP(X7, X, Sx);
    DEINTLV1(Y5, Y4, X5, X4);
    DEINTLV1(Y7, Y6, X7, X6);
    DEINTLV2(X5, X4, Y6, Y4);
    DEINTLV2(X7, X6, Y7, Y5);

    DEINTLV3(Y1, Y0, X4, X0);
    DEINTLV3(Y3, Y2, X5, X1);
    DEINTLV3(Y5, Y4, X6, X2);
    DEINTLV3(Y7, Y6, X7, X3);
    BBE_SVNX16_X(Y1, Y, off1);
    BBE_SVNX16_X(Y2, Y, off2);
    BBE_SVNX16_X(Y3, Y, off3);
    BBE_SVNX16_X(Y4, Y, off4);
    BBE_SVNX16_X(Y5, Y, off5);
    BBE_SVNX16_X(Y6, Y, off6);
    BBE_SVNX16_IP(Y0, Y, 2 * BBE_SIMD_WIDTH);
  }
} /* cbsmxn_odd() */

/* MN is a multiple of 8 */
void cbsmxn_even(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
  csbmxn_even(y, x, L, MN);
} /* cbsmxn_even() */

/* Sx=<see the description>, Sy=M*N */
void cbsmxn(complex_fract16 * restrict y, const complex_fract16 * restrict x, int M, int N, int L)
{
  int MN = M*N;
  /* virtual functions supporting sizes 2...8 */
  typedef void(*fnbs)(complex_fract16*, const complex_fract16*, int, int);
  static const fnbs fnsmall[] = { _cbs2x1, cbs3x1, _cbs4x1, cbs5_7x1, cbs5_7x1, cbs5_7x1, _cbs8x1 };
  fnbs fn;
  if (L <= 0 || M <= 0 || N <= 0) return;
  if (MN == 1)  { memcpy(y, x, L * 2 * 2); return; }
  fn = (MN <= 8) ? fnsmall[MN - 2] : ((MN & 7) ? cbsmxn_odd : cbsmxn_even);
  fn(y, x, MN, L);
    
} /* cbsmxn() */

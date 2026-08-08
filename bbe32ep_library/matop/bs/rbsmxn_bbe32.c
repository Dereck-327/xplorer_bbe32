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
/* redirectors */
static void _rbs2x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L)
{
    (void)MN;
    rbs2x1(y,x,L);
}
static void _rbs4x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L)
{
    (void)MN;
    rbs4x1(y,x,L);
}
static void _rbs8x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L)
{
    (void)MN;
    rbs8x1(y,x,L);
}
static void _rbs16x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L)
{
    (void)MN;
    rbs4x4(y,x,L);
}

#define BBE_MOVVW(hi,lo,w) { hi=BBE_MOVVWH(w);lo=BBE_MOVVWL(w);}

/*   Sx=M*N, Sy=M*N
Restrictions:
MN is a multiple of 16. If not it computes MN&(~15) columns
*/
void rbsmxn_even(int16_t* restrict y, const int16_t* restrict x, int MN, int L)
{
  const xb_vecNx16* restrict w = (const xb_vecNx16*)x;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  int k, l, _8L, _7L, back;
  int wstride, zstride;
  int wstridebig, zstridebig;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  NASSERT(MN >= 2);
  NASSERT((2 * MN) % BBE_SIMD_WIDTH == 0);
  _8L = (L << 4);
  _7L = (L << 1) - _8L;
  back = 2 * BBE_SIMD_WIDTH - _8L + _7L;
  wstridebig = 2 * BBE_SIMD_WIDTH - (L - 1)*(MN << 1);
  zstridebig = back + 2 * BBE_SIMD_WIDTH*(L - (L >> 4));
  __Pragma("loop_count min=1");
  for (k = l = 0; l<(L >> 4)*(MN / (2 * BBE_SIMD_WIDTH / 2)); l++)
  {
    k = BBE_ADDMOD16U(k, ((L >> 4) << 16) | 1);
    wstride = (MN << 1);
    zstride = back;
    XT_MOVEQZ(wstride, wstridebig, k);
    XT_MOVEQZ(zstride, zstridebig, k);

    BBE_LVNX16_XP(Y0, w, ((MN) << 1));
    BBE_LVNX16_XP(Y1, w, ((MN) << 1));
    BBE_LVNX16_XP(Y2, w, ((MN) << 1));
    BBE_LVNX16_XP(Y3, w, ((MN) << 1));
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y4, w, ((MN) << 1));
    BBE_LVNX16_XP(Y5, w, ((MN) << 1));
    BBE_LVNX16_XP(Y6, w, ((MN) << 1));
    BBE_LVNX16_XP(Y7, w, ((MN) << 1));
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y8, w, ((MN) << 1));
    BBE_LVNX16_XP(Y9, w, ((MN) << 1));
    BBE_LVNX16_XP(YA, w, ((MN) << 1));
    BBE_LVNX16_XP(YB, w, ((MN) << 1));
    BBE_DSELNX16I(X9, X8, Y9, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XB, XA, YB, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YA, Y8, XA, X8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YB, Y9, XB, X9, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(YC, w, ((MN) << 1));
    BBE_LVNX16_XP(YD, w, ((MN) << 1));
    BBE_LVNX16_XP(YE, w, ((MN) << 1));
    BBE_LVNX16_XP(YF, w, wstride);
    BBE_DSELNX16I(XD, XC, YD, YC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XE, YF, YE, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YE, YC, XE, XC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YF, YD, XF, XD, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(X4, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X1, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X6, X2, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X3, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(XC, X8, YC, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XD, X9, YD, Y9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XE, XA, YE, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XB, YF, YB, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(Y8, Y0, X8, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y0, z, _8L);
    BBE_SVNX16_XP(Y8, z, _7L);
    BBE_DSELNX16I(Y9, Y1, X9, X1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y1, z, _8L);
    BBE_SVNX16_XP(Y9, z, _7L);
    BBE_DSELNX16I(YA, Y2, XA, X2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y2, z, _8L);
    BBE_SVNX16_XP(YA, z, _7L);
    BBE_DSELNX16I(YB, Y3, XB, X3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y3, z, _8L);
    BBE_SVNX16_XP(YB, z, _7L);
    BBE_DSELNX16I(YC, Y4, XC, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y4, z, _8L);
    BBE_SVNX16_XP(YC, z, _7L);
    BBE_DSELNX16I(YD, Y5, XD, X5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y5, z, _8L);
    BBE_SVNX16_XP(YD, z, _7L);
    BBE_DSELNX16I(YE, Y6, XE, X6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y6, z, _8L);
    BBE_SVNX16_XP(YE, z, _7L);
    BBE_DSELNX16I(YF, Y7, XF, X7, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y7, z, _8L);
    BBE_SVNX16_XP(YF, z, zstride);
  }
} /* rbsmxn_even() */


/*   Sx=M*N, Sy=M*N
Restrictions:
N*M is not a multiple of 16 and >16
*/
void rbsmxn_odd(int16_t* restrict y, const int16_t* restrict x, int MN, int L)
{
  const xb_vecNx16* restrict w = (const xb_vecNx16*)x;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  int k, l, _8L, _MN;
  int wstride, wstridebig, zstride, zstridebig, back;
  unsigned int idx1_0, idx1_1, idx1_2, idx1_3, idx1_4, idx1_5, idx1_6;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
  xb_vecNx40 W0;
  vboolN b;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
  _MN = MN & 15;
  MN = (MN + 15)&(~15);
  _8L = (L << 4);
  back = 2 * BBE_SIMD_WIDTH + (L << 1) - (L << 4);
  wstridebig = 2 * BBE_SIMD_WIDTH - (L - 1)*(MN << 1);
  zstridebig = back + 2 * BBE_SIMD_WIDTH*(L - (L >> 4));

  for (k = l = 0; l<(L >> 4)*((MN >> 4) - 1); l++)
  {
    k = BBE_ADDMOD16U(k, ((L >> 4) << 16) | 1);
    wstride = (MN << 1);
    zstride = back;
    XT_MOVEQZ(wstride, wstridebig, k);
    XT_MOVEQZ(zstride, zstridebig, k);

    BBE_LVNX16_XP(Y0, w, (MN << 1));
    BBE_LVNX16_XP(Y1, w, (MN << 1));
    BBE_LVNX16_XP(Y2, w, (MN << 1));
    BBE_LVNX16_XP(Y3, w, (MN << 1));
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y4, w, (MN << 1));
    BBE_LVNX16_XP(Y5, w, (MN << 1));
    BBE_LVNX16_XP(Y6, w, (MN << 1));
    BBE_LVNX16_XP(Y7, w, (MN << 1));
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y8, w, (MN << 1));
    BBE_LVNX16_XP(Y9, w, (MN << 1));
    BBE_LVNX16_XP(YA, w, (MN << 1));
    BBE_LVNX16_XP(YB, w, (MN << 1));
    BBE_DSELNX16I(X9, X8, Y9, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XB, XA, YB, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YA, Y8, XA, X8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YB, Y9, XB, X9, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(YC, w, (MN << 1));
    BBE_LVNX16_XP(YD, w, (MN << 1));
    BBE_LVNX16_XP(YE, w, (MN << 1));
    BBE_LVNX16_XP(YF, w, wstride);
    BBE_DSELNX16I(XD, XC, YD, YC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XE, YF, YE, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YE, YC, XE, XC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YF, YD, XF, XD, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(X4, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X1, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X6, X2, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X3, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(XC, X8, YC, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XD, X9, YD, Y9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XE, XA, YE, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XB, YF, YB, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(Y8, Y0, X8, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(Y8, z, _8L);
    BBE_SVNX16_XP(Y0, z, (L << 1));
    BBE_DSELNX16I(Y9, Y1, X9, X1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(Y9, z, _8L);
    BBE_SVNX16_XP(Y1, z, (L << 1));
    BBE_DSELNX16I(YA, Y2, XA, X2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YA, z, _8L);
    BBE_SVNX16_XP(Y2, z, (L << 1));
    BBE_DSELNX16I(YB, Y3, XB, X3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YB, z, _8L);
    BBE_SVNX16_XP(Y3, z, (L << 1));
    BBE_DSELNX16I(YC, Y4, XC, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YC, z, _8L);
    BBE_SVNX16_XP(Y4, z, (L << 1));
    BBE_DSELNX16I(YD, Y5, XD, X5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YD, z, _8L);
    BBE_SVNX16_XP(Y5, z, (L << 1));
    BBE_DSELNX16I(YE, Y6, XE, X6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YE, z, _8L);
    BBE_SVNX16_XP(Y6, z, (L << 1));
    BBE_DSELNX16I(YF, Y7, XF, X7, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YF, z, _8L);
    BBE_SVNX16_XP(Y7, z, zstride);
  }
  /* last iteration with incomplete saves */

  /* compute offsets, settings 0 for position not to be saved */
  /*
  for (k=0; k<=15; k++)
  {
  off[k]=(_MN-1)>=k ? (L<<2)*k: 0;
  }
  */
  X0 = BBE_SEQNX16();
  X1 = BBE_MOVVA16(_MN);
  b = BBE_LTNX16(X0, X1);
  X2 = BBE_MOVVA16((int16_t)(L << 1));
  X3 = 0;
  X0 = BBE_MOVNX16T(X0, X3, b);
  W0 = BBE_MULNX16(X0, X2);
  X0 = BBE_MOVVWH(W0);
  idx1_0 = BBE_EXTRNX16C(X0, 0);
  idx1_1 = BBE_EXTRNX16C(X0, 1);
  idx1_2 = BBE_EXTRNX16C(X0, 2);
  idx1_3 = BBE_EXTRNX16C(X0, 3);
  idx1_4 = BBE_EXTRNX16C(X0, 4);
  idx1_5 = BBE_EXTRNX16C(X0, 5);
  idx1_6 = BBE_EXTRNX16C(X0, 6);

  __Pragma("loop_count min=1");
  for (l = 0; l<(L >> 4); l++)
  {
    xb_vecNx16 index0;
    unsigned idx;
    BBE_LVNX16_XP(Y0, w, (MN << 1));
    BBE_LVNX16_XP(Y1, w, (MN << 1));
    BBE_LVNX16_XP(Y2, w, (MN << 1));
    BBE_LVNX16_XP(Y3, w, (MN << 1));
    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y4, w, (MN << 1));
    BBE_LVNX16_XP(Y5, w, (MN << 1));
    BBE_LVNX16_XP(Y6, w, (MN << 1));
    BBE_LVNX16_XP(Y7, w, (MN << 1));
    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(Y8, w, (MN << 1));
    BBE_LVNX16_XP(Y9, w, (MN << 1));
    BBE_LVNX16_XP(YA, w, (MN << 1));
    BBE_LVNX16_XP(YB, w, (MN << 1));
    BBE_DSELNX16I(X9, X8, Y9, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XB, XA, YB, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YA, Y8, XA, X8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YB, Y9, XB, X9, BBE_DSELI_DEINTERLEAVE_1);

    BBE_LVNX16_XP(YC, w, (MN << 1));
    BBE_LVNX16_XP(YD, w, (MN << 1));
    BBE_LVNX16_XP(YE, w, (MN << 1));
    BBE_LVNX16_XP(YF, w, (MN << 1));
    BBE_DSELNX16I(XD, XC, YD, YC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XE, YF, YE, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YE, YC, XE, XC, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(YF, YD, XF, XD, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(X4, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X1, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X6, X2, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X3, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(XC, X8, YC, Y8, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XD, X9, YD, Y9, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XE, XA, YE, YA, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(XF, XB, YF, YB, BBE_DSELI_DEINTERLEAVE_1);

    index0 = BBE_MOVVWL(W0);
    W0 = BBE_MOVWV(index0, index0);
    BBE_DSELNX16I(YF, Y7, XF, X7, BBE_DSELI_DEINTERLEAVE_1);
    idx = BBE_EXTRNX16C(index0, 7);
    BBE_SVNX16_X(Y7, z, idx);
    BBE_DSELNX16I(YE, Y6, XE, X6, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YE, z, idx1_6);
    idx = BBE_EXTRNX16C(index0, 6);
    BBE_SVNX16_X(Y6, z, idx);
    BBE_DSELNX16I(YD, Y5, XD, X5, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YD, z, idx1_5);
    idx = BBE_EXTRNX16C(index0, 5);
    BBE_SVNX16_X(Y5, z, idx);
    BBE_DSELNX16I(YC, Y4, XC, X4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YC, z, idx1_4);
    idx = BBE_EXTRNX16C(index0, 4);
    BBE_SVNX16_X(Y4, z, idx);
    BBE_DSELNX16I(YB, Y3, XB, X3, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YB, z, idx1_3);
    idx = BBE_EXTRNX16C(index0, 3);
    BBE_SVNX16_X(Y3, z, idx);
    BBE_DSELNX16I(YA, Y2, XA, X2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(YA, z, idx1_2);
    idx = BBE_EXTRNX16C(index0, 2);
    BBE_SVNX16_X(Y2, z, idx);
    BBE_DSELNX16I(X9, X1, X9, X1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(X9, z, idx1_1);
    idx = BBE_EXTRNX16C(index0, 1);
    BBE_SVNX16_X(X1, z, idx);
    BBE_DSELNX16I(X8, X0, X8, X0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X(X8, z, idx1_0);
    BBE_SVNX16_IP(X0, z, 2 * BBE_SIMD_WIDTH);
  }
} /* rbsmxn_odd() */

/* Sx=<see the description>, Sy=M*N */
void rbsmxn ( int16_t * restrict y, const int16_t * restrict x, int M, int N, int L )
{
  int MN = M*N;
  /* virtual functions supporting sizes 2...16 */
  typedef void(*fnbs)(int16_t*, const int16_t*, int, int);
  static const fnbs fnsmall[] = { _rbs2x1, rbs3x1, _rbs4x1, rbs5_7x1, rbs5_7x1, rbs5_7x1,
    _rbs8x1, rbs9_15x1, rbs9_15x1, rbs9_15x1, rbs9_15x1,
    rbs9_15x1, rbs9_15x1, rbs9_15x1, _rbs16x1 };
  if (L <= 0 || M <= 0 || N <= 0) return;
  if (MN == 1)  { memcpy(y, x, L * 2); return; }
  fnbs fn;
  fn = (MN <= 16) ? fnsmall[MN - 2] : ((MN & 15) ? rbsmxn_odd : rbsmxn_even);
  fn(y, x, MN, L);
} /* rbsmxn() */

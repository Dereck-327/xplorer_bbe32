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
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, rcmatvmulnxmn,(void* pScr,
                  complex_fract16 * restrict z, 
            const int16_t * restrict x, 
            const complex_fract16 * restrict y, 
            int N, int M, int L, int Q))
size_t rcmatvmulnxmn_getScratchSize(int N, int M, int L) { (void)N; (void)M;  (void)L; return 0; }
#else
/* get allocated space per one real/complex matrix written in the block order */
static int getSpaceC(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, LOG2_BBE_SIMD_WIDTH - 1);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = (((S - 1) >> m) + 1) << m;
  return S;
}
static int getSpaceR(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, LOG2_BBE_SIMD_WIDTH);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = (((S - 1) >> m) + 1) << m;
  return S;
}

/* Block Order, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     N,M,L must be multiples of 4
*/
void rcmatvmulnxmn ( void * pScr,
                     complex_fract16 * restrict z, 
               const int16_t * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int M, int L, int Q )
{
  int l, k, m, n;

  int Sx = getSpaceR(N*M);
  int Sy = getSpaceC(N);
  int Sz = getSpaceC(M);

  int strage;

  int x_k = (M*N << 16) + 16;
  int _k = 0;

  xb_vecNx16 * restrict X_scr = (xb_vecNx16 *)((int32_t)pScr + 0 * (M + 0)*N*L);
  const xb_vecNx16 *          X_rd;
  xb_vecNx16 * restrict X_wr;
  xb_vecNx16 * restrict Y_scr = (xb_vecNx16 *)((int32_t)pScr + 4 * (M + 0)*N*L);
  const xb_vecNx16 *          Y_rd;
  xb_vecNx16 * restrict Y_wr;
  xb_vecNx16 * restrict Z_scr = (xb_vecNx16 *)((int32_t)pScr + 4 * (M + 1)*N*L);
  const xb_vecNx16 *          Z_rd;
  xb_vecNx16 * restrict Z_wr;

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;
  xb_vecNx16 zero = 0;

  xb_vecNx16 X0, X1, X2, X3;
  xb_vecNx16 X0l, X1l, X2l, X3l;
  xb_vecNx16 X0h, X1h, X2h, X3h;
  xb_vecNx16 X0_l, X1_l, X2_l, X3_l;
  xb_vecNx16 X0_h, X1_h, X2_h, X3_h;
  xb_vecNx16 Y0, Y1, Y2, Y3;
  xb_vecNx16 Y0_l, Y1_l, Y2_l, Y3_l, Y_l;
  xb_vecNx16 Y0_h, Y1_h, Y2_h, Y3_h, Y_h;
  xb_vecNx16 Z0_l, Z1_l, Z2_l, Z3_l, Z_l;
  xb_vecNx16 Z0_h, Z1_h, Z2_h, Z3_h, Z_h;

  xb_vecNx40 Acc_l, Acc_h;

  vsaN q = BBE_MOVVSA32(Q);

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L % 4 == 0);
  NASSERT(N % 4 == 0);
  NASSERT(M % 4 == 0);

  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, 2 * Sz*L*sizeof(int16_t));
    return;
  }
  for(k=0; k<Sz*L; k+=8)
  {
    BBE_SVNX16_IP(zero, pz, 2 * BBE_SIMD_WIDTH);
  }
  pz = (xb_vecNx16 *)z;
  /* Convert from block format to streaming format */
  X_wr = X_scr;
  __Pragma("loop_count min=1")
  for (k = 0; k<M*N*L; k += 4 * BBE_SIMD_WIDTH)
  {
    _k = BBE_ADDMOD16U(_k, x_k);
    strage = -3 * 2 * Sx + 2 * BBE_SIMD_WIDTH;
    XT_MOVEQZ(strage, 2 * Sx - 2 * M*N + 2 * BBE_SIMD_WIDTH, _k);

    /* Load input matrix */
    BBE_LVNX16_XP(X0, px, 2 * Sx);
    BBE_LVNX16_XP(X1, px, 2 * Sx);
    BBE_LVNX16_XP(X2, px, 2 * Sx);
    BBE_LVNX16_XP(X3, px, strage);

    /* Convert to streaming */
    BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_1);

    BBE_DSELNX16I(X0_h, X0_l, zero, X0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X1_h, X1_l, zero, X1, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X2_h, X2_l, zero, X2, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(X3_h, X3_l, zero, X3, BBE_DSELI_INTERLEAVE_1);

    /* Save to scr memory */
    BBE_SVNX16_IP(X0_l, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X0_h, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1_l, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1_h, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2_l, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2_h, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X3_l, X_wr, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X3_h, X_wr, 2 * BBE_SIMD_WIDTH);
  }

  Y_wr = Y_scr;
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += 4)
  {
    __Pragma("loop_count min=1")
    for (n = 0; n<N; n += 4)
    {
      /* Load input matrix Y */
      Y0_l = BBE_LV4X16_I(py, 0);
      Y0_h = BBE_LV4X16_I(py, 8);
      py = (const xb_vecNx16 *)XT_ADDX2(2 * Sy, (int32_t)py);
      Y1_l = BBE_LV4X16_I(py, 0);
      Y1_h = BBE_LV4X16_I(py, 8);
      py = (const xb_vecNx16 *)XT_ADDX2(2 * Sy, (int32_t)py);
      Y2_l = BBE_LV4X16_I(py, 0);
      Y2_h = BBE_LV4X16_I(py, 8);
      py = (const xb_vecNx16 *)XT_ADDX2(2 * Sy, (int32_t)py);
      Y3_l = BBE_LV4X16_I(py, 0);
      Y3_h = BBE_LV4X16_I(py, 8);
      py = (const xb_vecNx16 *)XT_ADDX2(-3 * 2 * Sy + 8, (int32_t)py);

      /* Select input matrix Y */
      Y0_l = BBE_SELNX16I(Y2_l, Y0_l, BBE_SELI_INTERLEAVE_2_LO);
      Y1_l = BBE_SELNX16I(Y3_l, Y1_l, BBE_SELI_INTERLEAVE_2_LO);
      Y_l = BBE_SELNX16I(Y1_l, Y0_l, BBE_SELI_INTERLEAVE_2_LO);
      Y0_h = BBE_SELNX16I(Y2_h, Y0_h, BBE_SELI_INTERLEAVE_2_LO);
      Y1_h = BBE_SELNX16I(Y3_h, Y1_h, BBE_SELI_INTERLEAVE_2_LO);
      Y_h = BBE_SELNX16I(Y1_h, Y0_h, BBE_SELI_INTERLEAVE_2_LO);

      /* Save matrix to scr memory */
      BBE_SVNX16_IP(Y_l, Y_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y_h, Y_wr, 2 * BBE_SIMD_WIDTH);
    }

    py = (const xb_vecNx16 *)XT_ADDX2(4 * 2 * Sy - 2 * N, (int32_t)py);
  }

  X_rd = (const xb_vecNx16 *)((int32_t)pScr + 0 * (M + 0)*N*L);
  Y_rd = (const xb_vecNx16 *)((int32_t)pScr + 4 * (M + 0)*N*L);
  Z_wr = Z_scr;
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += 4)
  {
    __Pragma("loop_count min=1")
    for (m = 0; m<M; m += 4)
    {
      /* Load matrix x from scr memory */
      X0_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X0_l, X_rd, 2 * 2 * 4 * N);
      X1_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X1_l, X_rd, 2 * 2 * 4 * N);
      X2_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X2_l, X_rd, 2 * 2 * 4 * N);
      X3_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X3_l, X_rd, -3 * 2 * 2 * 4 * N + 2 * 2 * BBE_SIMD_WIDTH);

      X0l = BBE_SELNX16I(X1_l, X0_l, BBE_SELI_EXTRACT_LO_HALVES);
      X1l = BBE_SELNX16I(X1_l, X0_l, BBE_SELI_EXTRACT_HI_HALVES);
      X2l = BBE_SELNX16I(X1_h, X0_h, BBE_SELI_EXTRACT_LO_HALVES);
      X3l = BBE_SELNX16I(X1_h, X0_h, BBE_SELI_EXTRACT_HI_HALVES);

      X0h = BBE_SELNX16I(X3_l, X2_l, BBE_SELI_EXTRACT_LO_HALVES);
      X1h = BBE_SELNX16I(X3_l, X2_l, BBE_SELI_EXTRACT_HI_HALVES);
      X2h = BBE_SELNX16I(X3_h, X2_h, BBE_SELI_EXTRACT_LO_HALVES);
      X3h = BBE_SELNX16I(X3_h, X2_h, BBE_SELI_EXTRACT_HI_HALVES);

      /* Load input matrix Y */
      BBE_LVNX16_IP(Y_l, Y_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Y_h, Y_rd, 2 * BBE_SIMD_WIDTH);

      Y0 = BBE_SHFLNX16I(Y_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
      Acc_l = BBE_MULRNX16C(X0l, Y0, q);
      Acc_h = BBE_MULRNX16C(X0h, Y0, q);
      Y1 = BBE_SHFLNX16I(Y_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
      BBE_MULANX16C(Acc_l, X1l, Y1);
      BBE_MULANX16C(Acc_h, X1h, Y1);
      Y2 = BBE_SHFLNX16I(Y_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
      BBE_MULANX16C(Acc_l, X2l, Y2);
      BBE_MULANX16C(Acc_h, X2h, Y2);
      Y3 = BBE_SHFLNX16I(Y_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
      BBE_MULANX16C(Acc_l, X3l, Y3);
      BBE_MULANX16C(Acc_h, X3h, Y3);

      for (n = 0; n<N - 4; n += 4)
      {
        /* Load matrix x from scr memory */
        X0_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X0_l, X_rd, 2 * 2 * 4 * N);
        X1_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X1_l, X_rd, 2 * 2 * 4 * N);
        X2_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X2_l, X_rd, 2 * 2 * 4 * N);
        X3_h = BBE_LVNX16_X(X_rd, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(X3_l, X_rd, -3 * 2 * 2 * 4 * N + 2 * 2 * BBE_SIMD_WIDTH);

        X0l = BBE_SELNX16I(X1_l, X0_l, BBE_SELI_EXTRACT_LO_HALVES);
        X1l = BBE_SELNX16I(X1_l, X0_l, BBE_SELI_EXTRACT_HI_HALVES);
        X2l = BBE_SELNX16I(X1_h, X0_h, BBE_SELI_EXTRACT_LO_HALVES);
        X3l = BBE_SELNX16I(X1_h, X0_h, BBE_SELI_EXTRACT_HI_HALVES);

        X0h = BBE_SELNX16I(X3_l, X2_l, BBE_SELI_EXTRACT_LO_HALVES);
        X1h = BBE_SELNX16I(X3_l, X2_l, BBE_SELI_EXTRACT_HI_HALVES);
        X2h = BBE_SELNX16I(X3_h, X2_h, BBE_SELI_EXTRACT_LO_HALVES);
        X3h = BBE_SELNX16I(X3_h, X2_h, BBE_SELI_EXTRACT_HI_HALVES);

        /* Load input matrix Y */
        BBE_LVNX16_IP(Y_l, Y_rd, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y_h, Y_rd, 2 * BBE_SIMD_WIDTH);

        Y0 = BBE_SHFLNX16I(Y_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
        BBE_MULANX16C(Acc_l, X0l, Y0);
        BBE_MULANX16C(Acc_h, X0h, Y0);
        Y1 = BBE_SHFLNX16I(Y_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
        BBE_MULANX16C(Acc_l, X1l, Y1);
        BBE_MULANX16C(Acc_h, X1h, Y1);
        Y2 = BBE_SHFLNX16I(Y_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
        BBE_MULANX16C(Acc_l, X2l, Y2);
        BBE_MULANX16C(Acc_h, X2h, Y2);
        Y3 = BBE_SHFLNX16I(Y_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
        BBE_MULANX16C(Acc_l, X3l, Y3);
        BBE_MULANX16C(Acc_h, X3h, Y3);
      }

      X_rd = (const xb_vecNx16 *)XT_ADDX2(3 * 2 * 4 * N, (int32_t)X_rd);
      Y_rd = (const xb_vecNx16 *)XT_ADDX2(-2 * 4 * N, (int32_t)Y_rd);

      /* Pack and save results */
      Z_l = BBE_PACKVNX40(Acc_l, q);
      Z_h = BBE_PACKVNX40(Acc_h, q);
      BBE_SVNX16_IP(Z_l, Z_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Z_h, Z_wr, 2 * BBE_SIMD_WIDTH);
    }

    Y_rd = (const xb_vecNx16 *)XT_ADDX2(2 * 4 * N, (int32_t)Y_rd);
  }

  /* Convert from streaming format to block format */
  Z_rd = (xb_vecNx16 *)((int32_t)pScr + 4 * (M + 1)*N*L);
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += 4)
  {
    __Pragma("loop_count min=1")
    for (m = 0; m<M; m += 4)
    {
      /* Load matrix Z from scr */
      BBE_LVNX16_IP(Z_l, Z_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Z_h, Z_rd, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Z1_l, Z0_l, Z_l, Z_l, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z2_l, Z0_l, Z0_l, Z0_l, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z3_l, Z1_l, Z1_l, Z1_l, BBE_DSELI_DEINTERLEAVE_2);

      BBE_DSELNX16I(Z1_h, Z0_h, Z_h, Z_h, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z2_h, Z0_h, Z0_h, Z0_h, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z3_h, Z1_h, Z1_h, Z1_h, BBE_DSELI_DEINTERLEAVE_2);

      BBE_SV4X16_I(Z0_l, pz, 0);
      BBE_SV4X16_I(Z0_h, pz, 8);
      pz = (xb_vecNx16 *)XT_ADDX2(2 * Sz, (int32_t)pz);
      BBE_SV4X16_I(Z1_l, pz, 0);
      BBE_SV4X16_I(Z1_h, pz, 8);
      pz = (xb_vecNx16 *)XT_ADDX2(2 * Sz, (int32_t)pz);
      BBE_SV4X16_I(Z2_l, pz, 0);
      BBE_SV4X16_I(Z2_h, pz, 8);
      pz = (xb_vecNx16 *)XT_ADDX2(2 * Sz, (int32_t)pz);
      BBE_SV4X16_I(Z3_l, pz, 0);
      BBE_SV4X16_I(Z3_h, pz, 8);
      pz = (xb_vecNx16 *)XT_ADDX2(-3 * 2 * Sz + 8, (int32_t)pz);
    }

    pz = (xb_vecNx16 *)XT_ADDX2(4 * 2 * Sz - 2 * M, (int32_t)pz);
  }
} /* rcmatvmulnxmn() */

/* Return the scratch area size, in bytes. */
size_t rcmatvmulnxmn_getScratchSize ( int N, int M, int L )
{
  if (M <= 0 || N <= 0 || L <= 0) return 0;
  return (M*N*L+N*L+M*L)*2*sizeof(int16_t);
} /* rcmatvmulnxmn_getScratchSize() */
#endif

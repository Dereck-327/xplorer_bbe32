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
#include "rcmatmulnxmn_common.h"
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, rcmatmulnxmn_gen, (void * pScr,
                complex_fract16 * restrict z,
                const int16_t * restrict x,
                const complex_fract16 * restrict y,
                int N, int M, int L, int Q))
size_t rcmatmulnxmn_gen_getScratchSize(int N, int M) { (void)N, (void)M; return 0; }
#else
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
y[2*L*Sy]   Sequence of right-hand complex matrices or vectors. Real and
            imaginary components are interleaved, with real parts stored at
            even indices.
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[2*L*Sz]   Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

static int getSpaceR(int S)
{
  int m;
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, (LOG2_BBE_SIMD_WIDTH - 0));
  /* round up to the  next multiple of 8 or lesser degree of 2 */
  S = (((S - 1) >> m) + 1) << m;
  return S;
}/*getSpaceR()*/

/* Block Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
      N must be multiples of 4
      M must be multiples of 4
*/
void rcmatmulnxmn_gen ( void * pScr,
                    complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q )
{
  int l, k, m, n;
  xb_vecNx40 W;
  int off1, off2, off3;
  vboolN b;
  vsaN q = BBE_MOVVSA32(Q);

  int Sx = getSpaceR(N*M);
  int Sy = N*M;
  int Sz = M*M;

        xb_vecNx16 * restrict X_scr = (xb_vecNx16 *)((int32_t)pScr + 0 * M*N * 4);
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
        xb_vecNx16 * restrict Y_scr = (xb_vecNx16 *)((int32_t)pScr + 4 * M*N * 4);
  const xb_vecNx16 *          Y_rd;
        xb_vecNx16 * restrict Y_wr;
        xb_vecNx16 * restrict Z_scr = (xb_vecNx16 *)((int32_t)pScr + 8 * M*N * 4);
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;

  const xb_vecNx16 * restrict px = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;

  xb_vecNx16 X0, X1, X2, X3;
  xb_vecNx16 X0_l, X1_l, X2_l, X3_l, X_l;
  xb_vecNx16 X0_h, X1_h, X2_h, X3_h, X_h;
  xb_vecNx16 Y0, Y1, Y2, Y3;
  xb_vecNx16 zero = 0;
  xb_vecNx16 Z0, Z1, Z2, Z3, Z;
  xb_vecNx40 Acc;

  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(N>0 && N % 4 == 0);
  NASSERT(M>0 && M % 4 == 0);
  __Pragma("loop_count min=1")
  for (l = 0; l<L; l += 4)
  {
    b = BBE_LTRN(XT_MIN(4, (L - l)));
    /* Convert to streaming format */
    X0 = BBE_SEQNX16(); X1 = BBE_MOVVA16(2 * Sx);
    X1 = BBE_MOVNX16T(X1, zero, b);
    W = BBE_MULNX16(X0, X1);
    X0 = BBE_MOVVWL(W);
    off1 = BBE_EXTRNX16C(X0, 1);
    off2 = BBE_EXTRNX16C(X0, 2);
    off3 = BBE_EXTRNX16C(X0, 3);
    X_wr = X_scr;
    __Pragma("loop_count min=1")

    for (k = 0; k<M*N; k += BBE_SIMD_WIDTH)
    {
      X1 = BBE_LVNX16_X(px, off1);
      X2 = BBE_LVNX16_X(px, off2);
      X3 = BBE_LVNX16_X(px, off3);
      BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_1);

      BBE_DSELNX16I(X0_h, X0_l, zero, X0, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X1_h, X1_l, zero, X1, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X2_h, X2_l, zero, X2, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(X3_h, X3_l, zero, X3, BBE_DSELI_INTERLEAVE_1);

      BBE_SVNX16_IP(X0_l, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X0_h, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X1_l, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X1_h, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X2_l, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X2_h, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X3_l, X_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X3_h, X_wr, 2 * BBE_SIMD_WIDTH);
    }
    px = (const xb_vecNx16 *)XT_ADDX2(4 * Sx - M*N, (int32_t)px);
    Y_wr = Y_scr;
    /* compute offsets */
    X0 = BBE_SEQNX16(); X1 = BBE_MOVVA16(Sy << 2);
    X1 = BBE_MOVNX16T(X1, zero, b);
    W = BBE_MULNX16(X0, X1);
    X0 = BBE_MOVVWL(W);
    off1 = BBE_EXTRNX16C(X0, 1);
    off2 = BBE_EXTRNX16C(X0, 2);
    off3 = BBE_EXTRNX16C(X0, 3);
    __Pragma("loop_count min=2")
    for (k = 0; k<2 * M*N; k += BBE_SIMD_WIDTH)
    {
      Y1 = BBE_LVNX16_X(py, off1);
      Y2 = BBE_LVNX16_X(py, off2);
      Y3 = BBE_LVNX16_X(py, off3);
      BBE_LVNX16_IP(Y0, py, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Y2, Y0, Y2, Y0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y3, Y1, Y3, Y1, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y1, Y0, Y1, Y0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y3, Y2, Y3, Y2, BBE_DSELI_INTERLEAVE_2);

      BBE_SVNX16_IP(Y0, Y_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y1, Y_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y2, Y_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y3, Y_wr, 2 * BBE_SIMD_WIDTH);
    }
    py = (const xb_vecNx16 *)XT_ADDX4(3 * Sy, (int32_t)py);

    /* Multiple in streaming format */
    X_rd = (const xb_vecNx16 *)pScr;
    Y_rd = (const xb_vecNx16 *)XT_ADDX2(2 * M*N * 4, (int32_t)pScr);
    Z_wr = Z_scr;
    for (k = 0; k<M; k++)
    {
      __Pragma("loop_count min=1")
      for (m = 0; m<M; m += 2)
      {
        /* load input matrix X */
        BBE_LVNX16_IP(X_l, X_rd, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(X_h, X_rd, 2 * BBE_SIMD_WIDTH);

        X0 = BBE_SHFLNX16I(X_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
        X1 = BBE_SHFLNX16I(X_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
        X2 = BBE_SHFLNX16I(X_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
        X3 = BBE_SHFLNX16I(X_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

        /* Load input matrix Y */
        BBE_LVNX16_XP(Y0, Y_rd, 2 * 2 * 4 * M);
        BBE_LVNX16_XP(Y1, Y_rd, 2 * 2 * 4 * M);
        BBE_LVNX16_XP(Y2, Y_rd, 2 * 2 * 4 * M);
        BBE_LVNX16_XP(Y3, Y_rd, 2 * 2 * 4 * M);

        Acc = BBE_MULRNX16C(Y0, X0, q);
        BBE_MULANX16C(Acc, Y1, X1);
        BBE_MULANX16C(Acc, Y2, X2);
        BBE_MULANX16C(Acc, Y3, X3);
        __Pragma("loop_count min=1")
          /* load input matrix X */
          BBE_LVNX16_IP(X_l, X_rd, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_IP(X_h, X_rd, 2 * BBE_SIMD_WIDTH);
          for (n = 0; n<N - 4; n += 4)
          {
           

            X0 = BBE_SHFLNX16I(X_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
            X1 = BBE_SHFLNX16I(X_l, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
            X2 = BBE_SHFLNX16I(X_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
            X3 = BBE_SHFLNX16I(X_h, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

            /* Load input matrix Y */
            BBE_LVNX16_XP(Y0, Y_rd, 2 * 2 * 4 * M);
            BBE_LVNX16_XP(Y1, Y_rd, 2 * 2 * 4 * M);
            BBE_LVNX16_XP(Y2, Y_rd, 2 * 2 * 4 * M);
            BBE_LVNX16_XP(Y3, Y_rd, 2 * 2 * 4 * M);

            BBE_MULANX16C(Acc, Y0, X0);
            BBE_MULANX16C(Acc, Y1, X1);
            BBE_MULANX16C(Acc, Y2, X2);
            BBE_MULANX16C(Acc, Y3, X3);

            /* load input matrix X */
            BBE_LVNX16_IP(X_l, X_rd, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(X_h, X_rd, 2 * BBE_SIMD_WIDTH);
          }

          X_rd = (const xb_vecNx16 *)XT_ADDX2(-2 * 4 * N - 2*BBE_SIMD_WIDTH, (int32_t)X_rd);
        Y_rd = (const xb_vecNx16 *)XT_ADDX2(-2 * 4 * M*N + BBE_SIMD_WIDTH, (int32_t)Y_rd);

        /* Pack and save results */
        Z = BBE_PACKVNX40(Acc, q);
        BBE_SVNX16_IP(Z, Z_wr, 2 * BBE_SIMD_WIDTH);
      }

      X_rd = (const xb_vecNx16 *)XT_ADDX2(2 * 4 * N, (int32_t)X_rd);
      Y_rd = (const xb_vecNx16 *)XT_ADDX2(-2 * 4 * M, (int32_t)Y_rd);
    }

    /* Convert from streaming format to block format */
    Z_rd = (const xb_vecNx16 *)((int32_t)pScr + 8 * M*N * 4);
    X0 = BBE_SEQNX16(); X1 = BBE_MOVVA16(Sz << 2);
    X1 = BBE_MOVNX16T(X1, zero, b);
    W = BBE_MULNX16(X0, X1);
    X0 = BBE_MOVVWL(W);
    off1 = BBE_EXTRNX16C(X0, 1);
    off2 = BBE_EXTRNX16C(X0, 2);
    off3 = BBE_EXTRNX16C(X0, 3);
    __Pragma("loop_count min=8")
    for (k = 0; k<2 * M*M; k += BBE_SIMD_WIDTH)
    {
      BBE_LVNX16_IP(Z0, Z_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Z1, Z_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Z2, Z_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(Z3, Z_rd, 2 * BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Z1, Z0, Z1, Z0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z3, Z2, Z3, Z2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z2, Z0, Z2, Z0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Z3, Z1, Z3, Z1, BBE_DSELI_DEINTERLEAVE_2);

      BBE_SVNX16_X(Z1, pz, off1);
      BBE_SVNX16_X(Z2, pz, off2);
      BBE_SVNX16_X(Z3, pz, off3);
      BBE_SVNX16_IP(Z0, pz, 2 * BBE_SIMD_WIDTH);
    }
    pz = (xb_vecNx16 *)XT_ADDX4(3 * Sz, (int32_t)pz);
  }
} /* rcmatmulnxmn_gen() */

/* Return the scratch area size, in bytes. */
size_t rcmatmulnxmn_gen_getScratchSize ( int N, int M )
{
  int Sx, Sz;
  Sx = getSpaceR(N*M);
  Sz = getSpaceR(M*M);
  return (Sx * 2 + Sz )*2 * 4*sizeof(int16_t);
} /* rcmatmulnxmn_gen_getScratchSize() */
#endif


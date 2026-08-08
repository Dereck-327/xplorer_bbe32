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
 * Matrix Hermitian Product
 * C code optimized for BBE32
 */
/*  
    Optimized code for matrix multiplication
	Integrit, 2006-2016

    specialized function for:
    L - multiple of 8
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "cmathermnxmn_common.h"

/* get allocated space per one matrix (complex) */
static int getSpaceC(int S)
{
    int m;
    /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
    m=30-XT_NSA(S);
    m=XT_MIN(m,3);
    /* round up to the  next multiple of 32 or lesser degree of 2 */
    S=2*((((S-1)>>m)+1)<<m);
    return S;
} /* getSpaceC() */

void cmathermnxmn_L8(void * pScr, complex_fract16 * restrict y, const complex_fract16 * restrict x, int N, int M, int L, int Q)
{
  xb_vecNx16 * restrict Y_scr;
  const xb_vecNx16 *          Y_rd;
  xb_vecNx16 * restrict Y_wr;
  xb_vecNx16 * restrict X_scr;
  const xb_vecNx16 *          X_rd;
  xb_vecNx16 * restrict X_wr;

  int l;

  int MN = M*N, NN = N*N;

  NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);

  NASSERT(!(M & 3) && !(N & 3) && (L>0));

  NASSERT(!(L & (BBE_SIMD_WIDTH / 2 - 1)));

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    Y_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)Y_scr + (BBE_SIMD_WIDTH / 2)*NN * 4);
    X_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)X_scr + (BBE_SIMD_WIDTH / 2)*MN * 4);
  }

  /*
  * Compute BBE_SIMD_WIDTH/2 matrix products at a time.
  */

  for (l = 0; l<L / (BBE_SIMD_WIDTH / 2); l++)
  {
    /*
    * Convert BBE_SIMD_WIDTH/2 MxN complex matrices X to streaming format
    * (out-of-place):
    * x[BBE_SIMD_WIDTH/2][M*N][2] => x[M*N][BBE_SIMD_WIDTH/2][2]
    */

    X_rd = (const xb_vecNx16*)((uintptr_t)x + l*(BBE_SIMD_WIDTH / 2)*MN * 4);
    X_wr = (xb_vecNx16*)X_scr;

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;

      int n;

      __Pragma("loop_count min=1")
      for (n = 0; n<MN / (BBE_SIMD_WIDTH / 2); n++)
      {
        BBE_LVNX16_XP(a0, X_rd, MN * 4);
        BBE_LVNX16_XP(a1, X_rd, MN * 4);
        BBE_LVNX16_XP(a2, X_rd, MN * 4);
        BBE_LVNX16_XP(a3, X_rd, MN * 4);
        BBE_LVNX16_XP(a4, X_rd, MN * 4);
        BBE_LVNX16_XP(a5, X_rd, MN * 4);
        BBE_LVNX16_XP(a6, X_rd, MN * 4);

        BBE_LVNX16_XP(a7, X_rd, -7 * MN * 4 + 4 * BBE_SIMD_WIDTH / 2);

        BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a5, a4, a5, a4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a6, a7, a6, BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a6, a4, a6, a4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a5, a7, a5, BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(a4, a0, a4, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a5, a1, a5, a1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a6, a2, a6, a2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a3, a7, a3, BBE_DSELI_DEINTERLEAVE_2);

        BBE_SVNX16_IP(a0, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a1, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a2, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a3, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a4, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a5, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a6, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a7, X_wr, 4 * BBE_SIMD_WIDTH / 2);
      }
    }

    /*
    * Each input matrix X is to be multiplied on the left by its conjugate
    * transpose: Y = conj(X')*X. Compute BBE_SIMD_WIDTH/2 matrix products in
    * streaming format.
    */

    {
      const xb_vecNx16 * X0;
      const xb_vecNx16 * X1;
      xb_vecNx16 * Y0;
      xb_vecNx16 * Y1;

      xb_vecNx16 x00, x01, x10, x11;
      xb_vecNx16 y00, y01, y02, y03;
      xb_vecNx16 y10, y11, y12, y13;

      xb_vecNx40 w0, w1, w2, w3;

      vsaN vs0;

      int i, j, k;

      vs0 = BBE_MOVVSA32(Q);

      __Pragma("loop_count min=1");
      for (i = 0; i<N / 2; i++)
      {
        __Pragma("loop_count min=1");
        for (j = i; j<N / 2; j++)
        {
          X0 = (const xb_vecNx16*)((uintptr_t)X_scr + 2 * i*(BBE_SIMD_WIDTH / 2) * 4);
          X1 = (const xb_vecNx16*)((uintptr_t)X_scr + 2 * j*(BBE_SIMD_WIDTH / 2) * 4);

          Y0 = (xb_vecNx16*)((uintptr_t)Y_scr + 2 * (i*N + j)*(BBE_SIMD_WIDTH / 2) * 4);
          Y1 = (xb_vecNx16*)((uintptr_t)Y_scr + 2 * (j*N + i)*(BBE_SIMD_WIDTH / 2) * 4);

          {
            x01 = BBE_LVNX16_I(X0, 4 * BBE_SIMD_WIDTH / 2);
            BBE_LVNX16_XP(x00, X0, N * 4 * BBE_SIMD_WIDTH / 2);

            x11 = BBE_LVNX16_I(X1, 4 * BBE_SIMD_WIDTH / 2);
            BBE_LVNX16_XP(x10, X1, N * 4 * BBE_SIMD_WIDTH / 2);

            w0 = BBE_MULRNX16J(x10, x00, vs0);
            w1 = BBE_MULRNX16J(x11, x00, vs0);
            w2 = BBE_MULRNX16J(x10, x01, vs0);
            w3 = BBE_MULRNX16J(x11, x01, vs0);
          }

          __Pragma("ymemory( X0 )");
          __Pragma("ymemory( X1 )");
          __Pragma("loop_count min=1");
          for (k = 0; k<M - 1; k++)
          {
            x01 = BBE_LVNX16_I(X0, 4 * BBE_SIMD_WIDTH / 2);
            BBE_LVNX16_XP(x00, X0, N * 4 * BBE_SIMD_WIDTH / 2);

            x11 = BBE_LVNX16_I(X1, 4 * BBE_SIMD_WIDTH / 2);
            BBE_LVNX16_XP(x10, X1, N * 4 * BBE_SIMD_WIDTH / 2);

            BBE_MULANX16J(w0, x10, x00);
            BBE_MULANX16J(w1, x11, x00);
            BBE_MULANX16J(w2, x10, x01);
            BBE_MULANX16J(w3, x11, x01);
          }

          y00 = BBE_PACKVNX40(w0, vs0);
          y01 = BBE_PACKVNX40(w1, vs0);
          y02 = BBE_PACKVNX40(w2, vs0);
          y03 = BBE_PACKVNX40(w3, vs0);

          y10 = BBE_CONJSNX16C(y00);
          y11 = BBE_CONJSNX16C(y02);
          y12 = BBE_CONJSNX16C(y01);
          y13 = BBE_CONJSNX16C(y03);

          BBE_SVNX16_X(y02, Y0, N * 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_IP(y00, Y0, 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_X(y03, Y0, N * 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_IP(y01, Y0, 4 * BBE_SIMD_WIDTH / 2);

          BBE_SVNX16_I(y11, Y1, 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_XP(y10, Y1, N * 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_I(y13, Y1, 4 * BBE_SIMD_WIDTH / 2);
          BBE_SVNX16_XP(y12, Y1, N * 4 * BBE_SIMD_WIDTH / 2);
        }
      }
    }

    /*
    * Convert BBE_SIMD_WIDTH/2 NxN complex matrices Y to block format
    * (out-of-place):
    * y[N*N][BBE_SIMD_WIDTH/2][2] => y[BBE_SIMD_WIDTH/2][N*N][2].
    */

    Y_rd = Y_scr;
    Y_wr = (xb_vecNx16*)((uintptr_t)y + l*(BBE_SIMD_WIDTH / 2)*NN * 4);

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;

      int n;

      __Pragma("loop_count min=1")
      for (n = 0; n<NN / (BBE_SIMD_WIDTH / 2); n++)
      {
        BBE_LVNX16_IP(a0, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a1, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a2, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a3, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a4, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a5, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a6, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a7, Y_rd, 4 * BBE_SIMD_WIDTH / 2);

        BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a5, a4, a5, a4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a6, a7, a6, BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a6, a4, a6, a4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a5, a7, a5, BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(a4, a0, a4, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a5, a1, a5, a1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a6, a2, a6, a2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(a7, a3, a7, a3, BBE_DSELI_DEINTERLEAVE_2);

        BBE_SVNX16_XP(a0, Y_wr, NN * 4);
        BBE_SVNX16_XP(a1, Y_wr, NN * 4);
        BBE_SVNX16_XP(a2, Y_wr, NN * 4);
        BBE_SVNX16_XP(a3, Y_wr, NN * 4);
        BBE_SVNX16_XP(a4, Y_wr, NN * 4);
        BBE_SVNX16_XP(a5, Y_wr, NN * 4);
        BBE_SVNX16_XP(a6, Y_wr, NN * 4);

        BBE_SVNX16_XP(a7, Y_wr, -7 * NN * 4 + 4 * BBE_SIMD_WIDTH / 2);
      }
    }
  }
} /* cmathermnxmn_L8() */
/* return scratch size for L8 family functions */
size_t cmathermnxmn_L8_getScratchSize(int N, int M)
{
    int Sx,Sz;
    Sx = getSpaceC(N*M);
    Sz = getSpaceC(N*N);
    return (Sx*8+Sz*8)*sizeof(int16_t);
} /* cmathermnxmn_L8_getScratchSize() */


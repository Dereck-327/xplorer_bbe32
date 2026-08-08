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
    L = 2
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

void cmathermnxmn_L2(void * pScr, complex_fract16 * restrict y, const complex_fract16 * restrict x, int N, int M, int L, int Q)
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

  NASSERT(!(L & (BBE_SIMD_WIDTH / 8 - 1)));

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    Y_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)Y_scr + (BBE_SIMD_WIDTH / 8)*NN * 4);
    X_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)X_scr + (BBE_SIMD_WIDTH / 8)*MN * 4);
  }

  /*
  * Compute BBE_SIMD_WIDTH/8 matrix products at a time.
  */

  for (l = 0; l<L / (BBE_SIMD_WIDTH / 8); l++)
  {
    /*
    * Convert BBE_SIMD_WIDTH/8 complex MxN matrices to streaming format: 
    * x[BBE_SIMD_WIDTH/8][M*N][2] => x[M*N][BBE_SIMD_WIDTH/8][2]
    */

    X_rd = (const xb_vecNx16*)((uintptr_t)x + l*(BBE_SIMD_WIDTH / 8)*MN * 4);
    X_wr = X_scr;

    {
      xb_vecNx16 a0, a1;

      int n;

      __Pragma("loop_count min=2");
      for (n = 0; n<MN / (BBE_SIMD_WIDTH / 2); n++)
      {
        a1 = BBE_LVNX16_X(X_rd, 4 * MN);
        BBE_LVNX16_IP(a0, X_rd, 4 * BBE_SIMD_WIDTH / 2);

        BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_2);

        BBE_SVNX16_IP(a0, X_wr, 4 * BBE_SIMD_WIDTH / 2);
        BBE_SVNX16_IP(a1, X_wr, 4 * BBE_SIMD_WIDTH / 2);
      }
    }

    /*
    * Each input matrix X is to be multiplied on the left by its conjugate
    * transpose: Y = conj(X')*X. Compute BBE_SIMD_WIDTH/8 matrix products in
    * streaming format.
    */

    {
      const xb_vecNx16 * X0;
      const xb_vecNx16 * X1;
      xb_vecNx16 * Y0;
      xb_vecNx16 * Y1;

      xb_vecNx16 x0, x1;
      xb_vecNx16 y0, y1, y2, y3;

      xb_vecNx16 x10, x11, x12, x13;
      xb_vecNx16 y00, y01, y02, y03;
      xb_vecNx16 y10, y11, y12, y13;

      xb_vecNx40 w0, w1, w2, w3;

      vsaN vsa0;

      int i, j, k;

      vsa0 = BBE_MOVVSA32(Q);

      __Pragma("loop_count min=1");
      for (i = 0; i<N / 4; i++)
      {
        __Pragma("loop_count min=1");
        for (j = i; j<N / 4; j++)
        {
          X0 = (const xb_vecNx16*)((uintptr_t)X_scr + 4 * i*(BBE_SIMD_WIDTH / 8) * 4);
          X1 = (const xb_vecNx16*)((uintptr_t)X_scr + 4 * j*(BBE_SIMD_WIDTH / 8) * 4);

          Y0 = (xb_vecNx16*)((uintptr_t)Y_scr + 4 * (i*N + j)*(BBE_SIMD_WIDTH / 8) * 4);
          Y1 = (xb_vecNx16*)((uintptr_t)Y_scr + 4 * (j*N + i)*(BBE_SIMD_WIDTH / 8) * 4);

          {
            BBE_LVNX16_XP(x0, X0, 4 * (BBE_SIMD_WIDTH / 8)*N);
            BBE_LVNX16_XP(x1, X1, 4 * (BBE_SIMD_WIDTH / 8)*N);

            x10 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_0X4);
            x11 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_1X4);
            x12 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_2X4);
            x13 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_3X4);

            w0 = BBE_MULRNX16J(x10, x0, vsa0);
            w1 = BBE_MULRNX16J(x11, x0, vsa0);
            w2 = BBE_MULRNX16J(x12, x0, vsa0);
            w3 = BBE_MULRNX16J(x13, x0, vsa0);

            BBE_LVNX16_XP(x0, X0, 4 * (BBE_SIMD_WIDTH / 8)*N);
            BBE_LVNX16_XP(x1, X1, 4 * (BBE_SIMD_WIDTH / 8)*N);

            x10 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_0X4);
            x11 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_1X4);
            x12 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_2X4);
            x13 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_3X4);

            BBE_MULANX16J(w0, x10, x0);
            BBE_MULANX16J(w1, x11, x0);
            BBE_MULANX16J(w2, x12, x0);
            BBE_MULANX16J(w3, x13, x0);
          }

          __Pragma("ymemory( X0 )");
          __Pragma("ymemory( X1 )");
          __Pragma("loop_count min=1");
          for (k = 0; k<M / 2 - 1; k++)
          {
            BBE_LVNX16_XP(x0, X0, 4 * (BBE_SIMD_WIDTH / 8)*N);
            BBE_LVNX16_XP(x1, X1, 4 * (BBE_SIMD_WIDTH / 8)*N);

            x10 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_0X4);
            x11 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_1X4);
            x12 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_2X4);
            x13 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_3X4);

            BBE_MULANX16J(w0, x10, x0);
            BBE_MULANX16J(w1, x11, x0);
            BBE_MULANX16J(w2, x12, x0);
            BBE_MULANX16J(w3, x13, x0);

            BBE_LVNX16_XP(x0, X0, 4 * (BBE_SIMD_WIDTH / 8)*N);
            BBE_LVNX16_XP(x1, X1, 4 * (BBE_SIMD_WIDTH / 8)*N);

            x10 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_0X4);
            x11 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_1X4);
            x12 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_2X4);
            x13 = BBE_SHFLNX16I(x1, BBE_SHFLI_REP_3X4);

            BBE_MULANX16J(w0, x10, x0);
            BBE_MULANX16J(w1, x11, x0);
            BBE_MULANX16J(w2, x12, x0);
            BBE_MULANX16J(w3, x13, x0);
          }

          y0 = BBE_PACKVNX40(w0, vsa0);
          y1 = BBE_PACKVNX40(w1, vsa0);
          y2 = BBE_PACKVNX40(w2, vsa0);
          y3 = BBE_PACKVNX40(w3, vsa0);

          BBE_DSELNX16I(y02, y00, y2, y0, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELNX16I(y03, y01, y3, y1, BBE_DSELI_INTERLEAVE_4);

          BBE_DSELNX16I(y01, y00, y01, y00, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELNX16I(y03, y02, y03, y02, BBE_DSELI_INTERLEAVE_4);

          BBE_SVNX16_XP(y00, Y0, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y01, Y0, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y02, Y0, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y03, Y0, 4 * (BBE_SIMD_WIDTH / 8)*(4 - 3 * N));

          y10 = BBE_CONJSNX16C(y0);
          y11 = BBE_CONJSNX16C(y1);
          y12 = BBE_CONJSNX16C(y2);
          y13 = BBE_CONJSNX16C(y3);

          BBE_SVNX16_XP(y10, Y1, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y11, Y1, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y12, Y1, 4 * (BBE_SIMD_WIDTH / 8)*N);
          BBE_SVNX16_XP(y13, Y1, 4 * (BBE_SIMD_WIDTH / 8)*N);
        }
      }
    }

    /*
    * Convert 4 complex NxN matrices to block format:
    * y[N*N][L][2] => y[L][N*N][2]
    */

    Y_rd = Y_scr;
    Y_wr = (xb_vecNx16*)((uintptr_t)y + l*(BBE_SIMD_WIDTH / 8)*NN * 4);

    {
      xb_vecNx16 a0, a1;

      int n;

      for (n = 0; n<NN / (BBE_SIMD_WIDTH / 2); n++)
      {
        BBE_LVNX16_IP(a0, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
        BBE_LVNX16_IP(a1, Y_rd, 4 * BBE_SIMD_WIDTH / 2);

        BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);

        BBE_SVNX16_X(a1, Y_wr, 4 * NN);
        BBE_SVNX16_IP(a0, Y_wr, 4 * BBE_SIMD_WIDTH / 2);
      }
    }
  }
} /* cmathermnxmn_L2() */
/* return scratch size for L2 family functions */
size_t cmathermnxmn_L2_getScratchSize(int N, int M)
{
    int Sx,Sz;
    Sx = getSpaceC(N*M);
    Sz = getSpaceC(N*N);
    return (Sx*2+Sz*2)*sizeof(int16_t);
} /* cmathermnxmn_L2_getScratchSize() */


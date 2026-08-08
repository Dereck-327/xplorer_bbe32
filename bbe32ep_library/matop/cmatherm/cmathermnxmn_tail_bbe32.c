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
    L = 1
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "cmathermnxmn_common.h"


#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmathermnxmn_tail, (void * pScr,
                   complex_fract16 * restrict y, 
             const complex_fract16 * restrict x, 
             int N, int M, int Q ))
size_t cmathermnxmn_tail_getScratchSize(int N, int M) { (void)N; (void)M; return 0; }
#else
/* get allocated space per one matrix (real) */
static int getSpaceC(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, 3);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = 2 * ((((S - 1) >> m) + 1) << m);
  return S;
} /* getSpaceC() */

void cmathermnxmn_tail(void* pScr, complex_fract16* y, const complex_fract16* x, int N, int M, int Q)
{
  xb_vecNx16 * restrict Y_scr;
  const xb_vecNx16 *          Y_rd;
  xb_vecNx16 * restrict Y_wr;
  xb_vecNx16 * restrict X_scr;
  const xb_vecNx16 *          X_rd;
  xb_vecNx16 * restrict X_wr;

  int MN = M*N, NN = N*N;

  NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);

  NASSERT(!(M & 3) && !(N & 3));

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    Y_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)Y_scr + 1 * NN * 4);
    X_scr = (xb_vecNx16*)ptr;
    ptr = (void*)((uintptr_t)X_scr + 1 * MN * 4);
  }

  /*
  * Interleave elements of upper and lower halves of the input matrix. This
  * is done to make the offset between any two rowa a multiple of vector
  * length. */

  X_rd = (const xb_vecNx16*)x;
  X_wr = (xb_vecNx16*)X_scr;

  {
    xb_vecNx16 a0, a1;

    int n;

    __Pragma("loop_count min=1");
    for (n = 0; n<(MN / 2) / (BBE_SIMD_WIDTH / 2); n++)
    {
      a1 = BBE_LVNX16_X(X_rd, MN / 2 * 4);
      BBE_LVNX16_IP(a0, X_rd, 4 * BBE_SIMD_WIDTH / 2);

      BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_2);

      BBE_SVNX16_IP(a0, X_wr, 4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(a1, X_wr, 4 * BBE_SIMD_WIDTH / 2);
    }
  }

  __Pragma("no_reorder");

  /*
  * Multiply the input matrix by its conjugate transpose: Y = conj(X')*X.
  * Resulting data are stored to the scratch memory. */

  {
    const xb_vecNx16 * X0;
    const xb_vecNx16 * X1;
    xb_vecNx16 * Y0;
    xb_vecNx16 * Y1;

    xb_vecNx40 w0, w1;
    xb_vecNx16 x0, x1, t0, t1;
    xb_vecNx16 x00, x01, x02, x03;
    xb_vecNx16 y00, y01, y10, y11;

    static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
      { 0, 1, 8, 9, 16, 17, 24, 25, 2, 3, 10, 11, 18, 19, 26, 27 },
      { 4, 5, 12, 13, 20, 21, 28, 29, 6, 7, 14, 15, 22, 23, 30, 31 }
    };

    vselN sel0, sel1;
    vsaN  vsa0;

    int i, j, k;

    t0 = BBE_LVNX16_I((const xb_vecNx16*)cst, 0 * 2 * BBE_SIMD_WIDTH);
    t1 = BBE_LVNX16_I((const xb_vecNx16*)cst, 1 * 2 * BBE_SIMD_WIDTH);

    sel0 = BBE_MOVVSV(t0, 0);
    sel1 = BBE_MOVVSV(t1, 0);

    vsa0 = BBE_MOVVSA32(Q);

    for (i = 0; i<N / 4; i++)
    {
      __Pragma("loop_count min=1");
      for (j = i; j<N / 4; j++)
      {
        X0 = (const xb_vecNx16*)((uintptr_t)X_scr + 4 * i * 2 * 4);
        X1 = (const xb_vecNx16*)((uintptr_t)X_scr + 4 * j * 2 * 4);

        Y0 = (xb_vecNx16*)((uintptr_t)Y_scr + 4 * (i*N + j * 2) * 4);
        Y1 = (xb_vecNx16*)((uintptr_t)Y_scr + 4 * (j*N + i * 2) * 4);

        w0 = w1 = BBE_MOVWA32((1UL << Q) >> 1);

        __Pragma("loop_count min=2");
        for (k = 0; k<M / 2; k++)
        {
          BBE_LVNX16_XP(x0, X0, 2 * N * 4);
          BBE_LVNX16_XP(x1, X1, 2 * N * 4);

          x1 = BBE_CONJSNX16C(x1);

          x00 = BBE_SHFLNX16I(x0, BBE_SHFLI_REP_0X4);
          x01 = BBE_SHFLNX16I(x0, BBE_SHFLI_REP_1X4);
          x02 = BBE_SHFLNX16I(x0, BBE_SHFLI_REP_2X4);
          x03 = BBE_SHFLNX16I(x0, BBE_SHFLI_REP_3X4);

          BBE_MULANX16PC_0(w0, x00, x1);
          BBE_MULANX16PC_0(w1, x01, x1);
          BBE_MULANX16PC_1(w0, x02, x1);
          BBE_MULANX16PC_1(w1, x03, x1);
        }

        t0 = BBE_PACKVNX40(w0, vsa0);
        t1 = BBE_PACKVNX40(w1, vsa0);

        y00 = BBE_CONJSNX16C(t0);
        y01 = BBE_CONJSNX16C(t1);

        /* 4x4 sub-block above the diagonal. */
        BBE_SVNX16_XP(y00, Y0, 2 * N * 4);
        BBE_SVNX16_XP(y01, Y0, 2 * N * 4);

        y10 = BBE_SELNX16(t1, t0, sel0);
        y11 = BBE_SELNX16(t1, t0, sel1);

        /* 4x4 sub-block below the diagonal. */
        BBE_SVNX16_XP(y10, Y1, 2 * N * 4);
        BBE_SVNX16_XP(y11, Y1, 2 * N * 4);
      }
    }
  }

  __Pragma("no_reorder");

  /*
  * Product matrix appears in the following format: Y[N/4][2*N][2][2]. That
  * is, for each group of 4 rows 2*N elements of the first two rows are
  * interleaved with 2*N elements of the second two rows. We have to copy data
  * to the output array in natural order: Y[N][N][2]. */

  Y_rd = (const xb_vecNx16*)Y_scr;
  Y_wr = (xb_vecNx16*)y;

  {
    xb_vecNx16 a0, a1;

    uint32_t ix, ix_upd;
    int    Y_upd;

    int i;

    ix = 0;

    ix_upd = ((uint32_t)((2 * N) / (BBE_SIMD_WIDTH / 2)) << 16 | 1);

    __Pragma("loop_count min=1");
    for (i = 0; i<(N*N / 2) / (BBE_SIMD_WIDTH / 2); i++)
    {
      BBE_LVNX16_IP(a0, Y_rd, 4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(a1, Y_rd, 4 * BBE_SIMD_WIDTH / 2);

      BBE_DSELNX16I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);

      ix = BBE_ADDMOD16U(ix, ix_upd);

      Y_upd = 4 * BBE_SIMD_WIDTH / 2;

      /*If we finish deinterleaving for the current block of 4 rows
       (4*N elements), then jump to the next block. Otherwise go to the
       next vectors of the current block. */
      XT_MOVEQZ(Y_upd, 2 * N * 4 + 4 * BBE_SIMD_WIDTH / 2, ix);

      BBE_SVNX16_X(a1, Y_wr, 2 * N * 4); 
      BBE_SVNX16_XP(a0, Y_wr, Y_upd); 
    }
  }
} /* cmathermnxmn_tail() */

/* Return the scratch area size, in bytes. */
size_t cmathermnxmn_tail_getScratchSize(int N, int M)
{
  int Sx, Sz;
  Sx = getSpaceC(N*M);
  Sz = getSpaceC(N*N);
  return (Sx * 1 + Sz * 1)*sizeof(int16_t);
} /* cmathermnxmn_tail_getScratchSize() */
#endif

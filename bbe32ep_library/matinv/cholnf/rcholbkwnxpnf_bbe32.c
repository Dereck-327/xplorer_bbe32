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
    Cholesky backward recursion, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
    */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "rcholnf_common.h"

#if (HAVE_VFPU)

#define VECLEN (BBE_SIMD_WIDTH/2)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-1)) m=LOG2_BBE_SIMD_WIDTH-1;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}
/*
    backward recursion: P!=1
    */
static void rbkwnxpf(
          float32_t* restrict x,
    const float32_t* restrict Rt,
    const float32_t* restrict D,
    const float32_t* restrict y,
    int N, int P, int L)
{
#if 0
    int m, k, p;
    int l;
    int SX = getSpace(N*P);
    int SD = getSpace(N);

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    for (k = N - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t* pRt = Rt + l*(N*(N - 1) >> 1) + (((N - k - 1)*(N - k - 2)) >> 1);
            // calculate y(m,:)-R(m,:)*X, 1xP
            for (p = 0; p < P; p++)
            {
                float32_t r_re;
                float32_t x_re;
                float32_t B_re;
                B_re = (y[l*SX + k*P + p]);

                for (m = 0; m < N - k - 1; m++)
                {
                    x_re = x[l*SX + (k + 1 + m)*P + p];
                    r_re = pRt[m + 0];
                    B_re -= (x_re*r_re);
                }
                x[l*SX + k*P + p] = B_re * D[l*SD + k];
            }
        }
    }
#endif // 0

    int m, k, p, l;
    int SX = getSpace(N*P);
    int SD = getSpace(N);

    const xb_vecN_2xf32 * restrict pY;
          xb_vecN_2xf32 * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX0;
    const xb_vecN_2xf32 * restrict pX_;
    const xtfloat       * restrict pR;
    const xtfloat       * restrict pD;

    valign vY, vX, vXw;
    xb_vecN_2xf32 Acc, Acc1, X0, X1, R0, D0;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    for (l = 0; l < L; l++)
    {
        pD = (const xtfloat *)(D + l*SD + (N - 1));
        for (k = N - 1; k >= 0; k--)
        {
            pX_ = (const xb_vecN_2xf32 *)(x + l*SX + (k + 1)*P);

            pY = (const xb_vecN_2xf32 *)(y + l*SX + k*P);
            pXw = (xb_vecN_2xf32 *)(x + l*SX + k*P);
            vY = BBE_LAN_2XF32_PP(pY);
            vXw = BBE_ZALIGN();

            BBE_LSN_2XF32_XP(D0, pD, -4);
            D0 = BBE_REPN_2XF32(D0, 0);

            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pR = (const xtfloat *)(Rt + l*(N*(N - 1) >> 1) + (((N - k - 1)*(N - k - 2)) >> 1));
                pX = pX_;

                BBE_LAVN_2XF32_XP(Acc, vY, pY, p * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(Acc1, vY, pY, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));

                for (m = 0; m < (N - k - 1); m++)
                {
                    pX0 = pX;
                    vX = BBE_LAN_2XF32_PP(pX0);
                    BBE_LAVN_2XF32_XP(X0, vX, pX0, 2 * p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(X1, vX, pX0, 2 * (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
                    BBE_LSN_2XF32_XP(R0, pR, sizeof(float32_t));
                    R0 = BBE_REPN_2XF32(R0, 0);

                    BBE_MULSN_2XF32(Acc, X0, R0);
                    BBE_MULSN_2XF32(Acc1, X1, R0);

                    pX = (const xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pX);
                }

                Acc = BBE_MULN_2XF32(Acc, D0);
                Acc1 = BBE_MULN_2XF32(Acc1, D0);
                BBE_SAVN_2XF32_XP(Acc, vXw, pXw, p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc1, vXw, pXw, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));

                pX_ = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH, (uintptr_t)pX_);
            }
            BBE_SAN_2XF32POS_FP(vXw, pXw);
        }
    }
}


/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SR,SD,SY,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SY = storage_size(N*P)
SX = storage_size(N*P)

Scratch size in bytes is defined by [r]cholbkw<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 N         Matrix dimension (number of columns and rows in matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal
 y[L][SY]  Sequence of intermediate decision matrices y
Output:         
 x[L][SX]  Sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4
---------------------------------------------------------------------------*/
void rcholbkwnxpnf(
            void * pScr,
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int N, int P, int L )
{
    float32_t* Rt = (float32_t*)pScr;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0 && P>0);
    if (L<0) return;

    rcholnfTransformR(Rt, R, N, L);
    if (P != 1)
    {
        rbkwnxpf(x, Rt, D, y, N, P, L);
    }
    else
    {
        rcholnfBkwnx1(x, Rt, D, y, N, L);
    }
}

/*
    return Scratch size
    Input:
    N:      matrix size
    P:      ignored
    L       number of matrices
    */
size_t rcholbkwnxpnf_getScratchSize(int N, int P, int L)
{
    L=XT_MAX(0,L);
    NASSERT(N>0 && N%4==0 && P>0);
    return (L*((N*(N - 1))>>1)*sizeof(float32_t));
}

#else
DISCARD_FUN(void, rcholbkwnxpnf,(
            void * pScr,
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int N, int P, int L ))

size_t rcholbkwnxpnf_getScratchSize(int N, int P, int L)
{
  (void)N; (void)P; (void)L;
  return 0;
}

#endif

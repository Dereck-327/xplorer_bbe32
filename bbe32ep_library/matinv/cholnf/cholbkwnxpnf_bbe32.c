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
#include "cholnf_common.h"

#if (HAVE_VFPU)

#define VECLEN (BBE_SIMD_WIDTH/4)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-2)) m=LOG2_BBE_SIMD_WIDTH-2;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}
/*
    backward recursion: P!=1
    */
static void bkwnxpf(
          float32_t* restrict x,
    const float32_t* restrict Rt,
    const float32_t* restrict D,
    const float32_t* restrict y,
    int N, int P, int L)
{
#if 0
    int m, k, p;
    int l;
    int SX = 2 * getSpace(N*P);
    int SD = 2 * getSpace(N);

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    for (k = N - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t* pRt = Rt + l*N*(N - 1) + (N - k - 1)*(N - k - 2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            for (p = 0; p < P; p++)
            {
                float32_t r_re, r_im;
                float32_t x_re, x_im;
                float32_t B_re, B_im;
                B_re = (y[l*SX + 2 * k*P + p * 2 + 0]);
                B_im = (y[l*SX + 2 * k*P + p * 2 + 1]);

                for (m = 0; m < N - k - 1; m++)
                {
                    x_re = x[l*SX + 2 * (k + 1 + m)*P + 2 * p + 0];
                    x_im = x[l*SX + 2 * (k + 1 + m)*P + 2 * p + 1];
                    r_re = pRt[2 * m + 0];
                    r_im = pRt[2 * m + 1];
                    B_re -= (x_re*r_re) - (x_im*r_im);
                    B_im -= (x_re*r_im) + (x_im*r_re);
                }
                x[l*SX + 2 * k*P + 2 * p + 0] = B_re * D[l*SD + 2 * k + 0];
                x[l*SX + 2 * k*P + 2 * p + 1] = B_im * D[l*SD + 2 * k + 1];
            }
        }
    }
#endif // 0

    int m, k, p, l;
    int SX = 2 * getSpace(N*P);
    int SD = 2 * getSpace(N);

    const xb_vecN_2xf32 * restrict pY;
          xb_vecN_2xf32 * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX0;
    const xb_vecN_2xf32 * restrict pX_;
    const long long     * restrict pR;
    const long long     * restrict pD;

    valign vY, vX, vXw;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, X0, X1, R0, D0;
    xb_vecN_4x64 temp;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    for (l = 0; l < L; l++)
    {
        pD = (const long long *)(D + l*SD + 2 * (N - 1));
        for (k = N - 1; k >= 0; k--)
        {
            pX_ = (const xb_vecN_2xf32 *)(x + l*SX + 2 * (k + 1)*P);

            pY = (const xb_vecN_2xf32 *)(y + l*SX + 2 * k*P);
            pXw = (xb_vecN_2xf32 *)(x + l*SX + 2 * k*P);
            vY = BBE_LAN_2XF32_PP(pY);
            vXw = BBE_ZALIGN();

            BBE_LSN_4X64_XP(temp, pD, -2 * 4);
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            D0 = BBE_SHFLN_2XF32I(D0, BBE_SHFLI_REP_0X4);

            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pR = (const long long *)(Rt + l*N*(N - 1) + (N - k - 1)*(N - k - 2));
                pX = pX_;

                BBE_LAVN_2XF32_XP(Acc, vY, pY, 2 * p * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(Acc2, vY, pY, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
                Acc1 = BBE_ZERON_2XF32();
                Acc3 = BBE_ZERON_2XF32();

                for (m = 0; m < (N - k - 1); m++)
                {
                    pX0 = pX;
                    vX = BBE_LAN_2XF32_PP(pX0);
                    BBE_LAVN_2XF32_XP(X0, vX, pX0, 2 * p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(X1, vX, pX0, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
                    BBE_LSN_4X64_XP(temp, pR, 2 * sizeof(float32_t));
                    R0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
                    R0 = BBE_SHFLN_2XF32I(R0, BBE_SHFLI_REP_0X4);

                    BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
                    BBE_MULMASN_2XF32(Acc1, X0, R0, 2, 11);

                    BBE_MULMASN_2XF32(Acc2, X1, R0, 3, 4);
                    BBE_MULMASN_2XF32(Acc3, X1, R0, 2, 11);

                    pX = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pX);
                }
                Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
                Acc = BBE_ADDN_2XF32(Acc, Acc1);

                Acc = BBE_MULN_2XF32(Acc, D0);
                Acc2 = BBE_MULN_2XF32(Acc2, D0);
                BBE_SAVN_2XF32_XP(Acc, vXw, pXw, 2 * p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc2, vXw, pXw, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));

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
void cholbkwnxpnf(
            void * pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
      int N, int P, int L )
{
          float32_t* restrict x=(      float32_t*)_x;
    const float32_t* restrict R=(const float32_t*)_R;
    const float32_t* restrict D=(const float32_t*)_D;
    const float32_t* restrict y=(const float32_t*)_y;
    float32_t* Rt = (float32_t*)pScr;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    cholnfTransformR(Rt, R, N, L);
    if (P != 1)
    {
        bkwnxpf(x, Rt, D, y, N, P, L);
    }
    else
    {
        cholnfBkwnx1(x, Rt,D,y, N,L);
    }
}

/*
    return Scratch size
    Input:
    N:      matrix size
    P:      ignored
    L       number of matrices
    */
size_t cholbkwnxpnf_getScratchSize(int N, int P, int L)
{
    return (L*N*(N - 1)*sizeof(float32_t));
}

#else
DISCARD_FUN(void, cholbkwnxpnf,(
            void * pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
      int N, int P, int L ))

size_t cholbkwnxpnf_getScratchSize(int N, int P, int L)
{
  (void)N; (void)P; (void)L;
  return 0;
}

#endif

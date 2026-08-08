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
    Cholesky forward recursion, floating point real data, block format
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
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void rcomputeABf(float32_t* Z,const float32_t* A,const float32_t* B, int M,int N,int P,int L)
{
#if 0
    float32_t B_re, B_im;
    int n, m, l, p;
    int SA = getSpace(M*N);
    int SB = getSpace(M*P);

    for (p = 0; p < P; p++)
        for (n = 0; n < N; n++)
            for (l = 0; l < L; l++)
            {
                B_re = B_im = 0;
                for (m = 0; m < M; m++)
                {
                    float32_t a_re, b_re;
                    a_re = A[l*SA + n + m*N + 0];
                    b_re = B[l*SB + m*P + p + 0];
                    B_re += (a_re*b_re);
                }
                Z[l*N*P + n*P + p] = B_re;
            }
#endif // 0

    int n, m, l, p;
    int SA = getSpace(M*N);
    int SB = getSpace(M*P);

    const xtfloat       * restrict pA;
    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pB0;
    const xb_vecN_2xf32 * restrict pB1;
    const xb_vecN_2xf32 * restrict pBl;
          xb_vecN_2xf32 * restrict pZ;

    valign vB, vB0, vZ;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, A0, A1, B0, B1, B2, B3;

    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
        {
            pBl = (const xb_vecN_2xf32 *)(B + l*SB);
            pZ = (xb_vecN_2xf32 *)(Z + l*N*P + n*P);
            vZ = BBE_ZALIGN();
            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pA = (const xtfloat *)(A + l*SA + n);
                pB = pBl;

                Acc = BBE_ZERON_2XF32();
                Acc1 = BBE_ZERON_2XF32();
                Acc2 = BBE_ZERON_2XF32();
                Acc3 = BBE_ZERON_2XF32();

                //__Pragma("loop_count min=2");
                for (m = 0; m < (M >> 1); m++)
                {
                    pB0 = pB;
                    vB = BBE_LAN_2XF32_PP(pB0);
                    BBE_LAVN_2XF32_XP(B0, vB, pB0, p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(B1, vB, pB0, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
                    BBE_LSN_2XF32_XP(A0, pA, N * sizeof(float32_t));
                    A0 = BBE_REPN_2XF32(A0, 0);

                    BBE_MULAN_2XF32(Acc, A0, B0);
                    BBE_MULAN_2XF32(Acc1, A0, B1);

                    pB1 = (const xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pB);


                    vB0 = BBE_LAN_2XF32_PP(pB1);
                    BBE_LAVN_2XF32_XP(B2, vB0, pB1, p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(B3, vB0, pB1, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
                    BBE_LSN_2XF32_XP(A1, pA, N * sizeof(float32_t));
                    A1 = BBE_REPN_2XF32(A1, 0);

                    BBE_MULAN_2XF32(Acc2, A1, B2);
                    BBE_MULAN_2XF32(Acc3, A1, B3);

                    pB = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pB);
                }
                Acc = BBE_ADDN_2XF32(Acc, Acc2);
                Acc1 = BBE_ADDN_2XF32(Acc1, Acc3);

                BBE_SAVN_2XF32_XP(Acc, vZ, pZ, p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc1, vZ, pZ, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));

                pBl = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH, (uintptr_t)pBl);
            }
            BBE_SAN_2XF32POS_FP(vZ, pZ);
        }
    }
}

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    MxNxP=MxNx1
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void rcomputeABmxnx1f(float32_t* Z,const float32_t* A,const float32_t* B, int M,int N,int L)
{
#if 0
    float32_t B_re;
    int n, m, l;
    int SA = getSpace(M*N);
    int SB = getSpace(M);

    for (n = 0; n < N; n++)
        for (l = 0; l < L; l++)
        {
            B_re = 0.f;
            for (m = 0; m < M; m++)
            {
                float32_t a_re, b_re;
                a_re = A[l*SA + n + m*N + 0];
                b_re = B[l*SB + m];
                B_re += (a_re*b_re);
            }
            Z[l*N + n] = B_re;
        }
#endif // 0

    int n, m, l;
    int SA = getSpace(M*N);
    int SB = getSpace(M);

    const xb_vecN_2xf32 * restrict pA;
    const xb_vecN_2xf32 * restrict pA0;
    const xb_vecN_2xf32 * restrict pAl;
    const xtfloat * restrict pB;
          xb_vecN_2xf32 * restrict pZ;

    valign vA, vZ;
    xb_vecN_2xf32 Acc, A0, B0;

    for (l = 0; l < L; l++)
    {
        pAl = (const xb_vecN_2xf32 *)(A + l*SA);
        pZ = (xb_vecN_2xf32 *)(Z + l*N);
        vZ = BBE_ZALIGN();
        for (n = N; n > 0; n -= VECLEN)
        {
            pA = pAl;
            pA0 = (const xb_vecN_2xf32 *)XT_ADDX4(N, (uintptr_t)pAl);
            pB = (const xtfloat *)(B + l*SB);

            Acc = BBE_ZERON_2XF32();

            __Pragma("loop_count min=2");
            for (m = 0; m < M >> 1; m++)
            {
                vA = BBE_LAN_2XF32_PP(pA);
                BBE_LAN_2XF32_IP(A0, vA, pA);
                pA = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N - 8, (uintptr_t)pA);
                BBE_LSN_2XF32_XP(B0, pB, sizeof(float32_t));
                B0 = BBE_REPN_2XF32(B0, 0);
                BBE_MULAN_2XF32(Acc, A0, B0);

                vA = BBE_LAN_2XF32_PP(pA0);
                BBE_LAN_2XF32_IP(A0, vA, pA0);
                pA0 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N - 8, (uintptr_t)pA0);
                BBE_LSN_2XF32_XP(B0, pB, sizeof(float32_t));
                B0 = BBE_REPN_2XF32(B0, 0);
                BBE_MULAN_2XF32(Acc, A0, B0);
            }
            BBE_SAVN_2XF32_XP(Acc, vZ, pZ, n * sizeof(float32_t));
            pAl += 1;
        }
        BBE_SAN_2XF32POS_FP(vZ, pZ);
    }
}

/*
    another algorithm for P!=1
*/
static void rfwdnxpf( 
          float32_t* restrict y, 
          const float32_t* restrict R, 
          const float32_t* restrict D, 
          const float32_t* restrict Z, 
          int N,int P,int L)
{
#if 0
    int SR = getSpace((N*(N + 1)) >> 1);
    int SY = getSpace(N*P);
    int SD = getSpace(N);
    int l, n, m, p;
    const float32_t *pR;
    float32_t B_re;
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 1 && N > 0);

    for (n = 0; n < N; n++)
    {
        for (p = 0; p < P; p++)
        {
            for (l = 0; l < L; l++)
            {
                pR = R + ((n*(n + 1)) >> 1) + l*SR;
                B_re = Z[l*N*P + n*P + p + 0];
                for (m = 0; m < n; m++)
                {
                    float32_t r_re;
                    float32_t y_re;
                    r_re = pR[m];
                    y_re = y[l*SY + m*P + p + 0];
                    B_re -= (y_re*r_re);
                }
                y[l*SY + m*P + p + 0] = B_re*D[l*SD + n + 0];
            }
        }
    }
#endif // 0

    int l, n, m, p;
    int SR = getSpace((N*(N + 1)) >> 1);
    int SY = getSpace(N*P);
    int SD = getSpace(N);

    const xtfloat       * restrict pD;
    const xtfloat       * restrict pR;
    const xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pY0;
    const xb_vecN_2xf32 * restrict pYl;
          xb_vecN_2xf32 * restrict pYw;
    const xb_vecN_2xf32 * restrict pZ;

    valign vY, vYw, vZ;
    xb_vecN_2xf32 Acc, Acc1, R0, Y0, Y1, D0;

    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 1 && N > 0);

    for (n = 0; n < N; n++)
    {
        pD = (const xtfloat *)(D + n);

        for (l = 0; l < L; l++)
        {
            pYw = (xb_vecN_2xf32 *)(y + l*SY + n*P);
            vYw = BBE_ZALIGN();
            pZ = (const xb_vecN_2xf32 *)(Z + l*N*P + n*P);
            vZ = BBE_LAN_2XF32_PP(pZ);
            pYl = (const xb_vecN_2xf32 *)(y + l*SY);

            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            D0 = BBE_REPN_2XF32(D0, 0);

            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pY0 = pYl;
                pR = (const xtfloat *)(R + l*SR + ((n*(n + 1)) >> 1));

                BBE_LAVN_2XF32_XP(Acc, vZ, pZ, p * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(Acc1, vZ, pZ, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));

                for (m = 0; m < n; m++)
                {
                    pY = pY0;
                    vY = BBE_LAN_2XF32_PP(pY);
                    BBE_LAVN_2XF32_XP(Y0, vY, pY, p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(Y1, vY, pY, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
                    BBE_LSN_2XF32_XP(R0, pR, sizeof(float32_t));
                    R0 = BBE_REPN_2XF32(R0, 0);

                    BBE_MULSN_2XF32(Acc, R0, Y0);
                    BBE_MULSN_2XF32(Acc1, R0, Y1);

                    pY0 = (const xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pY0);
                }
                Acc = BBE_MULN_2XF32(Acc, D0);
                Acc1 = BBE_MULN_2XF32(Acc1, D0);

                BBE_SAVN_2XF32_XP(Acc, vYw, pYw, p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc1, vYw, pYw, (p - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));

                pYl = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH, (uintptr_t)pYl);
            }

            BBE_SAN_2XF32POS_FP(vYw, pYw);
        }
    }
}
/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SA,SR,SD,SB,SY denote the number of data elements required to store
a matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next power
of two, otherwise it is matrix_size rounded up to the next multiple of the SIMD 
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SB = storage_size(M*P)
SY = storage_size(N*P)

Scratch size in bytes is defined by [r]cholfwd<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 M         Matrix dimension (number of rows in matrices A)
 N         Matrix dimension (number of columns and rows in 
           matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 A[L][SA]  Sequence of L complex matrices A
 B[L][SB]  Sequence of original right-side matrices B
 D[L][SD]  Reciprocal of main diagonal
Output:
 y[L][SY]   Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4 
4. M>=N
---------------------------------------------------------------------------*/
void rcholfwdmxnxpnf(
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int M, int N, int P, int L )
{
    float32_t* Z=(float32_t* )pScr;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    if (L<=0 || N<=0 || P<=0)  return;
    NASSERT(M % 4 == 0 && N % 4 == 0);

    {
        xb_vecN_2xf32 * pS = (xb_vecN_2xf32 *)(Z + L*N*P);
        valign vS;
        vS = BBE_ZALIGN();
        BBE_SAN_2XF32_IP(BBE_ZERON_2XF32(), vS, pS);
        BBE_SAN_2XF32POS_FP(vS, pS);
    }

    if (P==1)
    {
        int SR=getSpace((N*(N+1))>>1);
        int SY=getSpace(N);
        int SD=getSpace(N);
        rcomputeABmxnx1f(Z,A,B, M,N,L);
        rcholnfFwdrec(y,R,D,Z,N,L,SR,SD,SY,N);
    }
    else
    {
        // compute A'*B
        rcomputeABf(Z,A,B, M,N,P,L);
        rfwdnxpf(y,R,D,Z,N,P,L);
    }
}

size_t rcholfwdmxnxpnf_getScratchSize(int M,int N, int P,int L)
{
    size_t Z_size;
    NASSERT(M%4==0 && N%4==0);
    M=XT_MAX(0,M);
    N=XT_MAX(0,N);
    P=XT_MAX(0,P);
    L=XT_MAX(0,L);
    Z_size= (L*N*P*sizeof(float32_t) + 2 * BBE_SIMD_WIDTH);
    return (Z_size);
}

#else
DISCARD_FUN(void, rcholfwdmxnxpnf,(
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int M, int N, int P, int L ))

size_t rcholfwdmxnxpnf_getScratchSize(int M,int N, int P,int L)
{
  (void)M; (void)N; (void)P; (void)L;
  return 0;
}

#endif

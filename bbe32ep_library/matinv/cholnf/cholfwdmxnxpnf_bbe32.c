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
    Cholesky forward recursion, floating point complex data, 
    block format
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
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void computeABf(float32_t* Z,const float32_t* A,const float32_t* B, int M,int N,int P,int L)
{
#if 0
    float32_t B_re, B_im;
    int n, m, l, p;
    int SA = 2 * getSpace(M*N);
    int SB = 2 * getSpace(M*P);

    for (p = 0; p < P; p++)
        for (n = 0; n < N; n++)
            for (l = 0; l < L; l++)
            {
                B_re = B_im = 0;
                for (m = 0; m < M; m++)
                {
                    float32_t a_re, a_im, b_re, b_im;
                    a_re = A[l*SA + 2 * n + m*N * 2 + 0]; a_im = A[l*SA + 2 * n + m*N * 2 + 1];
                    b_re = B[l*SB + m*P * 2 + p * 2 + 0]; b_im = B[l*SB + m*P * 2 + p * 2 + 1];
                    B_re += (a_re*b_re) + (a_im*b_im);
                    B_im += (a_re*b_im) - (a_im*b_re);
                }
                Z[2 * l*N*P + 2 * n*P + 2 * p + 0] = B_re;
                Z[2 * l*N*P + 2 * n*P + 2 * p + 1] = B_im;
            }
#endif // 0

    int n, m, l, p;
    int SA = 2 * getSpace(M*N);
    int SB = 2 * getSpace(M*P);

    const long long     * restrict pA;
    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pB0;
    const xb_vecN_2xf32 * restrict pBl;
          xb_vecN_2xf32 * restrict pZ;

    valign vB, vZ;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, A0, B0, B1;
    xb_vecN_4x64 temp;

    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
        {
            pBl = (const xb_vecN_2xf32 *)(B + l*SB);
            pZ = (xb_vecN_2xf32 *)(Z + 2 * l*N*P + 2 * n*P);
            vZ = BBE_ZALIGN();
            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pA = (const long long *)(A + l*SA + 2 * n);
                pB0 = pBl;

                Acc = BBE_ZERON_2XF32();
                Acc1 = BBE_ZERON_2XF32();
                Acc2 = BBE_ZERON_2XF32();
                Acc3 = BBE_ZERON_2XF32();

                //__Pragma("loop_count min=2");
                for (m = 0; m < M; m++)
                {
                    pB = pB0;
                    vB = BBE_LAN_2XF32_PP(pB);
                    BBE_LAVN_2XF32_XP(B0, vB, pB, 2 * p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(B1, vB, pB, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
                    BBE_LSN_4X64_XP(temp, pA, 2 * N * sizeof(float32_t));
                    A0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
                    A0 = BBE_SHFLN_2XF32I(A0, BBE_SHFLI_REP_0X4);

                    BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
                    BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);

                    BBE_MULMASN_2XF32(Acc2, A0, B1, 0, 4);
                    BBE_MULMASN_2XF32(Acc3, A0, B1, 2, 11);

                    pB0 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pB0);
                }
                Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
                Acc = BBE_ADDN_2XF32(Acc, Acc1);

                BBE_SAVN_2XF32_XP(Acc, vZ, pZ, 2 * p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc2, vZ, pZ, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));

                pBl = (const xb_vecN_2xf32 *)XT_ADDX4(BBE_SIMD_WIDTH, (uintptr_t)pBl);
            }
            BBE_SAN_2XF32POS_FP(vZ, pZ);
        }
    }
}

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    P==1
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void computeABmxnx1f(float32_t* Z,const float32_t* A,const float32_t* B, int M,int N,int L)
{
#if 0
    float32_t B_re, B_im;
    int n, m, l;
    int SA = 2 * getSpace(M*N);
    int SB = 2 * getSpace(M);

    for (n = 0; n < N; n++)
        for (l = 0; l < L; l++)
        {
            B_re = B_im = 0;
            for (m = 0; m < M; m++)
            {
                float32_t a_re, a_im, b_re, b_im;
                a_re = A[l*SA + 2 * n + m*N * 2 + 0]; a_im = A[l*SA + 2 * n + m*N * 2 + 1];
                b_re = B[l*SB + m * 1 * 2 + 0]; b_im = B[l*SB + m * 1 * 2 + 1];
                B_re += (a_re*b_re) + (a_im*b_im);
                B_im += (a_re*b_im) - (a_im*b_re);
            }
            Z[2 * l*N + 2 * n + 0] = B_re;
            Z[2 * l*N + 2 * n + 1] = B_im;
        }
#endif // 0

    int n, m, l;
    int SA = 2 * getSpace(M*N);
    int SB = 2 * getSpace(M);

    const xb_vecN_2xf32 * restrict pA;
    const xb_vecN_2xf32 * restrict pAl;
    const long long     * restrict pB;
          xb_vecN_2xf32 * restrict pZ;

    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, A0, B0;
    xb_vecN_4x64 temp;

    for (l = 0; l < L; l++)
    {
        pAl = (const xb_vecN_2xf32 *)(A + l*SA);
        pZ = (xb_vecN_2xf32 *)(Z + 2 * l*N);
        for (n = 0; n < N; n += VECLEN)
        {
            pA = pAl;
            pB = (const long long *)(B + l*SB);
            Acc = BBE_ZERON_2XF32();
            Acc1 = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            Acc3 = BBE_ZERON_2XF32();

            __Pragma("loop_count min=2");
            for (m = 0; m < M >> 1; m++)
            {
                BBE_LVN_2XF32_XP(A0, pA, 2 * N * sizeof(float32_t));
                BBE_LSN_4X64_XP(temp, pB, 2 * sizeof(float32_t));
                B0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
                B0 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_0X4);
                BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);

                BBE_LVN_2XF32_XP(A0, pA, 2 * N * sizeof(float32_t));
                BBE_LSN_4X64_XP(temp, pB, 2 * sizeof(float32_t));
                B0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
                B0 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_0X4);
                BBE_MULMASN_2XF32(Acc2, A0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc3, A0, B0, 2, 11);
            }
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
            pAl += 1;
        }
    }
}
/*
    another algorithm for P!=1
*/
static void fwdnxpf( 
          float32_t* restrict y, 
          const float32_t* restrict R, 
          const float32_t* restrict D, 
          const float32_t* restrict Z, 
          int N,int P,int L)
{
#if 0
    int SR = 2 * getSpace((N*(N + 1)) >> 1);
    int SY = 2 * getSpace(N*P);
    int SD = 2 * getSpace(N);
    int l, n, m, p;
    const float32_t *pR;
    float32_t B_re, B_im;
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
                pR = R + (n*(n + 1)) + l*SR;
                B_re = Z[2 * l*N*P + 2 * n*P + 2 * p + 0];
                B_im = Z[2 * l*N*P + 2 * n*P + 2 * p + 1];
                for (m = 0; m < n; m++)
                {
                    float32_t r_re, r_im;
                    float32_t y_re, y_im;
                    r_re = pR[2 * m + 0];
                    r_im = pR[2 * m + 1];
                    y_re = y[l*SY + m*P * 2 + p * 2 + 0];
                    y_im = y[l*SY + m*P * 2 + p * 2 + 1];
                    B_re -= (y_re*r_re) + (y_im*r_im); // representation qA+qY
                    B_im -= (y_im*r_re) - (y_re*r_im);
                }
                y[l*SY + m*P * 2 + p * 2 + 0] = B_re*D[l*SD + 2 * n + 0];
                y[l*SY + m*P * 2 + p * 2 + 1] = B_im*D[l*SD + 2 * n + 1];
            }
        }
    }
#endif // 0

    int l, n, m, p;
    int SR = 2 * getSpace((N*(N + 1)) >> 1);
    int SY = 2 * getSpace(N*P);
    int SD = 2 * getSpace(N);

    const long long     * restrict pD;
    const long long     * restrict pR;
    const xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pY0;
    const xb_vecN_2xf32 * restrict pYl;
          xb_vecN_2xf32 * restrict pYw;
    const xb_vecN_2xf32 * restrict pZ;

    valign vY, vYw, vZ;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, R0, Y0, Y1, D0;
    xb_vecN_4x64 temp;

    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 1 && N > 0);

    for (n = 0; n < N; n++)
    {
        pD = (const long long *)(D + 2 * n);

        for (l = 0; l < L; l++)
        {
            pYw = (xb_vecN_2xf32 *)(y + l*SY + n*P * 2);
            vYw = BBE_ZALIGN();
            pZ = (const xb_vecN_2xf32 *)(Z + 2 * l*N*P + 2 * n*P);
            vZ = BBE_LAN_2XF32_PP(pZ);
            pYl = (const xb_vecN_2xf32 *)(y + l*SY);

            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            D0 = BBE_SHFLN_2XF32I(D0, BBE_SHFLI_REP_0X4);

            for (p = P; p > 0; p -= 2 * VECLEN)
            {
                pY0 = pYl;
                pR = (const long long *)(R + (n*(n + 1)) + l*SR);

                BBE_LAVN_2XF32_XP(Acc, vZ, pZ, 2 * p * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(Acc2, vZ, pZ, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
                Acc1 = BBE_ZERON_2XF32();
                Acc3 = BBE_ZERON_2XF32();

                for (m = 0; m < n; m++)
                {
                    pY = pY0;
                    vY = BBE_LAN_2XF32_PP(pY);
                    BBE_LAVN_2XF32_XP(Y0, vY, pY, 2 * p * sizeof(float32_t));
                    BBE_LAVN_2XF32_XP(Y1, vY, pY, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
                    BBE_LSN_4X64_XP(temp, pR, 2 * sizeof(float32_t));
                    R0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
                    R0 = BBE_SHFLN_2XF32I(R0, BBE_SHFLI_REP_0X4);

                    BBE_MULMASN_2XF32(Acc, R0, Y0, 3, 4);
                    BBE_MULMASN_2XF32(Acc1, R0, Y0, 1, 11);

                    BBE_MULMASN_2XF32(Acc2, R0, Y1, 3, 4);
                    BBE_MULMASN_2XF32(Acc3, R0, Y1, 1, 11);

                    pY0 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pY0);
                }
                Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
                Acc = BBE_ADDN_2XF32(Acc, Acc1);
                
                Acc = BBE_MULN_2XF32(Acc, D0);
                Acc2 = BBE_MULN_2XF32(Acc2, D0);
                BBE_SAVN_2XF32_XP(Acc, vYw, pYw, 2 * p * sizeof(float32_t));
                BBE_SAVN_2XF32_XP(Acc2, vYw, pYw, 2 * (p - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));

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
void cholfwdmxnxpnf(
            void * pScr,
            complex_float * restrict _y,
      const complex_float * restrict _R, 
      const complex_float * restrict _D,
      const complex_float * restrict _A, 
      const complex_float * restrict _B, 
      int M, int N, int P, int L )
{
          float32_t* restrict y=(      float32_t*)_y;
    const float32_t* restrict R=(const float32_t*)_R;
    const float32_t* restrict D=(const float32_t*)_D;
    const float32_t* restrict A=(const float32_t*)_A;
    const float32_t* restrict B=(const float32_t*)_B;
    float32_t* Z=(float32_t* )pScr;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    if (L<=0 || N<=0 || P<=0)  return;

    BBE_SVN_2XF32_I(BBE_ZERON_2XF32(), (xb_vecN_2xf32 *)(Z + 2 * L*N*P), 0);

    if (P==1)
    {
        int SR=2*getSpace((N*(N+1))>>1);
        int SY=2*getSpace(N);
        int SD=2*getSpace(N);
        // compute A'*B
        computeABmxnx1f(Z,A,B, M,N,L);
        cholnfFwdrec(y,R,D,Z,N,L,SR,SD,SY,2*N);
    }
    else
    {
        computeABf(Z,A,B, M,N,P,L);
        fwdnxpf(y,R,D,Z,N,P,L);
    }
}

size_t cholfwdmxnxpnf_getScratchSize(int M,int N, int P,int L)
{
    size_t Z_size;
    M=XT_MAX(0,M);
    N=XT_MAX(0,N);
    P=XT_MAX(0,P);
    L=XT_MAX(0,L);
    Z_size= (2*L*N*P*sizeof(float32_t) + 2 * BBE_SIMD_WIDTH);
    return (Z_size);
}

#else
DISCARD_FUN(void, cholfwdmxnxpnf,(
            void * pScr,
            complex_float * restrict _y,
      const complex_float * restrict _R, 
      const complex_float * restrict _D,
      const complex_float * restrict _A, 
      const complex_float * restrict _B, 
      int M, int N, int P, int L ))

size_t cholfwdmxnxpnf_getScratchSize(int M,int N, int P,int L)
{
  (void)M; (void)N; (void)P; (void)L;
  return 0;
}

#endif

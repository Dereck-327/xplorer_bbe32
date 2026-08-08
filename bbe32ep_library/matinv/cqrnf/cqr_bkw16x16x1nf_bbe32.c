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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, block format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#if HAVE_VFPU 

/*-------------------------------------------------------
    backward recursion: P==1

   Input:
    M, N, P      dimensional parameters
    L            number of matrices
   Input/output:
    x[L][SB][2]  at the input it is the sequence of L updated right parts Z=Q'B.
                 They will be replaced with MMSE solution vectors X (only N*P 
                 elements are used)
   Input:
    R[L][SA][2]  Upper triangle matrices R (only N*N 
                 elements of each matrix are used)
    D[L][SD][2]  reciprocal of main diagonal (mantissa, exponent) 
                 in the special format
-------------------------------------------------------*/
static void cqrnfBkw16x1(
                  float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            int L)
{
#if 0
    int m, k;
    int l;
    int SR = (16 * 16 << 1);
    int SX = (16 << 1);
    int SD = (16 << 1);

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);

    for (k = 16 - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t* pRt = R + l*SR + (k * 16 + (k + 1)) * 2;
            // calculate y(m,:)-R(m,:)*X, 1xP
            float32_t r_re, r_im;
            float32_t x_re, x_im;
            float32_t B_re, B_im;
            B_re = (x[l*SX + 2 * k + 0]);
            B_im = (x[l*SX + 2 * k + 1]);

            for (m = 0; m < 16 - k - 1; m++)
            {
                x_re = x[l*SX + 2 * (k + 1 + m) + 0];
                x_im = x[l*SX + 2 * (k + 1 + m) + 1];
                r_re = pRt[2 * m + 0];
                r_im = pRt[2 * m + 1];
                B_re -= (x_re*r_re) - (x_im*r_im);
                B_im -= (x_re*r_im) + (x_im*r_re);
            }
            x[l*SX + 2 * k + 0] = B_re*D[l*SD + 2 * k + 0];
            x[l*SX + 2 * k + 1] = B_im*D[l*SD + 2 * k + 1];
        }
    }
#endif // 0

    int k, M, N = 16;
    int l;
    int SR = (16 * 16 << 1);
    int SX = (16 << 1);
    int SD = (16 << 1);

    const long long     * restrict pY;
          long long     * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pR;
    const long long     * restrict pD;

    valign vR, vX;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, X0, R0, D0;
    xb_vecN_4x64 temp;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);

    k = N - 1;
    M = 0;    /* number of iterations in the innermost loop */
    for (; k >= XT_MAX(N - 4 - 1, 0); k--, M++)
    {
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(R + (k*N + (k + 1)) * 2);
        pY = (const long long *)(x + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
        __Pragma("loop_count min=1");
        for (l = 0; l < L; l++)
        {
            BBE_LSN_4X64_XP(temp, pY, SX * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * M * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * M * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - 2 * M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - 2 * M, (uintptr_t)pR);

            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pXw, SX * sizeof(float32_t));
        }
    }
    for (; k >= XT_MAX(N - 8 - 1, 0); k--, M++)
    {
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(R + (k*N + (k + 1)) * 2);
        pY = (const long long *)(x + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
        __Pragma("loop_count min=1");
        for (l = 0; l < L; l++)
        {
            BBE_LSN_4X64_XP(temp, pY, SX * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            Acc1 = BBE_MULMN_2XF32(X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc1, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - 2 * M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - 2 * M, (uintptr_t)pR);

            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pXw, SX * sizeof(float32_t));
        }
    }
    for (; k >= XT_MAX(N - 12 - 1, 0); k--, M++)
    {
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(R + (k*N + (k + 1)) * 2);
        pY = (const long long *)(x + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
        __Pragma("loop_count min=1");
        for (l = 0; l < L; l++)
        {
            BBE_LSN_4X64_XP(temp, pY, SX * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            Acc1 = BBE_MULMN_2XF32(X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            BBE_MULMASN_2XF32(Acc1, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - 2 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - 2 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc1, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - 2 * M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - 2 * M, (uintptr_t)pR);

            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pXw, SX * sizeof(float32_t));
        }
    }
    for (; k >= 0; k--, M++)
    {
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(R + (k*N + (k + 1)) * 2);
        pY = (const long long *)(x + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
        __Pragma("loop_count min=1");
        for (l = 0; l < L; l++)
        {
            BBE_LSN_4X64_XP(temp, pY, SX * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            Acc2 = BBE_MULMN_2XF32(X0, R0, 3, 4);
            Acc3 = BBE_MULMN_2XF32(X0, R0, 2, 11);
            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            Acc1 = BBE_MULMN_2XF32(X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            BBE_LAN_2XF32_IP(X0, vX, pX);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            BBE_MULMASN_2XF32(Acc2, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc3, X0, R0, 2, 11);
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - 3 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - 3 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc1, X0, R0, 3, 4);
            BBE_MULMASN_2XF32(Acc, X0, R0, 2, 11);
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - 2 * M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - 2 * M, (uintptr_t)pR);

            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pXw, SX * sizeof(float32_t));
        }
    }
}

/*-------------------------------------------------------------------------
Apply backward recursion process for QR decomposition for block ordered 
matrices.
Matrix sizes SA,SB are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SB=size(M*P)
SD=size(N)
Scratch size in bytes is defined by cqr_bkwmxnxpn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      Dimensional parameters
 L            Number of matrices
Input/output:
 X[L][SB]     On input it is the sequence of L updated right parts Z=Q'B.
              They will be replaced with MMSE solution vectors X (only N*P 
              elements are used)
Input:
 R[L][SA]     Upper triangular matrices R (only N*N 
              elements of each matrix are used)
 D[L*SD]      Reciprocals of main diagonal in a special format

Restrictions:
1. X, R, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void  cqr_bkw16x16x1nf(void *pScr,
                          complex_float* X,
                    const complex_float* R,
                    const complex_float* D,
                    int L)
{
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    (void)pScr;
    if (L<=0 ) return;
    cqrnfBkw16x1((float32_t*)X,(float32_t*)R,(float32_t*)D,L);
}

size_t cqr_bkw16x16x1nf_getScratchSize  (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==16 && M==16);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}

#else
DISCARD_FUN(void, cqr_bkw16x16x1nf, (void *pScr,
                    complex_float* X,
                    const complex_float* R,
                    const complex_float* D,
                    int L))
size_t cqr_bkw16x16x1nf_getScratchSize  (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==16 && M==16);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

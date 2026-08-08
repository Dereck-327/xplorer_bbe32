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
  QR decomposition, floating point, real data, block format
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
    M, N, P   dimensional parameters
    L         number of matrices
   Input/output:
    x[L][SB]  at the input it is the sequence of L updated right parts Z=Q'B.
              They will be replaced with MMSE solution vectors X (only N*P 
              elements are used)
   Input:
    R[L][SA]  Upper triangle matrices R (only N*N 
              elements of each matrix are used)
    D[L][SD]  reciprocal of main diagonal (mantissa, exponent) 
                 in the special format
-------------------------------------------------------*/
static void qrnfBkw8x1(
            float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            int L)
{
#if 0
    int m, k;
    int l;
    int SR = (8 * 8);
    int SX = (8);
    int SD = (8);

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);

    for (k = 8 - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t * pRt = R + l*SR + (k * 8 + (k + 1));
            // calculate y(m,:)-R(m,:)*X, 1xP
            float32_t r_re, x_re, B_re;
            B_re = (x[l*SX + k]);
            for (m = 0; m < 8 - k - 1; m++)
            {
                x_re = x[l*SX + (k + 1 + m)];
                r_re = pRt[m];
                B_re -= (x_re*r_re);
            }
            x[l*SX + k] = B_re*D[l*SD + k];
        }
    }
#endif // 0

    int k, M, N = 8;
    int l;
    int SR = (8 * 8);
    int SX = (8);
    int SD = (8);

    const xtfloat       * restrict pY;
          xtfloat       * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pR;
    const xtfloat       * restrict pD;

    valign vR, vX;
    xb_vecN_2xf32 Acc/*, Acc1, Acc2, Acc3*/, X0, R0, D0;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);

    k = N - 1;
    M = 0;    /* number of iterations in the innermost loop */
    for (; k >= 0; k--, M++)
    {
        pX = (const xb_vecN_2xf32 *)(x + (k + 1));
        pR = (const xb_vecN_2xf32 *)(R + (k*N + (k + 1)));
        pY = (const xtfloat *)(x + k);
        pD = (const xtfloat *)(D + k);
        pXw = (xtfloat *)(x + k);
        __Pragma("loop_count min=1");
        for (l = 0; l < L; l++)
        {
            BBE_LSN_2XF32_XP(Acc, pY, SX * sizeof(float32_t));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAVN_2XF32_XP(X0, vX, pX, M * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, M * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, X0, R0, 3, 12);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - M, (uintptr_t)pR);

            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_XP(Acc, pXw, SX * sizeof(float32_t));
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
void  qr_bkw8x8x1nf(void *pScr,
                    float32_t* X,
                    const float32_t* R,
                    const float32_t* D,
                    int L)
{
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    (void)pScr;
    if (L<=0 ) return;
    qrnfBkw8x1(X,R,D,L);
}

size_t qr_bkw8x8x1nf_getScratchSize    (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==8 && M==8);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}

#else
DISCARD_FUN(void, qr_bkw8x8x1nf, (void *pScr,
                    float32_t* X,
                    const float32_t* R,
                    const float32_t* D,
                    int L))

size_t qr_bkw8x8x1nf_getScratchSize    (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==8 && M==8);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

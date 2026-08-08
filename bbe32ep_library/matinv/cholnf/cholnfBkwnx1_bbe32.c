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
    BBE32 code for Cholesky backward recursion, block format complex 
    floating point
    IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cholnf_common.h"

#if HAVE_VFPU

#define VECLEN (BBE_SIMD_WIDTH/4)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-2));
    // round up to the  next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*--------------------------------------------------
    backward recursion: P==1
    Input:
    Rt
    D
    y
    N,L  size of matrices
    Output:
    X
--------------------------------------------------*/
void cholnfBkwnx1( float32_t* restrict x, 
            const float32_t* restrict Rt,
            const float32_t* restrict D,
            const float32_t* restrict y, 
                  int N,int L)
{
    int m,k,M;
    int l;
    int SX=2*getSpace(N);
    int SD=2*getSpace(N);
    //int SR=2*getSpace((N*(N+1))>>1);
    int SR=((N*(N-1)));
#if 0
    float32_t r_re, r_im;
    float32_t x_re, x_im;
    float32_t B_re, B_im;
#endif

    const long long     * restrict pY;
          long long     * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pR;
    const long long     * restrict pD;

    valign vR, vX;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, X0, R0, D0;
    xb_vecN_4x64 temp;

    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    k = N - 1;
    M = 0;    /* number of iterations in the innermost loop */
    for (; k >= XT_MAX(N - 4 - 1, 0); k--, M++)
    {
#if 0
        for (l = 0; l<L; l++)
        {
            const float32_t* pRt = Rt + l*SR + (N - k - 1)*(N - k - 2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re = (y[l*SX + 2 * k + 0]);
            B_im = (y[l*SX + 2 * k + 1]);
            NASSERT(M <= 4);
            for (m = 0; m<M; m++)
            {
                x_re = x[l*SX + 2 * (k + 1 + m) + 0];
                x_im = x[l*SX + 2 * (k + 1 + m) + 1];
                r_re = pRt[2 * m + 0];
                r_im = pRt[2 * m + 1];
                B_re -= (x_re* r_re) - (x_im* r_im);
                B_im -= (x_re* r_im) + (x_im* r_re);
            }
            x[l*SX + 2 * k + 0] = B_re*D[l*SD + 2 * k + 0];
            x[l*SX + 2 * k + 1] = B_im*D[l*SD + 2 * k + 1];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (N - k - 1)*(N - k - 2));
        pY = (const long long *)(y + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
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
#endif
    }
    for (; k >= XT_MAX(N - 8 - 1, 0); k--, M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); 
            B_im=(y[l*SX+2*k+1]); 
            NASSERT(M>4 && M<=8);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=(x_re* r_re)-(x_im* r_im);
                B_im-=(x_re* r_im)+(x_im* r_re);
            }
            x[l*SX+2*k+0]=B_re*D[l*SD+2*k+0];
            x[l*SX+2*k+1]=B_im*D[l*SD+2*k+1];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (N - k - 1)*(N - k - 2));
        pY = (const long long *)(y + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
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
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - VECLEN) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - VECLEN) * sizeof(float32_t));
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
#endif
    }
    for (; k >= XT_MAX(N - 12 - 1, 0); k--, M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); 
            B_im=(y[l*SX+2*k+1]); 
            NASSERT(M>8 && M<=12);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=(x_re* r_re)-(x_im* r_im);
                B_im-=(x_re* r_im)+(x_im* r_re);
            }
            x[l*SX+2*k+0]=B_re*D[l*SD+2*k+0];
            x[l*SX+2*k+1]=B_im*D[l*SD+2*k+1];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (N - k - 1)*(N - k - 2));
        pY = (const long long *)(y + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
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
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - 2 * VECLEN) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - 2 * VECLEN) * sizeof(float32_t));
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
#endif
    }
    for (; k >= XT_MAX(N - 16 - 1, 0); k--, M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); 
            B_im=(y[l*SX+2*k+1]); 
            NASSERT(M>12 && M<=16);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=(x_re* r_re)-(x_im* r_im);
                B_im-=(x_re* r_im)+(x_im* r_re);
            }
            x[l*SX+2*k+0]=B_re*D[l*SD+2*k+0];
            x[l*SX+2*k+1]=B_im*D[l*SD+2*k+1];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (N - k - 1)*(N - k - 2));
        pY = (const long long *)(y + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
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
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M - 3 * VECLEN) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M - 3 * VECLEN) * sizeof(float32_t));
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
#endif
    }

    for (; k >= 0; k--, M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); 
            B_im=(y[l*SX+2*k+1]); 
            NASSERT(M>16 );
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=(x_re* r_re)-(x_im* r_im);
                B_im-=(x_re* r_im)+(x_im* r_re);
            }
            x[l*SX+2*k+0]=B_re*D[l*SD+2*k+0];
            x[l*SX+2*k+1]=B_im*D[l*SD+2*k+1];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + 2 * (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (N - k - 1)*(N - k - 2));
        pY = (const long long *)(y + 2 * k);
        pD = (const long long *)(D + 2 * k);
        pXw = (long long *)(x + 2 * k);
        for (l = 0; l < L; l++)
        {
            BBE_LSN_4X64_XP(temp, pY, SX * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            
            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * (M & (2 * VECLEN - 1)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * (M & (2 * VECLEN - 1)) * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
            Acc1 = BBE_MULMN_2XF32(X0, R0, 2, 11);
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * ((M & (2 * VECLEN - 1)) - VECLEN) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, 2 * ((M & (2 * VECLEN - 1)) - VECLEN) * sizeof(float32_t));
            Acc2 = BBE_MULMN_2XF32(X0, R0, 3, 4);
            Acc3 = BBE_MULMN_2XF32(X0, R0, 2, 11);

            for (m = 0; m < M / (2 * VECLEN); m++)
            {
                BBE_LAN_2XF32_IP(X0, vX, pX);
                BBE_LAN_2XF32_IP(R0, vR, pR);
                BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
                BBE_MULMASN_2XF32(Acc1, X0, R0, 2, 11);

                BBE_LAN_2XF32_IP(X0, vX, pX);
                BBE_LAN_2XF32_IP(R0, vR, pR);
                BBE_MULMASN_2XF32(Acc2, X0, R0, 3, 4);
                BBE_MULMASN_2XF32(Acc3, X0, R0, 2, 11);
            }
            
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
#endif
    }
}
#endif

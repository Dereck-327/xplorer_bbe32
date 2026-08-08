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
    BBE32 code for Cholesky backward recursion, block format real 
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
#include "rcholnf_common.h"

#if HAVE_VFPU

#define VECLEN (BBE_SIMD_WIDTH/2)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the  next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/* --------------------------------------------------
    backward recursion: P==1
    Input:
    Rt
    D
    y
    Output:
    x
--------------------------------------------------*/
void rcholnfBkwnx1(float32_t* restrict x, 
            const float32_t* restrict Rt,
            const float32_t* restrict D,
            const float32_t* restrict y, 
                  int N,int L)
{
    int m,k,M;
    int l;
    int SX=getSpace(N);
    int SD=getSpace(N);
    //int SR=2*getSpace((N*(N+1))>>1);
    int SR=((N*(N-1))>>1);
#if 0
    float32_t r_re;
    float32_t x_re;
    float32_t B_re;
#endif

    const xtfloat       * restrict pY;
          xtfloat       * restrict pXw;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pR;
    const xtfloat       * restrict pD;

    valign vR, vX;
    xb_vecN_2xf32 Acc/*, Acc1, Acc2, Acc3*/, X0, R0, D0;

    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    k = N - 1;
    M = 0;    /* number of iterations in the innermost loop */
    for (; k>=XT_MAX(N-8-1,0); k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(((N-k-1)*(N-k-2))>>1);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+k+0]); 
            NASSERT(M<=8);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+(k+1+m)+0];
                r_re=pRt[m+0];
                B_re-=(x_re* r_re);
            }
            x[l*SX+k+0]=B_re*D[l*SD+k+0];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (((N - k - 1)*(N - k - 2)) >> 1));
        pY = (const xtfloat *)(y + k);
        pD = (const xtfloat *)(D + k);
        pXw = (xtfloat *)(x + k);
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

            //Acc = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(Acc));
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_XP(Acc, pXw, SX * sizeof(float32_t));
        }
#endif
    }
    for (; k>=XT_MAX(N-16-1,0); k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(((N-k-1)*(N-k-2))>>1);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+k+0]); 
            NASSERT(M>8 && M<=16);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+(k+1+m)+0];
                r_re=pRt[m+0];
                B_re-=(x_re* r_re);
            }
            x[l*SX+k+0]=B_re*D[l*SD+k+0];
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (((N - k - 1)*(N - k - 2)) >> 1));
        pY = (const xtfloat *)(y + k);
        pD = (const xtfloat *)(D + k);
        pXw = (xtfloat *)(x + k);
        for (l = 0; l < L; l++)
        {
            BBE_LSN_2XF32_XP(Acc, pY, SX * sizeof(float32_t));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            BBE_LAVN_2XF32_XP(X0, vX, pX, M * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, M * sizeof(float32_t));
            BBE_MULSN_2XF32(Acc, X0, R0);
            BBE_LAVN_2XF32_XP(X0, vX, pX, (M - VECLEN) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(R0, vR, pR, (M - VECLEN) * sizeof(float32_t));
            BBE_MULSN_2XF32(Acc, X0, R0);

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - M, (uintptr_t)pR);

            Acc = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(Acc));
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2), Acc);
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_XP(Acc, pXw, SX * sizeof(float32_t));
        }
#endif
    }

    for (; k>=0; k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
        for(l=0; l<L; l++)
        {
            const float32_t* pRt=Rt+l*SR+(((N-k-1)*(N-k-2))>>1);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+k+0]); 
            NASSERT(M>16);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+(k+1+m)+0];
                r_re=pRt[m+0];
                B_re-=(x_re* r_re);
            }
            x[l*SX+k+0]=B_re*D[l*SD+k+0];
        }
        }
#else
        pX = (const xb_vecN_2xf32 *)(x + (k + 1));
        pR = (const xb_vecN_2xf32 *)(Rt + (((N - k - 1)*(N - k - 2)) >> 1));
        pY = (const xtfloat *)(y + k);
        pD = (const xtfloat *)(D + k);
        pXw = (xtfloat *)(x + k);
        for (l = 0; l < L; l++)
        {
            BBE_LSN_2XF32_XP(Acc, pY, SX * sizeof(float32_t));

            vX = BBE_LAN_2XF32_PP(pX);
            vR = BBE_LAN_2XF32_PP(pR);

            for (m = M; m > 0; m -= VECLEN)
            {
                BBE_LAVN_2XF32_XP(X0, vX, pX, m * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(R0, vR, pR, m * sizeof(float32_t));
                BBE_MULSN_2XF32(Acc, X0, R0);
            }

            pX = (const xb_vecN_2xf32 *)XT_ADDX4(SX - M, (uintptr_t)pX);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - M, (uintptr_t)pR);

            Acc = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(Acc));
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2), Acc);
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            //Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);

            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_XP(Acc, pXw, SX * sizeof(float32_t));
        }
#endif
    }
}
#endif

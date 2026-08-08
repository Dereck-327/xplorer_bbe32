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
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "qrnf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    xt[L*M]   transposed x
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)
    temporary:
    pScr[]    scratch, defined by qrnfHousholder_getScratchSz()
-------------------------------------------------------*/
void qrnfHousholder(void* pScr,
                           float32_t* restrict v,
                           float32_t* restrict Fi,
                           float32_t *restrict D,
                           const float32_t* restrict x, 
                           const float32_t* restrict xt, 
                           int M, int SD, int L)
{
#if 0
    float32_t *tmp = (float32_t *)pScr;   //[L]
    float32_t *invd = tmp + L;             //[L]
    float32_t *norm = invd + L;             //[L]
    int l, m;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++) tmp[l] = 0.f;
    for (l = 0; l < L; l++)
        for (m = 0; m < M; m++)
        {
            float32_t  x_re;
            x_re = xt[m*L + l];
            tmp[l] += (x_re*x_re);
        }
    for (l = 0; l < L; l++)
    {
        tmp[l] = 1.f / sqrtf(tmp[l]);
        invd[l] = 1.f / tmp[l];
    }
    for (l = 0; l < L; l++)
    {
        D[l*SD] = tmp[l];
    }
    // compute fi=x0/sqrt(x0'x0)
    for (l = 0; l<L; l++)
    {
        Fi[l] = xt[l]>0.f ? -1.f : 1.f;
    }
    for (l = 0; l < L; l++) norm[l] = 0.f;
    for (l = 0; l < L; l++)
        for (m = 1; m < M; m++)
        {
            float32_t x_re;
            x_re = x[l*M + m];  // or xt[m*L+l]
            norm[l] += (x_re*x_re);
        }

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    for (l = 0; l < L; l++)
    {
        float32_t  f_re, x_re, v_re;
        float32_t d0, invd0;
        float32_t  A_re;
        f_re = Fi[l];
        invd0 = invd[l];
        x_re = xt[l];
        v_re = x_re - f_re*invd0;
        A_re = norm[l];
        A_re += v_re*v_re;
        d0 = 1.0f / sqrtf(A_re);
        invd[l] = d0;
        tmp[l] = v_re*d0;
    }
    // normalization
    for (l = 0; l < L; l++)
    {
        v[M*l] = tmp[l];
    }
    for (l = 0; l < L; l++)
    {
        float32_t d0, x_re;
        d0 = invd[l];
        for (m = 1; m < M; m++)
        {
            x_re = x[l*M + m];
            v[M*l + m] = x_re*d0;
        }
    }
#endif // 0

    int l, m;
    int L_ = ((L + (BBE_SIMD_WIDTH / 2 - 1))&(~(BBE_SIMD_WIDTH / 2 - 1)));

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pZr;
    const xb_vecN_2xf32 * restrict pZr1;
          xb_vecN_2xf32 * restrict pZw;
          xb_vecN_2xf32 * restrict pZw1;
          xtfloat       * restrict pD;
    const xtfloat       * restrict pZr_;
          xb_vecN_2xf32 * restrict pF;
          xtfloat       * restrict pV;
          xb_vecN_2xf32 * restrict pV1;
          xtfloat       * restrict pZw_;

    xb_vecN_2xf32 Acc, Acc0, X0, X1, F0, rsqrt;
    valign vX, vZr, vF, vV;
    vboolN_2 tmp;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    // compute d=1/sqrt(x'x)
    pX = (const xb_vecN_2xf32 *)x;
    pZw_ = (xtfloat *)pScr;
    for (l = 0; l < L; l++)
    {
        Acc = BBE_ZERON_2XF32();

        BBE_LAVN_2XF32_XP(X0, vX, pX, sizeof(float32_t));
        Acc0 = BBE_MULN_2XF32(X0, X0);
        for (m = M - 1; m > 0; m -= (BBE_SIMD_WIDTH / 2))
        {
            BBE_LAVN_2XF32_XP(X0, vX, pX, m*sizeof(float32_t));
            BBE_MULAN_2XF32(Acc, X0, X0);
        }
        Acc = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(Acc));
        Acc0 = BBE_ADDN_2XF32(Acc, Acc0);
        BBE_SSN_2XF32_X(Acc, pZw_, 2 * L_*sizeof(float32_t));
        BBE_SSN_2XF32_IP(Acc0, pZw_, sizeof(float32_t));
    }

    pZr = (const xb_vecN_2xf32 *)pScr;
    pZw = (xb_vecN_2xf32 *)((float32_t *)pScr + L_);
    pZw1 = (xb_vecN_2xf32 *)pScr;
    for (l = L; l > 0; l -= (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(Acc, pZr, 0);
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_LVN_2XF32_IP(Acc, pZr, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(rsqrt, pZw1, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Acc, pZw, 2 * BBE_SIMD_WIDTH);
    }

    pZr = (const xb_vecN_2xf32 *)pScr;
    pD = (xtfloat *)D;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(Acc, vZr, pZr, sizeof(float32_t));
        BBE_SSN_2XF32_XP(Acc, pD, SD*sizeof(float32_t));
    }

    // compute fi=x0/sqrt(x0'x0)
    pX = (const xb_vecN_2xf32 *)xt;
    pF = (xb_vecN_2xf32 *)Fi;
    vF = BBE_ZALIGN();
    for (l = 0; l < L / (BBE_SIMD_WIDTH / 2); l++)
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SAN_2XF32_IP(X0, vF, pF);
    }
    l = L % (BBE_SIMD_WIDTH / 2);
    if (l)
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SAVN_2XF32_XP(X0, vF, pF, l*sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vF, pF);

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)xt;
    pZr = (const xb_vecN_2xf32 *)((float32_t *)pScr + L_);
    pZw = (xb_vecN_2xf32 *)((float32_t *)pScr + L_);
    pZr1 = (const xb_vecN_2xf32 *)((float32_t *)pScr + 2 * L_);
    pZw1 = (xb_vecN_2xf32 *)pScr;
    pF = (xb_vecN_2xf32 *)Fi;
    vF = BBE_LAN_2XF32_PP(pF);
    for (l = L; l > 0; l -= (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(Acc, pZr1, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(X1, pZr, 2 * BBE_SIMD_WIDTH);
        BBE_LAN_2XF32_IP(F0, vF, pF);

        BBE_MULSN_2XF32(X0, X1, F0);
        BBE_MULAN_2XF32(Acc, X0, X0);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SVN_2XF32_IP(rsqrt, pZw, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pZw1, 2 * BBE_SIMD_WIDTH);
    }

    // normalization
    pZr = (const xb_vecN_2xf32 *)pScr;
    pV = (xtfloat *)v;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(Acc, vZr, pZr, sizeof(float32_t));
        BBE_SSN_2XF32_XP(Acc, pV, M*sizeof(float32_t));
    }

    //pX = (const xb_vecN_2xf32 *)x;
    pV1 = (xb_vecN_2xf32 *)v;
    pZr_ = (const xtfloat *)((float32_t *)pScr + L_);
    for (l = 0; l < L; l++)
    {
        pX = (const xb_vecN_2xf32 *)(x + M*l + 1);
        pV1 = (xb_vecN_2xf32 *)XT_ADDI((uintptr_t)pV1, 4);
        vX = BBE_LAN_2XF32_PP(pX);
        vV = BBE_ZALIGN();

        BBE_LSN_2XF32_IP(rsqrt, pZr_, sizeof(float32_t));
        rsqrt = BBE_REPN_2XF32(rsqrt, 0);

        for (m = M - 1; m > 0; m -= (BBE_SIMD_WIDTH / 2))
        {
            BBE_LAN_2XF32_IP(X0, vX, pX);
            Acc = BBE_MULN_2XF32(X0, rsqrt);
            BBE_SAVN_2XF32_XP(Acc, vV, pV1, m*sizeof(float32_t));
        }
        BBE_SAN_2XF32POS_FP(vV, pV1);
    }
}
#endif

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
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "cqrsf_common.h"
#include <float.h>
#if HAVE_VFPU

/*-----------------------------------------------------------
    find Householder vectors (V and Fi), diagonal element D
    Input:
    x[]       pointer to the diagonal element of original matrix A given 
              in streaming order
    NL,SV,SD  strides for x, V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Householder vectors (M elements each)
-----------------------------------------------------------*/
void cqrfsHousholder(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L)
{
#if 0
    float32_t *invd = (float32_t *)pScr; // [2L]
    float32_t *tmp = invd + 2 * L; // [2L]
    int l, m;
    int _NL = 2 * N*L;
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);

    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++)
    {
        float32_t  x_re, x_im;
        float32_t A_re;
        A_re = FLT_MIN;
        for (m = 0; m < M; m++)
        {
            x_re = x[l * 2 + m*_NL + 0];
            x_im = x[l * 2 + m*_NL + 1];
            A_re += (x_re*x_re) + (x_im*x_im);
        }
        D[l * 2 + 0] = D[l * 2 + 1] = 1.f / sqrtf(A_re);
        invd[l * 2 + 0] = A_re*D[l * 2 + 0];
        invd[l * 2 + 1] = A_re*D[l * 2 + 1];
    }
    // compute fi=x0/sqrt(x0'x0)
    for (l = 0; l < L; l++)
    {
        float32_t d0;
        float32_t x_re, x_im;
        float32_t A_re;
        x_re = x[l * 2 + 0];
        x_im = x[l * 2 + 1];
        A_re = (x_re*x_re) + (x_im*x_im);
        d0 = 1.f / sqrtf(A_re);
        if (fabsf(A_re) > FLT_MIN)
        {
            Fi[2 * l + 0] = -x_re*d0;
            Fi[2 * l + 1] = -x_im*d0;
        }
        else
        {
            Fi[2 * l + 0] = 1.f;
            Fi[2 * l + 1] = 0.f;
        }
    }
    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    for (l = 0; l < L; l++)
    {
        float32_t A_re, x_re, x_im;
        A_re = FLT_MIN;
        for (m = 1; m < M; m++)
        {
            x_re = x[l * 2 + m*_NL + 0];
            x_im = x[l * 2 + m*_NL + 1];
            A_re += (x_re*x_re) + (x_im*x_im);
        }
        tmp[2 * l + 0] = tmp[2 * l + 1] = A_re;
    }

    for (l = 0; l < L; l++)
    {
        float32_t f_re, f_im, x_re, x_im, v_re, v_im;
        float32_t d0, invd0;
        float32_t A_re;
        f_re = Fi[l * 2 + 0];
        f_im = Fi[l * 2 + 1];
        invd0 = invd[2 * l + 0];
        x_re = x[l * 2 + 0];
        x_im = x[l * 2 + 1];
        v_re = x_re - (f_re*invd0);
        v_im = x_im - (f_im*invd0);
        A_re = tmp[2 * l + 0];
        A_re += (v_re*v_re) + (v_im*v_im);
        d0 = invd[2 * l + 0] = invd[2 * l + 1] = 1.f / sqrtf(A_re);
        v[l * 2 + 0] = v_re*d0;
        v[l * 2 + 1] = v_im*d0;
    }
    for (l = 0; l < L; l++)
    {
        float32_t x_re, x_im, d0;
        d0 = invd[2 * l + 0];
        for (m = 1; m < M; m++)
        {
            x_re = x[l * 2 + m*_NL + 0];
            x_im = x[l * 2 + m*_NL + 1];
            v[l * 2 + m * 2 * L + 0] = x_re*d0;
            v[l * 2 + m * 2 * L + 1] = x_im*d0;
        }
    }
#endif // 0

    int l, m;

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX0;
    const xb_vecN_2xf32 * restrict pZr;
          xb_vecN_2xf32 * restrict pZ;
          xb_vecN_2xf32 * restrict pD;
          xb_vecN_2xf32 * restrict pF;
          xb_vecN_2xf32 * restrict pV;
          xb_vecN_2xf32 * restrict pV0;

    xb_vecN_2xf32 Acc, Acc0, Acc1, Acc2, Acc3, X0, X1, F0, rsqrt;
    xb_vecN_2xf32 r0;
    const float32_t ALIGN(32) rep_tbl[8] = { 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f };

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);

    // compute d=1/sqrt(x'x)
    pX0 = (const xb_vecN_2xf32 *)x;
    pZ = (xb_vecN_2xf32 *)pScr;
    pD = (xb_vecN_2xf32 *)D;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        pX = pX0;

        Acc0 = FLT_MIN;
        Acc1 = BBE_ZERON_2XF32();
        Acc2 = BBE_ZERON_2XF32();
        Acc3 = BBE_ZERON_2XF32();

        BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        for (m = 2; m < M; m += 2)
        {
            BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
            X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
            BBE_MULAN_2XF32(Acc0, X0, X0);
            BBE_MULAN_2XF32(Acc1, X1, X1);

            BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
            X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
            BBE_MULAN_2XF32(Acc2, X0, X0);
            BBE_MULAN_2XF32(Acc3, X1, X1);
        }
        if ((M - 1) & 1)
        {
            BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
            X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
            BBE_MULAN_2XF32(Acc0, X0, X0);
            BBE_MULAN_2XF32(Acc1, X1, X1);
        }
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc1 = BBE_ADDN_2XF32(Acc0, Acc2);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
#if 0
        Acc = FLT_MIN;
        Acc1 = BBE_ZERON_2XF32();

        BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        BBE_MULAN_2XF32(Acc, X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        for (m = 1; m < M; m++)
        {
            BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc1, X0, X0);
        }
        Acc2 = BBE_SHFLN_2XF32I(Acc1, BBE_SHFLI_SWAP_2);
        Acc1 = BBE_ADDN_2XF32(Acc1, Acc2);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc1 = BBE_ADDN_2XF32(Acc1, FLT_MIN);
#endif

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X(Acc1, pZ, 2 * L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
        pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pX0);
    }

    // compute fi=x0/sqrt(x0'x0)
    pX = (const xb_vecN_2xf32 *)x;
    pF = (xb_vecN_2xf32 *)Fi;
    r0 = BBE_LVN_2XF32_I((const xb_vecN_2xf32 *)rep_tbl, 0); /////////////
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        vboolN_2 tmp;

        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        Acc = BBE_NEGN_2XF32(Acc);
        Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(r0), BBE_MOVNX16_FROMN_2XF32(Acc), BBE_MOVN_FROMN_2 (tmp))); /////////////
        BBE_SVN_2XF32_IP(Acc, pF, 2 * BBE_SIMD_WIDTH);
    }

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)x;
    pZr = (const xb_vecN_2xf32 *)pScr;
    pZ = (xb_vecN_2xf32 *)pScr;
    pF = (xb_vecN_2xf32 *)Fi;
    pV = (xb_vecN_2xf32 *)v;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_LVN_2XF32_X(pZr, 2 * L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(X1, pZr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        BBE_MULSN_2XF32(X0, X1, F0);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        BBE_MULAN_2XF32(Acc, X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SVN_2XF32_IP(rsqrt, pZ, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);
    }

    pX0 = (const xb_vecN_2xf32 *)(x + 2 * N*L);
    pZr = (const xb_vecN_2xf32 *)pScr;
    pV0 = (xb_vecN_2xf32 *)(v + 2 * L);
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        pX = pX0;
        pV = pV0;
        BBE_LVN_2XF32_IP(rsqrt, pZr, 2 * BBE_SIMD_WIDTH);
        for (m = 1; m < M; m++)
        {
            BBE_LVN_2XF32_XP(X0, pX, 2 * N*L * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(X0, rsqrt);
            BBE_SVN_2XF32_XP(Acc, pV, 2 * L * sizeof(float32_t));
        }
        pX0 += 1;
        pV0 += 1;
    }
}
#endif

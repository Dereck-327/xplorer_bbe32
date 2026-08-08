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
  QR decomposition, floating point, real data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "qrsf_common.h"

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
void qrfsHousholder(void* pScr,
                                 float32_t* restrict v,
                                 float32_t* restrict Fi,
                                 float32_t *restrict D,
                           const float32_t* restrict x, 
                           int SV, int M, int N, int L)
{
#if 0
    int l, m;
    int _NL = N*L;
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++)
    {
        float32_t  x_re;
        float32_t A_re;
        A_re = 0;
        for (m = 0; m < M; m++)
        {
            x_re = x[l + m*_NL];
            A_re += (x_re*x_re);
        }
        D[l] = 1.f / sqrtf(A_re);
    }
    for (l = 0; l<L; l++)
    {
        Fi[l] = x[l]>0.f ? -1.f : 1.f;
    }
    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    for (l = 0; l < L; l++)
    {
        float32_t f_re, x_re, v_re;
        float32_t d0, invd0;
        float32_t A_re;
        f_re = Fi[l];
        d0 = D[l];
        invd0 = 1.f / d0;
        x_re = x[l];
        v_re = x_re - (f_re*invd0);

        A_re = 0.f;
        for (m = 1; m < M; m++)
        {
            x_re = x[l + m*_NL];
            A_re += (x_re*x_re);
        }
        A_re += (v_re*v_re);
        d0 = 1.f / sqrtf(A_re);
        // normalization
        v[l] = v_re*d0;
        for (m = 1; m < M; m++)
        {
            x_re = x[l + m*_NL];
            v[l + m*L] = x_re*d0;
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

    xb_vecN_2xf32 Acc, Acc1, /*Acc2, Acc3, */X0, X1, rsqrt;
    vboolN_2 tmp;

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    // compute d=1/sqrt(x'x)
    pX0 = (const xb_vecN_2xf32 *)x;
    pZ = (xb_vecN_2xf32 *)pScr;
    pD = (xb_vecN_2xf32 *)D;
    pF = (xb_vecN_2xf32 *)Fi;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        pX = pX0;
        Acc1 = BBE_ZERON_2XF32();

        BBE_LVN_2XF32_XP(X0, pX, N*L * sizeof(float32_t));
        Acc = BBE_MULN_2XF32(X0, X0);

        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SVN_2XF32_IP(X0, pF, 2 * BBE_SIMD_WIDTH);

        for (m = 1; m < M; m++)
        {
            BBE_LVN_2XF32_XP(X0, pX, N*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc1, X0, X0);
        }
        Acc = BBE_ADDN_2XF32(Acc, Acc1);


        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_X(Acc1, pZ, L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
        pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pX0);
    }

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)x;
    pZr = (const xb_vecN_2xf32 *)pScr;
    pZ = (xb_vecN_2xf32 *)pScr;
    pV = (xb_vecN_2xf32 *)v;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        xb_vecN_2xf32 Xp, Xm;

        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_LVN_2XF32_X(pZr, L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(X1, pZr, 2 * BBE_SIMD_WIDTH);

        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        BBE_ADDSUBN_2XF32(Xp, Xm, X0, X1);
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(Xm), BBE_MOVNX16_FROMN_2XF32(Xp), BBE_MOVN_FROMN_2(tmp)));
        BBE_MULAN_2XF32(Acc, X0, X0);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SVN_2XF32_IP(rsqrt, pZ, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);
    }

    pX0 = (const xb_vecN_2xf32 *)(x + N*L);
    pZr = (const xb_vecN_2xf32 *)pScr;
    pV0 = (xb_vecN_2xf32 *)(v + L);
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        pX = pX0;
        pV = pV0;
        BBE_LVN_2XF32_IP(rsqrt, pZr, 2 * BBE_SIMD_WIDTH);
        for (m = 1; m < M; m++)
        {
            BBE_LVN_2XF32_XP(X0, pX, N*L * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(X0, rsqrt);
            BBE_SVN_2XF32_XP(Acc, pV, L * sizeof(float32_t));
        }
        pX0 += 1;
        pV0 += 1;
    }
}
#endif

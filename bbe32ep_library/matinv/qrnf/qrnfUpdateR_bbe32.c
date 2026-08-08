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
    partial update of R matrix
    Fi[L][SV] diagonal rotation matrix (only 0-th element filled)
    v[L][SV]  Housholder vector (M elements filled)
    SV,SD     strides for V/Fi and D
    M,N       matrix size
    L         number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    Temporary:
    Z[N*L]
-------------------------------------------------------*/
void qrnfUpdateR(float32_t* restrict Z,
                float32_t* restrict R,
                const float32_t* restrict v,
                int SA, int M,int N, int N0, int L)
{
#if 0
    int m, l, n;
    float32_t A_re;
    // compute v'B first
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
        {
            A_re = 0.f;
            for (m = 0; m < M; m++)
            {
                float32_t v_re, b_re;
                v_re = v[M*l + m];
                b_re = R[m*N0 + n];
                A_re += 2.f*v_re*b_re;
            }
            Z[n] = A_re;
        }
        R += SA;
        Z += N;
    }
    R -= SA*L; Z -= N*L;
    for (l = 0; l < L; l++)
    {
        for (m = 0; m < M; m++)
        {
            for (n = 0; n < N; n++)
            {
                float32_t v_re, z_re;
                A_re = R[m*N0 + n];
                v_re = v[M*l + m];
                z_re = Z[n];
                A_re -= (z_re*v_re);
                R[m*N0 + n] = A_re;
            }
        }
        R += SA;
        Z += N;
    }
#endif // 0

    int m, l, n;
    int N_;

    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pR0;
    const xb_vecN_2xf32 * restrict pR1;
    const xb_vecN_2xf32 * restrict pR_ = (const xb_vecN_2xf32 *)R;
    const xtfloat       * restrict pV;
    const xtfloat       * restrict pV0 = (const xtfloat       *)v;
          xb_vecN_2xf32 * restrict pZ  = (      xb_vecN_2xf32 *)Z;
          xb_vecN_2xf32 * restrict pRw;

    xb_vecN_2xf32 Acc, Acc2, V0, R0, Z0;
    valign vR, vR1;

    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N0 % 4 == 0);

    if (N0 & 4)
    {
        for (l = 0; l < L; l++)
        {
            pR0 = pR_;
            for (n = N; n > 0; n -= (BBE_SIMD_WIDTH / 2))
            {
                pV = pV0;
                pR = pR0;
                pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(N0, (uintptr_t)pR0);
                Acc = BBE_ZERON_2XF32();
                Acc2 = BBE_ZERON_2XF32();
                for (m = 0; m < M >> 1; m++)
                {
                    BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
                    V0 = BBE_REPN_2XF32(V0, 0);
                    vR = BBE_LAN_2XF32_PP(pR);
                    BBE_LAN_2XF32_IP(R0, vR, pR);
                    BBE_MULAN_2XF32(Acc, V0, R0);

                    pR = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0 - 8, (uintptr_t)pR);

                    BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
                    V0 = BBE_REPN_2XF32(V0, 0);
                    vR1 = BBE_LAN_2XF32_PP(pR1);
                    BBE_LAN_2XF32_IP(R0, vR1, pR1);
                    BBE_MULAN_2XF32(Acc2, V0, R0);

                    pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0 - 8, (uintptr_t)pR1);
                }
                if (M & 1)
                {
                    BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                    V0 = BBE_REPN_2XF32(V0, 0);
                    vR = BBE_LAN_2XF32_PP(pR);
                    BBE_LAN_2XF32_IP(R0, vR, pR);
                    BBE_MULAN_2XF32(Acc, V0, R0);
                }
                Acc = BBE_ADDN_2XF32(Acc, Acc2);
                Z0 = BBE_MULN_2XF32(Acc, 2.f);
                BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);

                pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pR0);
            }
            pV0 = (const xtfloat *)XT_ADDX4(M, (uintptr_t)pV0);
            pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
        }

        N_ = (N + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
        pR_ = (const xb_vecN_2xf32 *)R;
        pV = (const xtfloat *)v;
        pZ = (xb_vecN_2xf32 *)Z;
        for (l = 0; l < L; l++)
        {
            pR0 = pR_;
            for (m = 0; m < M; m++)
            {
                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);

                pR = pR0;
                pRw = (xb_vecN_2xf32 *)pR0;
                vR = BBE_LAN_2XF32_PP(pR);
                vR1 = BBE_ZALIGN();

                for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
                {
                    BBE_LAN_2XF32_IP(Acc, vR, pR);
                    BBE_LVN_2XF32_XP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                    BBE_MULSN_2XF32(Acc, Z0, V0);
                    BBE_SAN_2XF32_IP(Acc, vR1, pRw);
                }
                n = N % (BBE_SIMD_WIDTH / 2);
                if (n)
                {
                    BBE_LAN_2XF32_IP(Acc, vR, pR);
                    BBE_LVN_2XF32_XP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                    BBE_MULSN_2XF32(Acc, Z0, V0);
                    BBE_SAVN_2XF32_XP(Acc, vR1, pRw, n*sizeof(float32_t));
                }
                BBE_SAN_2XF32POS_FP(vR1, pRw);

                pZ = (xb_vecN_2xf32 *)XT_ADDX4(-N_, (uintptr_t)pZ);
                pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(N0, (uintptr_t)pR0);
            }
            pZ = (xb_vecN_2xf32 *)XT_ADDX4(N_, (uintptr_t)pZ);
            pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
        }

        return;
    }

    pR_ = (const xb_vecN_2xf32 *)((uintptr_t)R & ~(2 * BBE_SIMD_WIDTH - 1));
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n += (BBE_SIMD_WIDTH / 2))
        {
            pV = (const xtfloat *)(v + M*l);
            pR = (const xb_vecN_2xf32 *)((float32_t*)pR_ + l*SA + n);
            pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(N0, (uintptr_t)pR);
            Acc = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            for (m = 0; m < M >> 1; m++)
            {
                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                BBE_LVN_2XF32_XP(R0, pR, 2 * N0 * sizeof(float32_t));
                BBE_MULAN_2XF32(Acc, V0, R0);

                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                BBE_LVN_2XF32_XP(R0, pR1, 2 * N0 * sizeof(float32_t));
                BBE_MULAN_2XF32(Acc2, V0, R0);
            }
            if (M & 1)
            {
                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                BBE_LVN_2XF32_XP(R0, pR, 2 * N0 * sizeof(float32_t));
                BBE_MULAN_2XF32(Acc, V0, R0);
            }
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Z0 = BBE_MULN_2XF32(Acc, 2.f);

            pV = (const xtfloat *)(v + M*l);
            pR = (const xb_vecN_2xf32 *)((float32_t*)pR_ + l*SA + n);
            pRw = (xb_vecN_2xf32 *)(pR);

            for (m = 0; m < M; m++)
            {
                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                BBE_LVN_2XF32_XP(Acc, pR, N0 * sizeof(float32_t));

                BBE_MULSN_2XF32(Acc, Z0, V0);

                BBE_SVN_2XF32_XP(Acc, pRw, N0 * sizeof(float32_t));
            }
        }
    }
}
#endif

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
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "cqrnf_common.h"
#include "common.h"
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
    Z[2*N*L]
-------------------------------------------------------*/
void cqrnfUpdateR(            float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L)
{
#if 0
    int m, l, n;
    float32_t A_re, A_im;
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N0 % 4 == 0);
    // compute v'B first
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
        {
            A_re = A_im = 0.f;
            for (m = 0; m < M; m++)
            {
                float32_t v_re, v_im, b_re, b_im;
                v_re = v[2 * M*l + 2 * m + 0];  v_im = v[2 * M*l + 2 * m + 1];
                b_re = R[2 * m*N0 + 2 * n + 0]; b_im = R[2 * m*N0 + 2 * n + 1];
                A_re += (v_re*b_re) + (v_im*b_im);
                A_im += (v_re*b_im) - (v_im*b_re);
            }
            Z[2 * n + 0] = 2.f*A_re;
            Z[2 * n + 1] = 2.f*A_im;
        }
        R += SA;
        Z += 2 * N;
    }
    R -= SA*L; Z -= 2 * N*L;
    for (l = 0; l < L; l++)
    {
        for (m = 0; m < M; m++)
        {
            for (n = 0; n < N; n++)
            {
                float32_t v_re, v_im, z_re, z_im;
                A_re = R[2 * m*N0 + 2 * n + 0]; A_im = R[2 * m*N0 + 2 * n + 1];
                v_re = v[2 * M*l + 2 * m + 0]; v_im = v[2 * M*l + 2 * m + 1];
                z_re = Z[2 * n + 0]; z_im = Z[2 * n + 1];
                A_re -= (z_re*v_re) - (z_im*v_im);
                A_im -= (z_re*v_im) + (z_im*v_re);
                R[2 * m*N0 + 2 * n + 0] = A_re;
                R[2 * m*N0 + 2 * n + 1] = A_im;
            }
        }
        R += SA;
        Z += 2 * N;
    }
#endif // 0

    int m, l, n;
    int N_;

    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pR0;
    const xb_vecN_2xf32 * restrict pR1;
    const xb_vecN_2xf32 * restrict pR_ = (const xb_vecN_2xf32 *)R;
    const long long     * restrict pV;
    const long long     * restrict pV0 = (const long long     *)v;
          xb_vecN_2xf32 * restrict pZ  = (      xb_vecN_2xf32 *)Z;
          xb_vecN_2xf32 * restrict pRw;

    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, R0, Z0;
    valign vR, vR1;

    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N0 % 4 == 0);

    if (N0 & 4)
    {
        for (l = 0; l < L; l++)
        {
            pR0 = pR_;
            for (n = N; n > 0; n -= (BBE_SIMD_WIDTH / 4))
            {
                pV = pV0;
                pR = pR0;
                pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0, (uintptr_t)pR0);
                Acc = BBE_ZERON_2XF32();
                Acc1 = BBE_ZERON_2XF32();
                Acc2 = BBE_ZERON_2XF32();
                Acc3 = BBE_ZERON_2XF32();
                for (m = 0; m < M >> 1; m++)
                {
                    BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                    V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                    vR = BBE_LAN_2XF32_PP(pR);
                    BBE_LAN_2XF32_IP(R0, vR, pR);
                    BBE_MULMASN_2XF32(Acc, V0, R0, 0, 4);
                    BBE_MULMASN_2XF32(Acc1, V0, R0, 2, 11);

                    pR = (const xb_vecN_2xf32 *)XT_ADDX4(4 * N0 - 8, (uintptr_t)pR);

                    BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                    V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                    vR1 = BBE_LAN_2XF32_PP(pR1);
                    BBE_LAN_2XF32_IP(R0, vR1, pR1);
                    BBE_MULMASN_2XF32(Acc2, V0, R0, 0, 4);
                    BBE_MULMASN_2XF32(Acc3, V0, R0, 2, 11);

                    pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(4 * N0 - 8, (uintptr_t)pR1);
                }
                if (M & 1)
                {
                    BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
                    V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                    vR = BBE_LAN_2XF32_PP(pR);
                    BBE_LAN_2XF32_IP(R0, vR, pR);
                    BBE_MULMASN_2XF32(Acc, V0, R0, 0, 4);
                    BBE_MULMASN_2XF32(Acc1, V0, R0, 2, 11);
                }
                Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
                Acc = BBE_ADDN_2XF32(Acc, Acc1);
                Acc = BBE_ADDN_2XF32(Acc, Acc2);
                Z0 = BBE_MULN_2XF32(Acc, 2.f);
                BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);

                pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pR0);
            }
            pV0 = (const long long *)XT_ADDX4(2 * M, (uintptr_t)pV0);
            pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
        }

        N_ = (N + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
        pR_ = (const xb_vecN_2xf32 *)R;
        pV = (const long long     *)v;
        pZ = (xb_vecN_2xf32 *)Z;
        for (l = 0; l < L; l++)
        {
            pR0 = pR_;
            for (m = 0; m < M; m++)
            {
                BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

                pR = pR0;
                pRw = (xb_vecN_2xf32 *)pR0;
                vR = BBE_LAN_2XF32_PP(pR);
                vR1 = BBE_ZALIGN();

                for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
                {
                    BBE_LAN_2XF32_IP(Acc, vR, pR);
                    BBE_LVN_2XF32_XP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                    BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
                    BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
                    BBE_SAN_2XF32_IP(Acc, vR1, pRw);
                }
                n = N % (BBE_SIMD_WIDTH / 4);
                if (n)
                {
                    BBE_LAN_2XF32_IP(Acc, vR, pR);
                    BBE_LVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                    BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
                    BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
                    BBE_SAVN_2XF32_XP(Acc, vR1, pRw, 2 * n*sizeof(float32_t));
                }

                BBE_SAN_2XF32POS_FP(vR1, pRw);

                pZ = (xb_vecN_2xf32 *)XT_ADDX4(-2 * N_, (uintptr_t)pZ);
                pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0, (uintptr_t)pR0);
            }
            pZ = (xb_vecN_2xf32 *)XT_ADDX4(2 * N_, (uintptr_t)pZ);
            pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
        }

        return;
    }

    pR_ = (const xb_vecN_2xf32 *)((uintptr_t)R & ~(2 * BBE_SIMD_WIDTH - 1));
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n += (BBE_SIMD_WIDTH / 4))
        {
            pV = (const long long *)(v + 2 * M*l);
            pR = (const xb_vecN_2xf32 *)((float32_t*)pR_ + l*SA + 2 * n);
            pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0, (uintptr_t)pR);
            Acc = BBE_ZERON_2XF32();
            Acc1 = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            Acc3 = BBE_ZERON_2XF32();
            for (m = 0; m < M >> 1; m++)
            {
                BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(R0, pR, 4 * N0 * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, V0, R0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, V0, R0, 2, 11);

                BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(R0, pR1, 4 * N0 * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc2, V0, R0, 0, 4);
                BBE_MULMASN_2XF32(Acc3, V0, R0, 2, 11);
            }
            if (M & 1)
            {
                BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(R0, pR, 4 * N0 * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, V0, R0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, V0, R0, 2, 11);
            }
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Z0 = BBE_MULN_2XF32(Acc, 2.f);

            pV = (const long long *)(v + 2 * M*l);
            pR = (const xb_vecN_2xf32 *)((float32_t*)pR_ + l*SA + 2 * n);
            pRw = (xb_vecN_2xf32 *)(pR);

            for (m = 0; m < M; m++)
            {
                BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(Acc, pR, 2 * N0 * sizeof(float32_t));

                BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
                BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);

                BBE_SVN_2XF32_XP(Acc, pRw, 2 * N0 * sizeof(float32_t));
            }
        }
    }
}
#endif

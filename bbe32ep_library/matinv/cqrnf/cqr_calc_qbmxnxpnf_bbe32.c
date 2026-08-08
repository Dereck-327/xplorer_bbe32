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
#include "cqrnf_common.h"
#if HAVE_VFPU

/*
Reference code:
% compute Q'B matrix
% input:
% V  - sequence of Housholder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% R  - upper triangle decomposition
function [B] = cqr_calcQB(B,V,Fi)
[M, P] = size(B); 
[N, t] = size(Fi);
Z=zeros(M,P);
im=1;
for m=1:N
    v=V(im:im+M-m);
    im=im+(M-m+1);
    Bm=B(m:end,:);
    Bm=(Bm-2*v*v'*Bm);
    B(m:end,:)=Bm;
end
B=diag([Fi;ones(M-N,1)])'*B;
*/

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    m=XT_MIN(m,LOG2_BBE_SIMD_WIDTH-1);
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}


/*------------------------------------------------
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    introduces additional shift right !
    Input:
    Fi[L][SV]
    Input/output:
    B[L][SB]
    Temporary:
    pScr    - size in bytes N*2*BBE_SIMD_WIDTH
------------------------------------------------*/
static void cqrnfUpdateB(float32_t *Z, float32_t* B,const float32_t* V,int M,int N,int P,int SB,int L)
{
#if 0
    int m, l, p;
    float32_t A_re, A_im;
    // compute v'B first
    for (l = 0; l < L; l++)
        for (p = 0; p < P; p++)
        {
            A_re = A_im = 0.f;
            for (m = 0; m < M; m++)
            {
                float32_t v_re, v_im, b_re, b_im;
                v_re = V[2 * M*l + 2 * m + 0];        v_im = V[2 * M*l + 2 * m + 1];
                b_re = B[SB*l + 2 * m*P + 2 * p + 0];  b_im = B[SB*l + 2 * m*P + 2 * p + 1];
                A_re += (v_re*b_re) + (v_im*b_im);
                A_im += (v_re*b_im) - (v_im*b_re);
            }
            Z[2 * P*l + 2 * p + 0] = 2.f*A_re;
            Z[2 * P*l + 2 * p + 1] = 2.f*A_im;
        }

    for (l = 0; l < L; l++)
        for (p = 0; p < P; p++)
            for (m = 0; m < M; m++)
            {
                float32_t v_re, v_im, z_re, z_im;
                A_re = B[SB*l + 2 * m*P + 2 * p + 0]; A_im = B[SB*l + 2 * m*P + 2 * p + 1];
                v_re = V[2 * M*l + 2 * m + 0]; v_im = V[2 * M*l + 2 * m + 1];
                z_re = Z[2 * P*l + 2 * p + 0]; z_im = Z[2 * P*l + 2 * p + 1];
                A_re -= (z_re*v_re) - (z_im*v_im);
                A_im -= (z_re*v_im) + (z_im*v_re);
                B[SB*l + 2 * m*P + 2 * p + 0] = A_re;
                B[SB*l + 2 * m*P + 2 * p + 1] = A_im;
            }
#endif // 0

    int m, l, p;
    int count, count2, P_, M_;
    int test;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pB0;
    const xb_vecN_2xf32 * restrict pB1;
    const xb_vecN_2xf32 * restrict pB_;
    const xb_vecN_2xf32 * restrict pV_;
    const long long     * restrict pV;
    const long long     * restrict pV0;
          xb_vecN_2xf32 * restrict pZ;
          xb_vecN_2xf32 * restrict pBw;
    
    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, B0, Z0;
    valign vB, vB1, vV;

    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);

    if (P == 1)
    {
        count = (M % (BBE_SIMD_WIDTH / 2)) * 2 * sizeof(float32_t);
        pZ = (xb_vecN_2xf32 *)Z;
        pB_ = (const xb_vecN_2xf32 *)B;
        pV_ = (const xb_vecN_2xf32 *)V;
        vV = BBE_LAN_2XF32_PP(pV_);
        for (l = 0; l < L; l++)
        {
            pB = pB_;
            vB = BBE_LAN_2XF32_PP(pB);

            BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
            BBE_LAVN_2XF32_XP(B0, vB, pB, count);
            Acc = BBE_MULMN_2XF32(V0, B0, 0, 4);
            Acc1 = BBE_MULMN_2XF32(V0, B0, 2, 11);

            BBE_LAVN_2XF32_XP(V0, vV, pV_, count - 2 * BBE_SIMD_WIDTH);
            BBE_LAVN_2XF32_XP(B0, vB, pB, count - 2 * BBE_SIMD_WIDTH);
            Acc2 = BBE_MULMN_2XF32(V0, B0, 0, 4);
            Acc3 = BBE_MULMN_2XF32(V0, B0, 2, 11);
            for (m = 0; m < M / (BBE_SIMD_WIDTH / 2); m++)
            {
                BBE_LAN_2XF32_IP(V0, vV, pV_);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULMASN_2XF32(Acc, V0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, V0, B0, 2, 11);

                BBE_LAN_2XF32_IP(V0, vV, pV_);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULMASN_2XF32(Acc2, V0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc3, V0, B0, 2, 11);
            }


            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
            Acc = BBE_MULN_2XF32(Acc, 2.f);
            Z0 = BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_REP_0X4);
            BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);

            pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
        }

        /*pZ = (xb_vecN_2xf32 *)Z;
        pB_ = (const xb_vecN_2xf32 *)B;
        pV_ = (const xb_vecN_2xf32 *)V;
        vV = BBE_LAN_2XF32_PP(pV_);
        for (l = 0; l < L; l++)
        {
            pB = pB_;
            vB = BBE_LAN_2XF32_PP(pB);
            pBw = (xb_vecN_2xf32 *)pB_;
            vB1 = BBE_ZALIGN();
            BBE_LVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
            for (m = M; m > 0; m -= (BBE_SIMD_WIDTH / 4))
            {
                BBE_LAN_2XF32_IP(Acc, vB, pB);
                BBE_LAVN_2XF32_XP(V0, vV, pV_, 2 * m*sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
                BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
                BBE_SAVN_2XF32_XP(Acc, vB1, pBw, 2 * m*sizeof(float32_t));
            }
            BBE_SAN_2XF32POS_FP(vB1, pBw);


            pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
        }*/
        M_ = (M + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
        count = M * 2 * sizeof(float32_t);
        count2 = M_ * 2 * sizeof(float32_t);
        pZ = (xb_vecN_2xf32 *)Z;
        pB = (const xb_vecN_2xf32 *)B;
        pBw = (xb_vecN_2xf32 *)B;
        pV_ = (const xb_vecN_2xf32 *)V;
        vV = BBE_LAN_2XF32_PP(pV_);
        for (l = 0; l < L * M_ / (BBE_SIMD_WIDTH / 4); l++)
        {
            vB = BBE_LAN_2XF32_PP(pB);
            vB1 = BBE_ZALIGN();

            Z0 = BBE_LVN_2XF32_I(pZ, 0);
            BBE_LAVN_2XF32_XP(Acc, vB, pB, count);
            BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, count);
            BBE_SAN_2XF32POS_FP(vB1, pBw);

            count = XT_ADDX4(-8, count);
            count2 = XT_ADDX4(-8, count2);
            test = 0;
            XT_MOVEQZ(test, 8, count2);
            pZ = (xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pZ);
            test = 0;
            XT_MOVEQZ(test, SB - 2 * M, count2);
            pB = (const xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pB);
            pBw = (xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pBw);
            XT_MOVEQZ(count, M * 2 * sizeof(float32_t), count2);
            XT_MOVEQZ(count2, M_ * 2 * sizeof(float32_t), count2);
        }

        return;
    }


    pB_ = (const xb_vecN_2xf32 *)B;
    pV0 = (const long long *)V;
    pZ = (xb_vecN_2xf32 *)Z;
    for (l = 0; l < L; l++)
    {
        pB0 = pB_;
        for (p = P; p > 0; p -= (BBE_SIMD_WIDTH / 4))
        {
            pV = pV0;
            pB = pB0;
            pB1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pB0);
            Acc = BBE_ZERON_2XF32();
            Acc1 = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            Acc3 = BBE_ZERON_2XF32();
            for (m = 0; m < M >> 1; m++)
            {
                BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                vB = BBE_LAN_2XF32_PP(pB);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULMASN_2XF32(Acc, V0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, V0, B0, 2, 11);

                pB = (const xb_vecN_2xf32 *)XT_ADDX4(4 * P - 8, (uintptr_t)pB);

                BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                vB1 = BBE_LAN_2XF32_PP(pB1);
                BBE_LAN_2XF32_IP(B0, vB1, pB1);
                BBE_MULMASN_2XF32(Acc2, V0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc3, V0, B0, 2, 11);

                pB1 = (const xb_vecN_2xf32 *)XT_ADDX4(4 * P - 8, (uintptr_t)pB1);
            }
            if (M & 1)
            {
                BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
                V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
                vB = BBE_LAN_2XF32_PP(pB);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULMASN_2XF32(Acc, V0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc1, V0, B0, 2, 11);
            }
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Z0 = BBE_MULN_2XF32(Acc, 2.f);
            BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);

            pB0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pB0);
        }
        pV0 = (const long long *)XT_ADDX4(2 * M, (uintptr_t)pV0);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }

    /*P_ = (P + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
    pB_ = (const xb_vecN_2xf32 *)B;
    pV = (const long long *)V;
    pZ = (xb_vecN_2xf32 *)Z;
    for (l = 0; l < L; l++)
    {
        pB0 = pB_;

        pB = pB0;
        pBw = (xb_vecN_2xf32 *)pB0;
        vB = BBE_LAN_2XF32_PP(pB);
        vB1 = BBE_ZALIGN();
        for (m = 0; m < M; m++)
        {
            BBE_LSN_4X64_IP(vTmp, pV, 2 * sizeof(float32_t));
            V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            for (p = 0; p < P / (BBE_SIMD_WIDTH / 4); p++)
            {
                BBE_LAN_2XF32_IP(Acc, vB, pB);
                BBE_LVN_2XF32_XP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
                BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
                BBE_SAN_2XF32_IP(Acc, vB1, pBw);
            }
            BBE_LAVN_2XF32_XP(Acc, vB, pB, 2 * (P % (BBE_SIMD_WIDTH / 4))*sizeof(float32_t));
            BBE_LVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, 2 * (P % (BBE_SIMD_WIDTH / 4))*sizeof(float32_t));

            pZ = (xb_vecN_2xf32 *)XT_ADDX4(-2 * P_, (uintptr_t)pZ);
        }
        BBE_SAN_2XF32POS_FP(vB1, pBw);
        pZ = (xb_vecN_2xf32 *)XT_ADDX4(2 * P_, (uintptr_t)pZ);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }*/
    P_ = (P + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
    count = P * 2 * sizeof(float32_t);
    count2 = P_ * 2 * sizeof(float32_t);
    pB_ = (const xb_vecN_2xf32 *)B;
    pV = (const long long *)V;
    pZ = (xb_vecN_2xf32 *)Z;
    for (l = 0; l < L; l++)
    {
        WUR_CBEGIN((uintptr_t)pZ);
        WUR_CEND((uintptr_t)pZ + 2 * P_ * sizeof(float32_t));

        pB = pB_;
        pBw = (xb_vecN_2xf32 *)pB_;
        vB = BBE_LAN_2XF32_PP(pB);
        vB1 = BBE_ZALIGN();
        for (m = 0; m < M * P_ / (BBE_SIMD_WIDTH / 4); m++)
        {
            vTmp = BBE_LSN_4X64_I(pV, 0);
            V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            BBE_LAVN_2XF32_XP(Acc, vB, pB, count);
            BBE_LVN_2XF32_IC(Z0, pZ);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, count);

            count = XT_ADDX4(-8, count);
            count2 = XT_ADDX4(-8, count2);
            test = 0;
            XT_MOVEQZ(test, 2, count2);
            pV = (const long long *)XT_ADDX4(test, (uintptr_t)pV);
            XT_MOVEQZ(count, P * 2 * sizeof(float32_t), count2);
            XT_MOVEQZ(count2, P_ * 2 * sizeof(float32_t), count2);
        }
        BBE_SAN_2XF32POS_FP(vB1, pBw);
        pZ = (xb_vecN_2xf32 *)XT_ADDX4(2 * P_, (uintptr_t)pZ);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}

/*
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
*/
static void cqrnfRotateB(float32_t* B,const float32_t* Fi,int N,int P,int SB,int L)
{
#if 0
    int n, p, l;
    float32_t A_re, A_im;

    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
            for (p = 0; p < P; p++)
            {
                float32_t f_re, f_im, b_re, b_im;
                f_re = Fi[2 * l + 2 * n*L + 0]; f_im = Fi[2 * l + 2 * n*L + 1];
                b_re = B[2 * n*P + 2 * p + 0]; b_im = B[2 * n*P + 2 * p + 1];
                A_re = (f_re*b_re) + (f_im*b_im);
                A_im = (f_re*b_im) - (f_im*b_re);
                B[2 * n*P + 2 * p + 0] = A_re;
                B[2 * n*P + 2 * p + 1] = A_im;
            }
        B += SB;
    }
#endif // 0

    int n, p, l;
    int count = (P % (BBE_SIMD_WIDTH / 4)) * 2 * sizeof(float32_t);

    const long long     * restrict pF;
    const long long     * restrict pF_;
    const xb_vecN_2xf32 * restrict pBr0;
    const xb_vecN_2xf32 * restrict pBr1;
          xb_vecN_2xf32 * restrict pBw0;
          xb_vecN_2xf32 * restrict pBw1;
          xb_vecN_2xf32 * restrict pB_;

    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 Acc, F0, F1, B0;
    valign vBr0, vBr1, vBw0, vBw1;

    NASSERT(N % 4 == 0);

    pF_ = (const long long *)Fi;
    pB_ = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        pF = pF_;
        pBw0 = pB_;
        pBw1 = (xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pB_);
        pBr0 = pB_;
        pBr1 = pBw1;
        for (n = 0; n < N >> 1; n++)
        {
            vBw0 = BBE_ZALIGN();
            vBw1 = BBE_ZALIGN();
            vBr0 = BBE_LAN_2XF32_PP(pBr0);
            vBr1 = BBE_LAN_2XF32_PP(pBr1);

            BBE_LSN_4X64_XP(vTmp, pF, 2 * L*sizeof(float32_t));
            F0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            BBE_LSN_4X64_XP(vTmp, pF, 2 * L*sizeof(float32_t));
            F1 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            for (p = 0; p < P / (BBE_SIMD_WIDTH / 4); p++)
            {
                BBE_LAN_2XF32_IP(B0, vBr0, pBr0);
                Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
                BBE_SAN_2XF32_IP(Acc, vBw0, pBw0);

                BBE_LAN_2XF32_IP(B0, vBr1, pBr1);
                Acc = BBE_MULMN_2XF32(F1, B0, 0, 4);
                BBE_MULMASN_2XF32(Acc, F1, B0, 2, 11);
                BBE_SAN_2XF32_IP(Acc, vBw1, pBw1);
            }
            BBE_LAVN_2XF32_XP(B0, vBr0, pBr0, count);
            Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
            BBE_SAVN_2XF32_XP(Acc, vBw0, pBw0, count);

            BBE_LAVN_2XF32_XP(B0, vBr1, pBr1, count);
            Acc = BBE_MULMN_2XF32(F1, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F1, B0, 2, 11);
            BBE_SAVN_2XF32_XP(Acc, vBw1, pBw1, count);

            BBE_SAN_2XF32POS_FP(vBw0, pBw0);
            BBE_SAN_2XF32POS_FP(vBw1, pBw1);

            pBw0 = (xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pBw0);
            pBw1 = (xb_vecN_2xf32 *)XT_ADDX4(2 * P, (uintptr_t)pBw1);
            pBr0 = pBw0;
            pBr1 = pBw1;
        }

        pF_ = (const long long *)XT_ADDX4(2, (uintptr_t)pF_);
        pB_ = (xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}

/*-------------------------------------------------------------------------
Update right side of equations for QR process for block ordered matrices.
Matrix sizes SB,SV are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2.  
SB=size(M*P)
SV=size(((2*M-N+1)*N/2+N)*L)
Scratch size in bytes is defined by cqr_calc_qbmxnn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      dimensional parameters
 L            Number of matrices
Input/output:
 B[L][SB]     On input it is the sequence of L matrices B. 
              At the end of the process, matrices Z replace input
              matrices A. In a case of non-square matrices (N!=M), 
              only N*P elements of each output matrix will be valid.
Input:
 V[SV]        Sequence of L Housholder rotation vectors 

Restrictions:
1. B, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_calc_qbmxnxpnf(void *pScr,
                    complex_float* _B,const complex_float* _V,int M, int N, int P,int L)
{
    float32_t*       B=(float32_t*      )_B;
    const float32_t* V=(const float32_t*)_V;
    float32_t* Z=(float32_t*)pScr;
    int m;
    int SB=getSpace(M*P<<1);
    const float32_t* pV;

    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(M%4==0 && N%4==0);
    NASSERT(M>0 && N>0);
    NASSERT(N<=M);
    if (L<=0 || M<=0 || N<=0 || N>M) return;
    for (pV=V,m=0; m<N; m++)
    {
        cqrnfUpdateB(Z,B+2*m*P,pV,(M-m),N,P,SB,L); 
        pV+=2*(M-m)*L;
    }
    if (P==1) cqrnfRotateB1(pScr,B,V+(2*M-N+1)*N*L,N,SB,L);
    else      cqrnfRotateB(B,V+(2*M-N+1)*N*L,N,P,SB,L);
}

/* scratch memory needed for calc_qb functions */
size_t cqr_calc_qbmxnxpnf_getScratchSize   (int M, int N,int P,int L)
{
    size_t Zsize,Bsize;
    NASSERT(L>0);
    NASSERT(M%4==0 && N%4==0);
    NASSERT(M>0 && N>0);
    NASSERT(N<=M);
    (void)M,(void)N,(void)P,(void)L;
    L=XT_MAX(L,0);
    L = (L + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
    P = (P + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1);
    Zsize=P*L*sizeof(complex_float);
    Bsize=2*N* (2*BBE_SIMD_WIDTH);
    return XT_MAX(Zsize,Bsize);
}

#else
DISCARD_FUN(void, cqr_calc_qbmxnxpnf, (void *pScr,
                    complex_float* _B,const complex_float* _V,int M, int N, int P,int L))

size_t cqr_calc_qbmxnxpnf_getScratchSize   (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(M%4==0 && N%4==0);
    NASSERT(M>0 && N>0);
    NASSERT(N<=M);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

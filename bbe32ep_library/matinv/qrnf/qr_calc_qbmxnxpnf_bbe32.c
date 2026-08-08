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
#include "qrnf_common.h"
#include "common.h"

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

/*------------------------------------------------
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV]
    Input/output:
    B[L][SB]
    Temporary:
    pScr    - size in bytes N*2*BBE_SIMD_WIDTH
------------------------------------------------*/
static void qrnfUpdateB(float32_t *Z, float32_t* B,const float32_t* V,int M,int N,int P,int SB,int L)
{
#if 0
    int m, l, p;
    float32_t A_re;
    // compute v'B first
    for (l = 0; l < L; l++)
    {
        for (p = 0; p < P; p++)
        {
            A_re = 0.f;
            for (m = 0; m < M; m++)
            {
                float32_t v_re, b_re;
                v_re = V[m];
                b_re = B[m*P + p];
                A_re += 2.f*v_re*b_re;
            }
            Z[p] = A_re;
        }
        V += M;
        B += SB;
        Z += P;
    }
    V -= M*L; B -= SB*L; Z -= P*L;
    for (l = 0; l < L; l++)
    {
        for (m = 0; m < M; m++)
        {
            for (p = 0; p < P; p++)
            {
                float32_t v_re, z_re;
                A_re = B[m*P + p];
                v_re = V[m];
                z_re = Z[p];
                A_re -= z_re*v_re;
                B[m*P + p] = A_re;
            }
        }
        V += M;
        B += SB;
        Z += P;
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
    const xtfloat       * restrict pV;
    const xtfloat       * restrict pV0;
          xb_vecN_2xf32 * restrict pZ;
          xb_vecN_2xf32 * restrict pBw;
    
    xb_vecN_2xf32 Acc, Acc2, V0, B0, Z0;
    valign vB, vB1, vV;

    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);

    if (P == 1)
    {
        count = (M % (BBE_SIMD_WIDTH)) * sizeof(float32_t);
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
            Acc = BBE_MULN_2XF32(V0, B0);

            BBE_LAVN_2XF32_XP(V0, vV, pV_, count - 2 * BBE_SIMD_WIDTH);
            BBE_LAVN_2XF32_XP(B0, vB, pB, count - 2 * BBE_SIMD_WIDTH);
            Acc2 = BBE_MULN_2XF32(V0, B0);
            for (m = 0; m < M / (BBE_SIMD_WIDTH); m++)
            {
                BBE_LAN_2XF32_IP(V0, vV, pV_);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULAN_2XF32(Acc, V0, B0);

                BBE_LAN_2XF32_IP(V0, vV, pV_);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULAN_2XF32(Acc2, V0, B0);
            }


            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
            Acc = BBE_MULN_2XF32(Acc, 2.f);
            Z0 = BBE_REPN_2XF32(Acc, 0);
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
            for (m = M; m > 0; m -= (BBE_SIMD_WIDTH / 2))
            {
                BBE_LAN_2XF32_IP(Acc, vB, pB);
                BBE_LAVN_2XF32_XP(V0, vV, pV_, m*sizeof(float32_t));
                BBE_MULSN_2XF32(Acc, Z0, V0);
                BBE_SAVN_2XF32_XP(Acc, vB1, pBw, m*sizeof(float32_t));
            }
            BBE_SAN_2XF32POS_FP(vB1, pBw);


            pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
        }*/
        M_ = (M + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
        count = M * sizeof(float32_t);
        count2 = M_ * sizeof(float32_t);
        pZ = (xb_vecN_2xf32 *)Z;
        pB = (const xb_vecN_2xf32 *)B;
        pBw = (xb_vecN_2xf32 *)B;
        pV_ = (const xb_vecN_2xf32 *)V;
        vV = BBE_LAN_2XF32_PP(pV_);
        for (l = 0; l < L * M_ / (BBE_SIMD_WIDTH / 2); l++)
        {
            vB = BBE_LAN_2XF32_PP(pB);
            vB1 = BBE_ZALIGN();

            Z0 = BBE_LVN_2XF32_I(pZ, 0);
            BBE_LAVN_2XF32_XP(Acc, vB, pB, count);
            BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
            BBE_MULSN_2XF32(Acc, Z0, V0);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, count);
            BBE_SAN_2XF32POS_FP(vB1, pBw);

            count = XT_ADDX4(-8, count);
            count2 = XT_ADDX4(-8, count2);
            test = 0;
            XT_MOVEQZ(test, 8, count2);
            pZ = (xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pZ);
            test = 0;
            XT_MOVEQZ(test, SB - M, count2);
            pB = (const xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pB);
            pBw = (xb_vecN_2xf32 *)XT_ADDX4(test, (uintptr_t)pBw);
            XT_MOVEQZ(count, M * sizeof(float32_t), count2);
            XT_MOVEQZ(count2, M_ * sizeof(float32_t), count2);
        }

        return;
    }


    pB_ = (const xb_vecN_2xf32 *)B;
    pV0 = (const xtfloat *)V;
    pZ = (xb_vecN_2xf32 *)Z;
    for (l = 0; l < L; l++)
    {
        pB0 = pB_;
        for (p = P; p > 0; p -= (BBE_SIMD_WIDTH / 2))
        {
            pV = pV0;
            pB = pB0;
            pB1 = (const xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pB0);
            Acc = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            for (m = 0; m < M >> 1; m++)
            {
                BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                vB = BBE_LAN_2XF32_PP(pB);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULAN_2XF32(Acc, V0, B0);

                pB = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P - 8, (uintptr_t)pB);

                BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                vB1 = BBE_LAN_2XF32_PP(pB1);
                BBE_LAN_2XF32_IP(B0, vB1, pB1);
                BBE_MULAN_2XF32(Acc2, V0, B0);

                pB1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * P - 8, (uintptr_t)pB1);
            }
            if (M & 1)
            {
                BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
                V0 = BBE_REPN_2XF32(V0, 0);
                vB = BBE_LAN_2XF32_PP(pB);
                BBE_LAN_2XF32_IP(B0, vB, pB);
                BBE_MULAN_2XF32(Acc, V0, B0);
            }
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            Z0 = BBE_MULN_2XF32(Acc, 2.f);
            BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);

            pB0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pB0);
        }
        pV0 = (const xtfloat *)XT_ADDX4(M, (uintptr_t)pV0);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }

    /*P_ = (P + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
    pB_ = (const xb_vecN_2xf32 *)B;
    pV = (const xtfloat *)V;
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
            BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
            V0 = BBE_REPN_2XF32(V0, 0);

            for (p = 0; p < P / (BBE_SIMD_WIDTH / 2); p++)
            {
                BBE_LAN_2XF32_IP(Acc, vB, pB);
                BBE_LVN_2XF32_XP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
                BBE_MULSN_2XF32(Acc, Z0, V0);
                BBE_SAN_2XF32_IP(Acc, vB1, pBw);
            }
            BBE_LAVN_2XF32_XP(Acc, vB, pB, (P % (BBE_SIMD_WIDTH / 2))*sizeof(float32_t));
            BBE_LVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
            BBE_MULSN_2XF32(Acc, Z0, V0);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, (P % (BBE_SIMD_WIDTH / 2))*sizeof(float32_t));

            pZ = (xb_vecN_2xf32 *)XT_ADDX4(-P_, (uintptr_t)pZ);
        }
        BBE_SAN_2XF32POS_FP(vB1, pBw);
        pZ = (xb_vecN_2xf32 *)XT_ADDX4(P_, (uintptr_t)pZ);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }*/
    P_ = (P + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
    count = P * sizeof(float32_t);
    count2 = P_ * sizeof(float32_t);
    pB_ = (const xb_vecN_2xf32 *)B;
    pV = (const xtfloat *)V;
    pZ = (xb_vecN_2xf32 *)Z;
    for (l = 0; l < L; l++)
    {
        WUR_CBEGIN((uintptr_t)pZ);
        WUR_CEND((uintptr_t)pZ + P_ * sizeof(float32_t));

        pB = pB_;
        pBw = (xb_vecN_2xf32 *)pB_;
        vB = BBE_LAN_2XF32_PP(pB);
        vB1 = BBE_ZALIGN();
        for (m = 0; m < M * P_ / (BBE_SIMD_WIDTH / 2); m++)
        {
            V0 = BBE_LSN_2XF32_I(pV, 0);
            V0 = BBE_REPN_2XF32(V0, 0);

            BBE_LAVN_2XF32_XP(Acc, vB, pB, count);
            BBE_LVN_2XF32_IC(Z0, pZ);
            BBE_MULSN_2XF32(Acc, Z0, V0);
            BBE_SAVN_2XF32_XP(Acc, vB1, pBw, count);

            count = XT_ADDX4(-8, count);
            count2 = XT_ADDX4(-8, count2);
            test = 0;
            XT_MOVEQZ(test, 1, count2);
            pV = (const xtfloat *)XT_ADDX4(test, (uintptr_t)pV);
            XT_MOVEQZ(count, P * sizeof(float32_t), count2);
            XT_MOVEQZ(count2, P_ * sizeof(float32_t), count2);
        }
        BBE_SAN_2XF32POS_FP(vB1, pBw);
        pZ = (xb_vecN_2xf32 *)XT_ADDX4(P_, (uintptr_t)pZ);
        pB_ = (const xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}

/*
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
*/
static void qrnfRotateB(float32_t* B,const float32_t* Fi,int N,int P,int SB,int L)
{
#if 0
    int n, p, l;
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
            for (p = 0; p < P; p++)
            {
                float32_t f_re, b_re;
                f_re = Fi[l + n*L];
                b_re = B[n*P + p];
                B[n*P + p] = f_re*b_re;
            }
        B += SB;
    }
#endif // 0

    int n, p, l;
    int count = (P % (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t);

    const xtfloat       * restrict pF;
    const xtfloat       * restrict pF_;
    const xb_vecN_2xf32 * restrict pBr0;
    const xb_vecN_2xf32 * restrict pBr1;
          xb_vecN_2xf32 * restrict pBw0;
          xb_vecN_2xf32 * restrict pBw1;
          xb_vecN_2xf32 * restrict pB_;

    xb_vecN_2xf32 Acc, F0, F1, B0;
    valign vBr0, vBr1, vBw0, vBw1;

    NASSERT(N % 4 == 0);

    pF_ = (const xtfloat *)Fi;
    pB_ = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        pF = pF_;
        pBw0 = pB_;
        pBw1 = (xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pB_);
        pBr0 = pB_;
        pBr1 = pBw1;
        for (n = 0; n < N >> 1; n++)
        {
            vBw0 = BBE_ZALIGN();
            vBw1 = BBE_ZALIGN();
            vBr0 = BBE_LAN_2XF32_PP(pBr0);
            vBr1 = BBE_LAN_2XF32_PP(pBr1);

            BBE_LSN_2XF32_XP(F0, pF, L*sizeof(float32_t));
            F0 = BBE_REPN_2XF32(F0, 0);

            BBE_LSN_2XF32_XP(F1, pF, L*sizeof(float32_t));
            F1 = BBE_REPN_2XF32(F1, 0);

            for (p = 0; p < P / (BBE_SIMD_WIDTH / 2); p++)
            {
                BBE_LAN_2XF32_IP(B0, vBr0, pBr0);
                Acc = BBE_MULN_2XF32(F0, B0);
                BBE_SAN_2XF32_IP(Acc, vBw0, pBw0);

                BBE_LAN_2XF32_IP(B0, vBr1, pBr1);
                Acc = BBE_MULN_2XF32(F1, B0);
                BBE_SAN_2XF32_IP(Acc, vBw1, pBw1);
            }
            BBE_LAVN_2XF32_XP(B0, vBr0, pBr0, count);
            Acc = BBE_MULN_2XF32(F0, B0);
            BBE_SAVN_2XF32_XP(Acc, vBw0, pBw0, count);

            BBE_LAVN_2XF32_XP(B0, vBr1, pBr1, count);
            Acc = BBE_MULN_2XF32(F1, B0);
            BBE_SAVN_2XF32_XP(Acc, vBw1, pBw1, count);

            BBE_SAN_2XF32POS_FP(vBw0, pBw0);
            BBE_SAN_2XF32POS_FP(vBw1, pBw1);

            pBw0 = (xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pBw0);
            pBw1 = (xb_vecN_2xf32 *)XT_ADDX4(P, (uintptr_t)pBw1);
            pBr0 = pBw0;
            pBr1 = pBw1;
        }

        pF_ = (const xtfloat *)XT_ADDX4(1, (uintptr_t)pF_);
        pB_ = (xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}


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
void qr_calc_qbmxnxpnf(void *pScr,
                    float32_t* B,const float32_t* V,int M, int N, int P,int L)
{
    float32_t*       Z=(float32_t*)pScr;
    int m;
    int SB=getSpace(M*P);
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
        qrnfUpdateB(Z,B+m*P,pV,(M-m),N,P,SB,L); 
        pV+=(M-m)*L;
    }
    if (P==1)    qrnfRotateB1(pScr,B,V+(((2*M-N+1)*N)>>1)*L,N,SB,L);
    else         qrnfRotateB(B,V+(((2*M-N+1)*N)>>1)*L,N,P,SB,L);
}

/* scratch memory needed for calc_qb functions */
size_t qr_calc_qbmxnxpnf_getScratchSize   (int M, int N,int P,int L)
{
    size_t Zsize,Bsize;
    NASSERT(L>0);
    NASSERT(M%4==0 && N%4==0);
    NASSERT(M>0 && N>0);
    NASSERT(N<=M);
    (void)M,(void)N,(void)P,(void)L;
    L=XT_MAX(L,0);
    L = (L + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
    N = (N + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
    P = (P + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1);
    Zsize = P*L*sizeof(float32_t);
    Bsize=N* (2*BBE_SIMD_WIDTH);
    return XT_MAX(Zsize,Bsize);
}

#else
DISCARD_FUN(void, qr_calc_qbmxnxpnf, (void *pScr,
                    float32_t* B,const float32_t* V,int M, int N, int P,int L))
size_t qr_calc_qbmxnxpnf_getScratchSize   (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(M%4==0 && N%4==0);
    NASSERT(M>0 && N>0);
    NASSERT(N<=M);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

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

/*
Reference Matlab code:
% build R matrix, Housholder vectors v and diagonal rotation matrix FI
% output:
% V  - sequence of Housholder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% D  - reciprocals of main diagonal elements
% R  - upper triangle decomposition
function [V,Fi,D,Q,R] = cqr_buildr(A)
A0=A;
[M, N] = size(A); 
if (M==N) Ncolumns=N;
else      Ncolumns=N;
end

V =[];
Q = eye(M); 
Fi = zeros(N,1);
D  = zeros(N,1);

for m=1: Ncolumns
    tmpA = A(m:end, m:end);
    e1 = zeros(M-m+1,1);
    e1(1) = 1;
    x = tmpA(:,1);

    fi = x(1)/(1E-10+abs(x(1)));
    D(m)=  1/sqrt(x'*x);
    alpha = - 1/D(m) * fi;
    v = x - alpha*e1;
    Fi(m) = fi;

    v = v/sqrt(v'*v);
    V=[V;v];
    P = eye(M-m+1) - 2*v*v';
    A(m:end, m:end) = P*tmpA;
    A(m,m)=alpha;
end
Fi=-Fi; % make diagonals positive 
F = diag([Fi;ones(M-N,1)]);
Q = Q*F;
R = F'*A;
*/

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
    pScr[]    scratch, defined by cqrnHousholder_getScratchSz()
-------------------------------------------------------*/
void cqrnfHousholder(void* pScr,
                                 float32_t* restrict v,
                                 float32_t* restrict Fi,
                                 float32_t* restrict D,
                           const float32_t* restrict x, 
                           const float32_t* restrict xt, 
                           int M, int SD, int L)
{
#if 0
    float32_t *tmp = (float32_t *)pScr;   //[2*L]
    float32_t *invd = tmp + 2 * L;            //[L]
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
            float32_t  x_re, x_im;
            x_re = xt[m * 2 * L + 2 * l + 0];
            x_im = xt[m * 2 * L + 2 * l + 1];
            tmp[l] += (x_re*x_re) + (x_im*x_im);
        }
    for (l = 0; l < L; l++)
    {
        tmp[l] = 1.f / sqrtf(tmp[l]);
        invd[l] = 1.f / tmp[l];
    }
    for (l = 0; l < L; l++)
    {
        D[l*SD + 0] = D[l*SD + 1] = tmp[l];
    }
    // compute fi=x0/sqrt(x0'x0)
    for (l = 0; l < L; l++)
    {
        float32_t d;
        float32_t x_re, x_im;
        float32_t A_re;
        x_re = xt[2 * l + 0];
        x_im = xt[2 * l + 1];

        A_re = (x_re*x_re) + (x_im*x_im);
        d = 1.f / sqrtf(A_re);
        if (fabsf(A_re) > FLT_MIN)
        {
            Fi[2 * l + 0] = -x_re*d;
            Fi[2 * l + 1] = -x_im*d;
        }
        else
        {
            Fi[2 * l + 0] = 1.f;
            Fi[2 * l + 1] = 0.f;
        }
    }
    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    for (l = 0; l < L; l++) norm[l] = 0.f;
    for (l = 0; l < L; l++)
        for (m = 1; m < M; m++)
        {
            float32_t x_re, x_im;
            x_re = x[l * 2 * M + 2 * m + 0];
            x_im = x[l * 2 * M + 2 * m + 1];
            norm[l] += (x_re*x_re) + (x_im*x_im);
        }
    for (l = 0; l < L; l++)
    {
        float32_t  f_re, f_im, x_re, x_im, v_re, v_im;
        float32_t d0, invd0;
        float32_t  A_re;
        f_re = Fi[2 * l + 0];
        f_im = Fi[2 * l + 1];
        invd0 = invd[l];
        x_re = xt[2 * l + 0];
        x_im = xt[2 * l + 1];
        v_re = x_re - (f_re*invd0);
        v_im = x_im - (f_im*invd0);
        A_re = norm[l];
        A_re += (v_re*v_re) + (v_im*v_im);
        d0 = 1.f / sqrtf(A_re);
        invd[l] = d0;
        tmp[2 * l + 0] = v_re*d0;
        tmp[2 * l + 1] = v_im*d0;
    }
    for (l = 0; l < L; l++)
    {
        v[2 * M*l + 0] = tmp[2 * l + 0];
        v[2 * M*l + 1] = tmp[2 * l + 1];
    }
    for (l = 0; l < L; l++)
    {
        float32_t d0, x_re, x_im;
        d0 = invd[l];
        for (m = 1; m < M; m++)
        {
            x_re = x[l * 2 * M + 2 * m + 0];
            x_im = x[l * 2 * M + 2 * m + 1];
            v[2 * M*l + 2 * m + 0] = x_re*d0;
            v[2 * M*l + 2 * m + 1] = x_im*d0;
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
          long long     * restrict pD;
    const xtfloat       * restrict pZr_;
          xb_vecN_2xf32 * restrict pF;
          long long     * restrict pV;
          xb_vecN_2xf32 * restrict pV1;
          xtfloat       * restrict pZw_;

    xb_vecN_2xf32 Acc, Acc0, X0, X1, F0, rsqrt;
    valign vX, vZr, vZr1, vZw, vF, vV;
    xb_vecN_2xf32 r0;
    const float32_t ALIGN(32) rep_tbl[8] = { 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f };

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

        BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * sizeof(float32_t));
        Acc0 = BBE_MULN_2XF32(X0, X0);
        Acc0 = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc0, BBE_SHFLI_SWAP_2), Acc0);
        for (m = M - 1; m > 0; m -= (BBE_SIMD_WIDTH / 4))
        {
            BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * m*sizeof(float32_t));
            BBE_MULAN_2XF32(Acc, X0, X0);
        }
        Acc = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(Acc));
        Acc0 = BBE_ADDN_2XF32(Acc, Acc0);
        BBE_SSN_2XF32_X(Acc, pZw_, 3 * L_*sizeof(float32_t));
        BBE_SSN_2XF32_IP(Acc0, pZw_, sizeof(float32_t));
    }

    pZr = (const xb_vecN_2xf32 *)pScr;
    pZw = (xb_vecN_2xf32 *)((float32_t *)pScr + 2 * L_);
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
    pD = (long long *)D;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(Acc, vZr, pZr, sizeof(float32_t));
        Acc = BBE_REPN_2XF32(Acc, 0);
        BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pD, SD*sizeof(float32_t));
    }

    // compute fi=x0/sqrt(x0'x0)
    pX = (const xb_vecN_2xf32 *)xt;
    pF = (xb_vecN_2xf32 *)Fi;
    vF = BBE_ZALIGN();
    r0 = BBE_LVN_2XF32_I((const xb_vecN_2xf32 *)rep_tbl, 0); /////////////
    for (l = 0; l < L / (BBE_SIMD_WIDTH / 4); l++)
    {
        vboolN_2 tmp;

        BBE_LAN_2XF32_IP(X0, vX, pX);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        Acc = BBE_NEGN_2XF32(Acc);
        Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(r0), BBE_MOVNX16_FROMN_2XF32(Acc), BBE_MOVN_FROMN_2(tmp))); /////////////
        BBE_SAN_2XF32_IP(Acc, vF, pF);
    }
    l = L % (BBE_SIMD_WIDTH / 4);
    if (l)
    {
        vboolN_2 tmp;

        BBE_LAVN_2XF32_XP(X0, vX, pX, 2 * l*sizeof(float32_t));
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        Acc = BBE_NEGN_2XF32(Acc);
        Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(r0), BBE_MOVNX16_FROMN_2XF32(Acc), BBE_MOVN_FROMN_2(tmp))); /////////////
        BBE_SAVN_2XF32_XP(Acc, vF, pF, 2 * l*sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vF, pF);

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)xt;
    pZr = (const xb_vecN_2xf32 *)((float32_t *)pScr + 2 * L_);
    pZw = (xb_vecN_2xf32 *)((float32_t *)pScr + 2 * L_);
    pZr1 = (const xb_vecN_2xf32 *)((float32_t *)pScr + 3 * L_);
    pZw1 = (xb_vecN_2xf32 *)pScr;
    pF = (xb_vecN_2xf32 *)Fi;
    vZr = BBE_LAN_2XF32_PP(pZr);
    vZw = BBE_ZALIGN();
    vZr1 = BBE_LAN_2XF32_PP(pZr1);
    vF = BBE_LAN_2XF32_PP(pF);
    for (l = L; l > 0; l -= (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LAVN_2XF32_XP(Acc, vZr1, pZr1, BBE_SIMD_WIDTH);
        BBE_LAVN_2XF32_XP(X1, vZr, pZr, BBE_SIMD_WIDTH);
        BBE_LAN_2XF32_IP(F0, vF, pF);
        Acc = BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_DOUBLE_2_LO);
        X1 = BBE_SHFLN_2XF32I(X1, BBE_SHFLI_DOUBLE_2_LO);

        BBE_MULSN_2XF32(X0, X1, F0);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        BBE_MULAN_2XF32(Acc, X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SAVN_2XF32_XP(BBE_SELN_2XF32I(rsqrt, rsqrt, BBE_SELI_EXTRACT_2_OF_4_OFF_0), vZw, pZw, BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pZw1, 2 * BBE_SIMD_WIDTH);
    }
    BBE_SAN_2XF32POS_FP(vZw, pZw);


    pZr = (const xb_vecN_2xf32 *)pScr;
    pV = (long long *)v;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(Acc, vZr, pZr, 2 * sizeof(float32_t));
        BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pV, 2 * M*sizeof(float32_t));
    }

    //pX = (const xb_vecN_2xf32 *)x;
    pV1 = (xb_vecN_2xf32 *)v;
    pZr_ = (const xtfloat *)((float32_t *)pScr + 2 * L_);
    for (l = 0; l < L; l++)
    {
        pX = (const xb_vecN_2xf32 *)(x + 2 * M*l + 2);
        pV1 = (xb_vecN_2xf32 *)XT_ADDI((uintptr_t)pV1, 8);
        vX = BBE_LAN_2XF32_PP(pX);
        vV = BBE_ZALIGN();

        BBE_LSN_2XF32_IP(rsqrt, pZr_, sizeof(float32_t));
        rsqrt = BBE_REPN_2XF32(rsqrt, 0);

        for (m = M - 1; m > 0; m -= (BBE_SIMD_WIDTH / 4))
        {
            BBE_LAN_2XF32_IP(X0, vX, pX);
            Acc = BBE_MULN_2XF32(X0, rsqrt);
            BBE_SAVN_2XF32_XP(Acc, vV, pV1, 2 * m*sizeof(float32_t));
        }
        BBE_SAN_2XF32POS_FP(vV, pV1);
    }
}

#endif

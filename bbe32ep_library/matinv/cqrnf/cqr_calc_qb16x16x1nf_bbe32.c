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
/*------------------------------------------------
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV][2]
    Input/output:
    B[L][SB][2]
------------------------------------------------*/
#if 0
static void cqrnfUpdateB(float32_t* B, const float32_t* V, int M, int L)
{
    int m, l;
    float32_t A_re, A_im;
    const int SB = 32;
    NASSERT(M <= 16);
    for (l = 0; l < L; l++)
    {
        float32_t b_re, b_im, v_re, v_im, z_re, z_im;
        z_re = z_im = 0.f;
        for (m = 0; m < M; m++)
        {
            v_re = V[2 * m + 0]; v_im = V[2 * m + 1];
            b_re = B[2 * m + 0]; b_im = B[2 * m + 1];
            z_re += (v_re*b_re) + (v_im*b_im);
            z_im += (v_re*b_im) - (v_im*b_re);
        }
        z_re *= 2.f;
        z_im *= 2.f;
        for (m = 0; m < M; m++)
        {
            A_re = B[2 * m + 0]; A_im = B[2 * m + 1];
            v_re = V[2 * m + 0]; v_im = V[2 * m + 1];
            A_re -= (z_re*v_re) - (z_im*v_im);
            A_im -= (z_re*v_im) + (z_im*v_re);
            B[2 * m + 0] = A_re;
            B[2 * m + 1] = A_im;
        }
        V += 2 * M;
        B += SB;
    }
}
#endif // 0

static void cqrnfUpdateB16(float32_t* B, const float32_t* V, int M, int L)
{
    int l;
    int count;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV_;
          xb_vecN_2xf32 * restrict pBw;
    
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, V1, V2, V3, B0, B1, B2, B3, Z0;
    valign vV;
    const uint16_t ALIGN(32) p_mask[] = { 
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
        24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7,
        28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vselN mask_;
    vboolN_2 mask_2;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(M > 12 && M <= 16);

    count = (M - 12) * 2 * sizeof(float32_t);
    mask_ = vselN_loadx((const vselN *)p_mask, count * 4);
    mask_2 = BBE_LTRN_2(2 * (16 - M));
    mask_2 = BBE_NOTBN_2(mask_2);

    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    pV_ = (const xb_vecN_2xf32 *)V;
    vV = BBE_LAN_2XF32_PP(pV_);
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(V0), BBE_MOVNX16_FROMN_2XF32(V0), mask_));
        BBE_LVN_2XF32T_IP(B0, pB, 2 * BBE_SIMD_WIDTH, mask_2);
        Acc = BBE_MULMN_2XF32(V0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(V0, B0, 2, 11);

        BBE_LAN_2XF32_IP(V1, vV, pV_);
        BBE_LVN_2XF32_IP(B1, pB, 2 * BBE_SIMD_WIDTH);
        Acc2 = BBE_MULMN_2XF32(V1, B1, 0, 4);
        Acc3 = BBE_MULMN_2XF32(V1, B1, 2, 11);

        BBE_LAN_2XF32_IP(V2, vV, pV_);
        BBE_LVN_2XF32_IP(B2, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, V2, B2, 0, 4);
        BBE_MULMASN_2XF32(Acc1, V2, B2, 2, 11);

        BBE_LAN_2XF32_IP(V3, vV, pV_);
        BBE_LVN_2XF32_XP(B3, pB, -3 * 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, V3, B3, 0, 4);
        BBE_MULMASN_2XF32(Acc3, V3, B3, 2, 11);

        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);

        BBE_LVN_2XF32T_IP(B0, pB, 2 * BBE_SIMD_WIDTH, mask_2);
        BBE_MULMASN_2XF32(B0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(B0, Z0, V0, 2, 11);
        BBE_SVN_2XF32T_IP(B0, pBw, 2 * BBE_SIMD_WIDTH, mask_2);

        BBE_LVN_2XF32_IP(B1, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(B1, Z0, V1, 3, 4);
        BBE_MULMASN_2XF32(B1, Z0, V1, 2, 11);
        BBE_SVN_2XF32_IP(B1, pBw, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_IP(B2, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(B2, Z0, V2, 3, 4);
        BBE_MULMASN_2XF32(B2, Z0, V2, 2, 11);
        BBE_SVN_2XF32_IP(B2, pBw, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_IP(B3, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(B3, Z0, V3, 3, 4);
        BBE_MULMASN_2XF32(B3, Z0, V3, 2, 11);
        BBE_SVN_2XF32_IP(B3, pBw, 2 * BBE_SIMD_WIDTH);
    }
}

static void cqrnfUpdateB12(float32_t* B, const float32_t* V, int M, int L)
{
    int l;
    int count;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV_;
          xb_vecN_2xf32 * restrict pBw;
    
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, V1, V2, B0, B1, B2, Z0;
    valign vV;
    const uint16_t ALIGN(32) p_mask[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
        24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7,
        28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vselN mask_;
    vboolN_2 mask_2;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(M > 8 && M <= 12);

    count = (M - 8) * 2 * sizeof(float32_t);
    mask_ = vselN_loadx((const vselN *)p_mask, count * 4);
    mask_2 = BBE_LTRN_2(2 * (12 - M));
    mask_2 = BBE_NOTBN_2(mask_2);

    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    pV_ = (const xb_vecN_2xf32 *)V;
    vV = BBE_LAN_2XF32_PP(pV_);
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(V0), BBE_MOVNX16_FROMN_2XF32(V0), mask_));
        BBE_LVN_2XF32T_IP(B0, pB, 2 * BBE_SIMD_WIDTH, mask_2);
        Acc = BBE_MULMN_2XF32(V0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(V0, B0, 2, 11);

        BBE_LAN_2XF32_IP(V1, vV, pV_);
        BBE_LVN_2XF32_IP(B1, pB, 2 * BBE_SIMD_WIDTH);
        Acc2 = BBE_MULMN_2XF32(V1, B1, 0, 4);
        Acc3 = BBE_MULMN_2XF32(V1, B1, 2, 11);

        BBE_LAN_2XF32_IP(V2, vV, pV_);
        BBE_LVN_2XF32_IP(B2, pB, 4 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, V2, B2, 0, 4);
        BBE_MULMASN_2XF32(Acc1, V2, B2, 2, 11);

        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);

        BBE_MULMASN_2XF32(B0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(B0, Z0, V0, 2, 11);
        BBE_SVN_2XF32T_IP(B0, pBw, 2 * BBE_SIMD_WIDTH, mask_2);

        BBE_MULMASN_2XF32(B1, Z0, V1, 3, 4);
        BBE_MULMASN_2XF32(B1, Z0, V1, 2, 11);
        BBE_SVN_2XF32_IP(B1, pBw, 2 * BBE_SIMD_WIDTH);

        BBE_MULMASN_2XF32(B2, Z0, V2, 3, 4);
        BBE_MULMASN_2XF32(B2, Z0, V2, 2, 11);
        BBE_SVN_2XF32_IP(B2, pBw, 4 * BBE_SIMD_WIDTH);
    }
}

static void cqrnfUpdateB8(float32_t* B, const float32_t* V, int M, int L)
{
    int l;
    int count;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV_;
          xb_vecN_2xf32 * restrict pBw;
    
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, V1, B0, B1, Z0;
    valign vV;
    const uint16_t ALIGN(32) p_mask[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
        24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7,
        28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vselN mask_;
    vboolN_2 mask_2;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(M > 4 && M <= 8);

    count = (M - 4) * 2 * sizeof(float32_t);
    mask_ = vselN_loadx((const vselN *)p_mask, count * 4);
    mask_2 = BBE_LTRN_2(2 * (8 - M));
    mask_2 = BBE_NOTBN_2(mask_2);

    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    pV_ = (const xb_vecN_2xf32 *)V;
    vV = BBE_LAN_2XF32_PP(pV_);
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(V0), BBE_MOVNX16_FROMN_2XF32(V0), mask_));
        BBE_LVN_2XF32T_IP(B0, pB, 2 * BBE_SIMD_WIDTH, mask_2);
        Acc = BBE_MULMN_2XF32(V0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(V0, B0, 2, 11);

        BBE_LAN_2XF32_IP(V1, vV, pV_);
        BBE_LVN_2XF32_IP(B1, pB, 6 * BBE_SIMD_WIDTH);
        Acc2 = BBE_MULMN_2XF32(V1, B1, 0, 4);
        Acc3 = BBE_MULMN_2XF32(V1, B1, 2, 11);

        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);

        BBE_MULMASN_2XF32(B0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(B0, Z0, V0, 2, 11);
        BBE_SVN_2XF32T_IP(B0, pBw, 2 * BBE_SIMD_WIDTH, mask_2);

        BBE_MULMASN_2XF32(B1, Z0, V1, 3, 4);
        BBE_MULMASN_2XF32(B1, Z0, V1, 2, 11);
        BBE_SVN_2XF32_IP(B1, pBw, 6 * BBE_SIMD_WIDTH);
    }
}

static void cqrnfUpdateB4(float32_t* B, const float32_t* V, int M, int L)
{
    int l;
    int count;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV_;
    xb_vecN_2xf32 * restrict pBw;

    xb_vecN_2xf32 Acc, Acc1, V0, B0, Z0;
    valign vV;
    const uint16_t ALIGN(32) p_mask[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
        24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7,
        28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vselN mask_;
    vboolN_2 mask_2;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(M <= 4);

    count = M * 2 * sizeof(float32_t);
    mask_ = vselN_loadx((const vselN *)p_mask, count * 4);
    mask_2 = BBE_LTRN_2(2 * (4 - M));
    mask_2 = BBE_NOTBN_2(mask_2);

    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    pV_ = (const xb_vecN_2xf32 *)V;
    vV = BBE_LAN_2XF32_PP(pV_);
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(V0), BBE_MOVNX16_FROMN_2XF32(V0), mask_));
        BBE_LVN_2XF32T_IP(B0, pB, 8 * BBE_SIMD_WIDTH, mask_2);
        Acc = BBE_MULMN_2XF32(V0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(V0, B0, 2, 11);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);

        BBE_MULMASN_2XF32(B0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(B0, Z0, V0, 2, 11);
        BBE_SVN_2XF32T_IP(B0, pBw, 8 * BBE_SIMD_WIDTH, mask_2);
    }
}

static void cqrnfRotateB1_16(void *pScr, float32_t* B, const float32_t* Fi, int N, int SB, int L)
{
    int l;

    const long long     * restrict pF;
    const xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pBw;

    xb_vecN_4x64 t0, t1, t2, t3;
    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(N == 16);

    pF = (const long long *)Fi;
    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        BBE_LSN_4X64_XP(t0, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t1, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t2, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t3, pF, 2 * L*sizeof(float32_t));
        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        F0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));

        BBE_LVN_2XF32_XP(B0, pB, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
        BBE_SVN_2XF32_IP(Acc, pBw, 2 * BBE_SIMD_WIDTH);


        BBE_LSN_4X64_XP(t0, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t1, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t2, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t3, pF, 2 * L*sizeof(float32_t));
        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        F0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));

        BBE_LVN_2XF32_XP(B0, pB, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
        BBE_SVN_2XF32_IP(Acc, pBw, 2 * BBE_SIMD_WIDTH);


        BBE_LSN_4X64_XP(t0, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t1, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t2, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t3, pF, 2 * L*sizeof(float32_t));
        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        F0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));

        BBE_LVN_2XF32_XP(B0, pB, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
        BBE_SVN_2XF32_IP(Acc, pBw, 2 * BBE_SIMD_WIDTH);


        BBE_LSN_4X64_XP(t0, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t1, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t2, pF, 2 * L*sizeof(float32_t));
        BBE_LSN_4X64_XP(t3, pF, 2 * sizeof(float32_t) - 15 * 2 * L*sizeof(float32_t));
        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        F0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));

        BBE_LVN_2XF32_XP(B0, pB, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
        BBE_SVN_2XF32_IP(Acc, pBw, 2 * BBE_SIMD_WIDTH);
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
void cqr_calc_qb16x16x1nf(void *pScr,
                    complex_float* _B,const complex_float* _V,int L)
{
#if 0
    float32_t*       B = (float32_t*)_B;
    const float32_t* V = (const float32_t*)_V;
    int m;
    const float32_t* pV;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    if (L <= 0) return;
    for (pV = V, m = 0; m < 16; m++)
    {
        cqrnfUpdateB(B + 2 * m, pV, (16 - m), L);
        pV += 2 * (16 - m)*L;
    }
    cqrnfRotateB1_16(pScr, B, V + (2 * 16 - 16 + 1) * 16 * L, 16, 32, L);
#endif // 0

    int m;
          float32_t * B = (      float32_t *)_B;
    const float32_t * V = (const float32_t *)_V;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    if (L <= 0) return;

    for (m = 16; m > 12; m--)
    {
        cqrnfUpdateB16(B, V, m, L);
        V += 2 * m*L;
    }
    B += 8;
    for (; m > 8; m--)
    {
        cqrnfUpdateB12(B, V, m, L);
        V += 2 * m*L;
    }
    B += 8;
    for (; m > 4; m--)
    {
        cqrnfUpdateB8(B, V, m, L);
        V += 2 * m*L;
    }
    B += 8;
    for (; m > 0; m--)
    {
        cqrnfUpdateB4(B, V, m, L);
        V += 2 * m*L;
    }
    B -= 3 * 8;

    cqrnfRotateB1_16(pScr, B, V, 16, 32, L);
}

size_t cqr_calc_qb16x16x1nf_getScratchSize (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==16 && M==16);
    (void)M,(void)N,(void)P,(void)L;
    return 2*N* (2*BBE_SIMD_WIDTH);
}

#else
DISCARD_FUN(void, cqr_calc_qb16x16x1nf, (void *pScr,
                    complex_float* _B,const complex_float* _V,int L))

size_t cqr_calc_qb16x16x1nf_getScratchSize (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==16 && M==16);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

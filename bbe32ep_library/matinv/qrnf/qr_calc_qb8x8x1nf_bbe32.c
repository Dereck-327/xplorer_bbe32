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
#include "common.h"
#include "qrnf_common.h"

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
------------------------------------------------*/
#if 0
static void qrnfUpdateB(float32_t* B, const float32_t* V, int M, int L)
{
    int m, l;
    const int SB = 8;
    NASSERT(M <= 8);
    for (l = 0; l < L; l++)
    {
        float32_t z_re;
        z_re = 0.f;
        for (m = 0; m < M; m++) z_re += V[m] * B[m];
        z_re *= 2.f;
        for (m = 0; m < M; m++) B[m] -= z_re*V[m];
        V += M;
        B += SB;
    }
}
#endif // 0

static void qrnfUpdateB8(float32_t* B, const float32_t* V, int M, int L)
{
    int l;
    int count;

    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV_;
    xb_vecN_2xf32 * restrict pBw;

    xb_vecN_2xf32 Acc, V0, B0, Z0;
    valign vV;
    const uint16_t ALIGN(32) p_mask[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
        22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5,
        24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7,
        26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
        28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vselN mask_;
    vboolN_2 mask_2;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(M <= 8);

    count = M * sizeof(float32_t);
    mask_ = vselN_loadx((const vselN *)p_mask, count * 8);
    mask_2 = BBE_LTRN_2((8 - M));
    mask_2 = BBE_NOTBN_2(mask_2);

    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    pV_ = (const xb_vecN_2xf32 *)V;
    vV = BBE_LAN_2XF32_PP(pV_);
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(V0, vV, pV_, count);
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(V0), BBE_MOVNX16_FROMN_2XF32(V0), mask_));
        BBE_LVN_2XF32T_XP(B0, pB, 2 * BBE_SIMD_WIDTH, mask_2);
        Acc = BBE_MULN_2XF32(V0, B0);

        BBE_MULMASN_2XF32(Acc, Acc, 1.f, 0, 6);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
        Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);

        BBE_MULSN_2XF32(B0, Z0, V0);
        BBE_SVN_2XF32T_XP(B0, pBw, 2 * BBE_SIMD_WIDTH, mask_2);
    }
}

static void qrnfRotateB1_8(void *pScr, float32_t* B, const float32_t* Fi, int N, int SB, int L)
{
    int l;

    const xtfloat       * restrict pF;
    const xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pBw;

    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;
    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % 4 == 0);

    pF = (const xtfloat *)Fi;
    pB = (const xb_vecN_2xf32 *)B;
    pBw = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        t7 = BBE_LSN_2XF32_X(pF, 7 * L*sizeof(float32_t));
        t6 = BBE_LSN_2XF32_X(pF, 6 * L*sizeof(float32_t));
        t5 = BBE_LSN_2XF32_X(pF, 5 * L*sizeof(float32_t));
        t4 = BBE_LSN_2XF32_X(pF, 4 * L*sizeof(float32_t));
        t3 = BBE_LSN_2XF32_X(pF, 3 * L*sizeof(float32_t));
        t2 = BBE_LSN_2XF32_X(pF, 2 * L*sizeof(float32_t));
        t1 = BBE_LSN_2XF32_X(pF, L*sizeof(float32_t));
        BBE_LSN_2XF32_IP(t0, pF, sizeof(float32_t));
        t0 = BBE_SELN_2XF32I(t1, t0, BBE_SELI_INTERLEAVE_2_LO);
        t2 = BBE_SELN_2XF32I(t3, t2, BBE_SELI_INTERLEAVE_2_LO);
        t4 = BBE_SELN_2XF32I(t5, t4, BBE_SELI_INTERLEAVE_2_LO);
        t6 = BBE_SELN_2XF32I(t7, t6, BBE_SELI_INTERLEAVE_2_LO);
        t0 = BBE_SELN_2XF32I(t2, t0, BBE_SELI_INTERLEAVE_4_LO);
        t4 = BBE_SELN_2XF32I(t6, t4, BBE_SELI_INTERLEAVE_4_LO);
        F0 = BBE_SELN_2XF32I(t4, t0, BBE_SELI_EXTRACT_LO_HALVES);

        BBE_LVN_2XF32_XP(B0, pB, 2 * BBE_SIMD_WIDTH);

        Acc = BBE_MULN_2XF32(F0, B0);

        BBE_SVN_2XF32_XP(Acc, pBw, 2 * BBE_SIMD_WIDTH);
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
void qr_calc_qb8x8x1nf(void *pScr,
                    float32_t* B,const float32_t* V,int L)
{
#if 0
    int m;
    const float32_t* pV;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    if (L <= 0) return;
    for (pV = V, m = 0; m < 8; m++)
    {
        qrnfUpdateB(B + m, pV, (8 - m), L);
        pV += (8 - m)*L;
    }
    qrnfRotateB1_8(pScr, B, V + (((2 * 8 - 8 + 1) * 8) >> 1)*L, 8, 8, L);
#endif // 0

    int m;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    if (L <= 0) return;

    for (m = 8; m > 0; m--)
    {
        qrnfUpdateB8(B, V, m, L);
        V += m*L;
    }

    qrnfRotateB1_8(pScr, B, V, 8, 8, L);
}

size_t qr_calc_qb8x8x1nf_getScratchSize   (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==8 && M==8);
    (void)M,(void)N,(void)P,(void)L;
    return N*(2*BBE_SIMD_WIDTH);
}
#else
DISCARD_FUN(void, qr_calc_qb8x8x1nf, (void *pScr,
                    float32_t* B,const float32_t* V,int L))

size_t qr_calc_qb8x8x1nf_getScratchSize   (int M, int N,int P,int L)
{
    NASSERT(L>0);
    NASSERT(N==8 && M==8);
    (void)M,(void)N,(void)P,(void)L;
    return 0;
}
#endif

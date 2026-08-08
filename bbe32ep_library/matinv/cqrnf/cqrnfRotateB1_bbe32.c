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
#include "common.h"
#include "cqrnf_common.h"

#if HAVE_VFPU

/*------------------------------------------------
    rotate complex B[L][SB] by complex diagonal matrix Fi'[L][SV]
    Input:
    B[L][SB][2]
    Fi'[L][SV][2]
    Output:
    B[L][SB][2]
    Temporary:
    pScr - 2N vectors

    Note, Fi might be non-aligned!
------------------------------------------------*/
void cqrnfRotateB1(void *pScr,float32_t* B,const float32_t* Fi,int N,int SB,int L)
{
#if 0
    int k, n, l;
    float32_t * tmp = (float32_t *)pScr;
    float32_t A_re[BBE_SIMD_WIDTH / 4], A_im[BBE_SIMD_WIDTH / 4];
    float32_t f_re[BBE_SIMD_WIDTH / 4], f_im[BBE_SIMD_WIDTH / 4], b_re[BBE_SIMD_WIDTH / 4], b_im[BBE_SIMD_WIDTH / 4];
    int L0, Lmod;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);

    L0 = L & ~(BBE_SIMD_WIDTH / 4 - 1);
    Lmod = L - L0;
    if (L0)
    {
        for (n = 0; n < N; n += BBE_SIMD_WIDTH / 4)
        {
            for (l = 0; l < L0; l++)
            {
                for (k = 0; k < BBE_SIMD_WIDTH / 4; k++)
                {
                    f_re[k] = Fi[2 * l + 2 * (n + k)*L + 0];
                    f_im[k] = Fi[2 * l + 2 * (n + k)*L + 1];
                    b_re[k] = B[2 * (n + k) + SB*l + 0];
                    b_im[k] = B[2 * (n + k) + SB*l + 1];
                    A_re[k] = (f_re[k] * b_re[k]) + (f_im[k] * b_im[k]);
                    A_im[k] = (f_re[k] * b_im[k]) - (f_im[k] * b_re[k]);
                    B[2 * (n + k) + SB*l + 0] = A_re[k];
                    B[2 * (n + k) + SB*l + 1] = A_im[k];
                }
            }
        }
    }
    if (Lmod == 0) return;
    // process last portion, first load with transpose, next 
    // process as usual before in SIMD manner
    for (n = 0; n < N; n += BBE_SIMD_WIDTH / 4)
    {
        for (l = L0; l < L; l++)
        {
            for (k = 0; k < BBE_SIMD_WIDTH / 4; k++)
            {
                tmp[2 * (l - L0)*N + 2 * (n + k) + 0] = Fi[2 * l + 2 * (n + k)*L + 0];
                tmp[2 * (l - L0)*N + 2 * (n + k) + 1] = Fi[2 * l + 2 * (n + k)*L + 1];
            }
        }
    }

    for (n = 0; n < N; n += BBE_SIMD_WIDTH / 4)
    {
        for (l = L0; l < L; l++)
        {
            for (k = 0; k < BBE_SIMD_WIDTH / 4; k++)
            {
                f_re[k] = tmp[2 * (l - L0)*N + 2 * (n + k) + 0];
                f_im[k] = tmp[2 * (l - L0)*N + 2 * (n + k) + 1];
                b_re[k] = B[2 * (n + k) + SB*l + 0];
                b_im[k] = B[2 * (n + k) + SB*l + 1];
                A_re[k] = (f_re[k] * b_re[k]) + (f_im[k] * b_im[k]);
                A_im[k] = (f_re[k] * b_im[k]) - (f_im[k] * b_re[k]);
                B[2 * (n + k) + SB*l + 0] = A_re[k];
                B[2 * (n + k) + SB*l + 1] = A_im[k];
            }
        }
    }
#endif // 0

    int n, l;

    const long long     * restrict pF;
    const long long     * restrict pF_;
    const xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pBw;
          xb_vecN_2xf32 * restrict pB_;

    xb_vecN_4x64 t0, t1, t2, t3;
    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % 4 == 0);

    pF_ = (const long long *)Fi;
    pB_ = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        pB = pB_;
        pBw = pB_;
        pF = pF_;
        for (n = 0; n < N; n += BBE_SIMD_WIDTH / 4)
        {
            BBE_LSN_4X64_XP(t0, pF, 2 * L*sizeof(float32_t));
            BBE_LSN_4X64_XP(t1, pF, 2 * L*sizeof(float32_t));
            BBE_LSN_4X64_XP(t2, pF, 2 * L*sizeof(float32_t));
            BBE_LSN_4X64_XP(t3, pF, 2 * L*sizeof(float32_t));
            t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
            t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
            F0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));
                
            BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);

            Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);

            BBE_SVN_2XF32_IP(Acc, pBw, 2 * BBE_SIMD_WIDTH);
        }
        pF_ = (const long long *)XT_ADDX4(2, (uintptr_t)pF_);
        pB_ = (xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}

#endif

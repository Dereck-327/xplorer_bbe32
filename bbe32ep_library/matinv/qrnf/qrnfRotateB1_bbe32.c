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
#include "common.h"
#include "qrnf_common.h"

#if HAVE_VFPU

/*------------------------------------------------
    rotate real B[L][SB] by real diagonal matrix Fi'[L][SV]
    Input:
    B[L][SB]
    Fi'[L][SV]
    Output:
    B[L][SB]
    Temporary:
    pScr - N vectors

    Note, Fi might be non-aligned!
------------------------------------------------*/
void qrnfRotateB1(void *pScr,float32_t* B,const float32_t* Fi,int N,int SB,int L)
{
#if 0
    int k, n, l;
    float32_t * tmp = (float32_t *)pScr;
    float32_t A_re[BBE_SIMD_WIDTH / 2], f_re[BBE_SIMD_WIDTH / 2], b_re[BBE_SIMD_WIDTH / 2];
    int L0, Lmod;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);

    L0 = L & ~(BBE_SIMD_WIDTH / 2 - 1);
    Lmod = L - L0;
    if (L0)
    {
        for (n = 0; n < N; n += BBE_SIMD_WIDTH / 2)
        {
            int K = XT_MIN(N - n, BBE_SIMD_WIDTH / 2);
            for (l = 0; l < L0; l++)
            {
                for (k = 0; k < K; k++)
                {
                    f_re[k] = Fi[l + (n + k)*L + 0];
                    b_re[k] = B[(n + k) + SB*l + 0];
                    A_re[k] = f_re[k] * b_re[k];
                    B[(n + k) + SB*l + 0] = A_re[k];
                }
            }
        }
    }
    if (Lmod == 0) return;
    // process last portion, first load with transpose, next 
    // process as usual before in SIMD manner
    for (n = 0; n < N; n += BBE_SIMD_WIDTH / 2)
    {
        int K = XT_MIN(N - n, BBE_SIMD_WIDTH / 2);
        for (l = L0; l < L; l++)
        {
            for (k = 0; k < K; k++)
            {
                tmp[(l - L0)*N + (n + k)] = Fi[l + (n + k)*L];
            }
        }
    }

    for (n = 0; n < N; n += BBE_SIMD_WIDTH / 2)
    {
        int K = XT_MIN(N - n, BBE_SIMD_WIDTH / 2);
        for (l = L0; l < L; l++)
        {
            for (k = 0; k < K; k++)
            {
                f_re[k] = tmp[(l - L0)*N + (n + k) + 0];
                b_re[k] = B[(n + k) + SB*l + 0];
                A_re[k] = f_re[k] * b_re[k];
                B[(n + k) + SB*l + 0] = A_re[k];
            }
        }
    }
#endif // 0

    int n, l;

    const xtfloat       * restrict pF;
    const xtfloat       * restrict pF_;
    const xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pBw;
          xb_vecN_2xf32 * restrict pB_;

    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;
    xb_vecN_2xf32 Acc, F0, B0;
    valign vB, vBw;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % 4 == 0);

    pF_ = (const xtfloat *)Fi;
    pB_ = (xb_vecN_2xf32 *)B;
    for (l = 0; l < L; l++)
    {
        pB = pB_;
        pBw = pB_;
        pF = pF_;
        vB = BBE_LAN_2XF32_PP(pB);
        vBw = BBE_ZALIGN();
        for (n = 0; n < (N&(~(BBE_SIMD_WIDTH / 2 - 1))); n += (BBE_SIMD_WIDTH / 2))
        {
            BBE_LSN_2XF32_XP(t0, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t1, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t2, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t3, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t4, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t5, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t6, pF, L*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t7, pF, L*sizeof(float32_t));
            t0 = BBE_SELN_2XF32I(t1, t0, BBE_SELI_INTERLEAVE_2_LO);
            t2 = BBE_SELN_2XF32I(t3, t2, BBE_SELI_INTERLEAVE_2_LO);
            t4 = BBE_SELN_2XF32I(t5, t4, BBE_SELI_INTERLEAVE_2_LO);
            t6 = BBE_SELN_2XF32I(t7, t6, BBE_SELI_INTERLEAVE_2_LO);
            t0 = BBE_SELN_2XF32I(t2, t0, BBE_SELI_INTERLEAVE_4_LO);
            t4 = BBE_SELN_2XF32I(t6, t4, BBE_SELI_INTERLEAVE_4_LO);
            F0 = BBE_SELN_2XF32I(t4, t0, BBE_SELI_EXTRACT_LO_HALVES);
                
            BBE_LAN_2XF32_IP(B0, vB, pB);

            Acc = BBE_MULN_2XF32(F0, B0);

            BBE_SAN_2XF32_IP(Acc, vBw, pBw);
        }
        for (; n < N; n++)
        {
            BBE_LSN_2XF32_XP(t0, pF, L*sizeof(float32_t));

            BBE_LAVN_2XF32_XP(B0, vB, pB, sizeof(float32_t));

            Acc = BBE_MULN_2XF32(t0, B0);

            BBE_SAVN_2XF32_XP(Acc, vBw, pBw, sizeof(float32_t));
        }
        BBE_SAN_2XF32POS_FP(vBw, pBw);

        pF_ = (const xtfloat *)XT_ADDX4(1, (uintptr_t)pF_);
        pB_ = (xb_vecN_2xf32 *)XT_ADDX4(SB, (uintptr_t)pB_);
    }
}

#endif

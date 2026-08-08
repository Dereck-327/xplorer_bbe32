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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    Cholesky decomposition, floating point complex data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cholnf_common.h"

#if (HAVE_VFPU)

/*---------------------------------------------------
   compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   sigma2[L]    regularization term
   n            number of column
   Output:
   Z[L][n+1][2] results
---------------------------------------------------*/
void cholnfConv8_11(float32_t* Z,const float32_t* A,const float32_t * restrict sigma2,int n,int M,int N,int L,int SA)
{
    int l, m;
          xb_vecN_2xf32 * restrict pZ;
    const xb_vecN_2xf32 * restrict pAk;
    const xb_vecN_2xf32 * restrict pA;
    const xtfloat       * restrict pAn;
    const xtfloat       * restrict pS;

    valign vz;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, Acc4, Acc5, sigma;
    xb_vecN_2xf32 A0, A1;

    NASSERT(n >= 8 && n <= 11);
    NASSERT(N % 4 == 0 && M % 4 == 0 && SA == 2 * (M*N));
    NASSERT_ALIGN(A, 2 * BBE_SIMD_WIDTH);

    pS = (const xtfloat *)sigma2;
    vz = BBE_ZALIGN();
    pZ = (xb_vecN_2xf32 *)(Z);
    pA = (const xb_vecN_2xf32 *)(A);

    for (l = 0; l < L; l++)
    {
        BBE_LSN_2XF32_IP(sigma, pS, 4);
        Acc = BBE_SHFLN_2XF32I(sigma, BBE_SHFLI_REP_0X4);
        Acc1 = BBE_ZERON_2XF32();
        Acc2 = Acc;
        Acc3 = BBE_ZERON_2XF32();
        Acc4 = Acc;
        Acc5 = BBE_ZERON_2XF32();
        pAn = (const xtfloat *)XT_ADDX8(n, (uintptr_t)pA);
        pAk = (const xb_vecN_2xf32 *)(pA);

        __Pragma("loop_count min=4");
        for (m = 0; m < M; m++)
        {
            A1 = BBE_LSN_2XF32_I(pAn, 4);
            BBE_LSN_2XF32_XP(A0, pAn, 2 * N * sizeof(float32_t));
            A1 = BBE_SELN_2XF32I(A1, A0, BBE_SELI_PACK_2);
            A1 = BBE_SHFLN_2XF32I(A1, BBE_SHFLI_REP_0X4);
            A0 = BBE_LVN_2XF32_I(pAk, 4 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc1, A0, A1, 2, 11);

            A0 = BBE_LVN_2XF32_I(pAk, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc2, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);

            BBE_LVN_2XF32_XP(A0, pAk, 2 * N * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc4, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc5, A0, A1, 2, 11);
        }
        Acc4 = BBE_ADDN_2XF32(Acc4, Acc5);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        BBE_SAVN_2XF32_XP(Acc4, vz, pZ, 2 * (n + 1) * sizeof(float32_t));
        BBE_SAVN_2XF32_XP(Acc2, vz, pZ, 2 * (n + 1 - (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
        BBE_SAVN_2XF32_XP(Acc, vz, pZ, 2 * (n + 1 - 2 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));

        pA = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pA);
    }
    BBE_SAN_2XF32POS_FP(vz, pZ);
}
#endif

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
    take column from sequence of block ordered matrices
    and put it to the linear array
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
void cqrnfTakeColumn(float32_t* restrict x, float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L)
{
#if 0
    int l, m;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < L; l++)
        for (m = 0; m < M; m++)
        {
            xt[m * 2 * L + 2 * l + 0] = x[l * 2 * M + 2 * m + 0] = A[l*SA + m * 2 * N + 0];
            xt[m * 2 * L + 2 * l + 1] = x[l * 2 * M + 2 * m + 1] = A[l*SA + m * 2 * N + 1];
        }
#endif // 0

    int l, m;

    const long long  * restrict pA  = (const long long  *)A;
    const long long  * restrict pA0 = (const long long  *)A;
          long long  * restrict pXt = (      long long  *)xt;
          long long  * restrict pXt0= (      long long  *)xt;
          xb_vecNx16 * restrict pX  = (      xb_vecNx16 *)x;

    xb_vecN_4x64 t0, t1, t2, t3;
    xb_vecNx16 vTmp;
    valign vX;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);

    for (l = 0; l < L; l++)
    {
        pA = pA0;
        pXt = pXt0;

        for (m = 0; m < (M&(~(BBE_SIMD_WIDTH / 4 - 1))); m += (BBE_SIMD_WIDTH / 4))
        {
            BBE_LSN_4X64_XP(t0, pA, 2 * N * sizeof(float32_t));
            BBE_SSN_4X64_XP(t0, pXt, 2 * L * sizeof(float32_t));
            BBE_LSN_4X64_XP(t1, pA, 2 * N * sizeof(float32_t));
            BBE_SSN_4X64_XP(t1, pXt, 2 * L * sizeof(float32_t));
            BBE_LSN_4X64_XP(t2, pA, 2 * N * sizeof(float32_t));
            BBE_SSN_4X64_XP(t2, pXt, 2 * L * sizeof(float32_t));
            BBE_LSN_4X64_XP(t3, pA, 2 * N * sizeof(float32_t));
            BBE_SSN_4X64_XP(t3, pXt, 2 * L * sizeof(float32_t));

            t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
            t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
            vTmp = BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES);
            BBE_SANX16_IP(vTmp, vX, pX);
        }
        for (; m < M; m++)
        {
            BBE_LSN_4X64_XP(t0, pA, 2 * N * sizeof(float32_t));
            BBE_SSN_4X64_XP(t0, pXt, 2 * L * sizeof(float32_t));

            vTmp = BBE_MOVNX16_FROMN_4X64(t0);
            BBE_SAVNX16_XP(vTmp, vX, pX, 2 * sizeof(float32_t));
        }

        pXt0 = (long long *)XT_ADDX4(2, (uintptr_t)pXt0);
        pA0 = (const long long *)XT_ADDX4(SA, (uintptr_t)pA0);
    }
    BBE_SANX16POS_FP(vX, pX);
}
#endif

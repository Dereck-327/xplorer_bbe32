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
// N=13...16
void cqrnfUpdateR16(          float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L)
{
    int m, l;

    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pR_;
    const long long     * restrict pV;
          xb_vecN_2xf32 * restrict pRw;
          xb_vecN_2xf32 * restrict pZ = (xb_vecN_2xf32 *)Z;

    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, R0, Z0, Z1, Z2, Z3;

    int K = ((N + (BBE_SIMD_WIDTH / 4 - 1)) &~(BBE_SIMD_WIDTH / 4 - 1)) - N;
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N0 % 4 == 0);
    NASSERT(N >= 13 && N <= 16);
    /* this is trick: we may simply update bit wider matrix because
       we always updating upper triangle part of original matrix, so
       we have zeroes in K columns from the left side !
       */
    R -= 2 * K;
    N += K;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);

    pZ = (xb_vecN_2xf32 *)Z;
    pR_ = (const xb_vecN_2xf32 *)R;
    pV = (const long long *)v;
    for (l = 0; l < L; l++)
    {
        pR = pR_;

        Acc = BBE_ZERON_2XF32();
        Acc1 = BBE_ZERON_2XF32();
        Acc2 = BBE_ZERON_2XF32();
        Acc3 = BBE_ZERON_2XF32();
        for (m = 0; m < M; m++)
        {
            BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
            V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, V0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc, V0, R0, 2, 11);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc1, V0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc1, V0, R0, 2, 11);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc2, V0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc2, V0, R0, 2, 11);

            BBE_LVN_2XF32_XP(R0, pR, 2 * N0*sizeof(float32_t) - 3 * 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc3, V0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc3, V0, R0, 2, 11);
        }
        Z0 = BBE_MULN_2XF32(Acc, 2.f);
        Z1 = BBE_MULN_2XF32(Acc1, 2.f);
        Z2 = BBE_MULN_2XF32(Acc2, 2.f);
        Z3 = BBE_MULN_2XF32(Acc3, 2.f);
        BBE_SVN_2XF32_IP(Z0, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Z1, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Z2, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Z3, pZ, 0);


        pV = (const long long *)XT_ADDX4(-2 * M, (uintptr_t)pV);
        pR = pR_;
        pRw = (xb_vecN_2xf32 *)pR_;
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_IP(Z3, pZ, -2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(Z2, pZ, -2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(Z1, pZ, -2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(Z0, pZ, 3 * 2 * BBE_SIMD_WIDTH);

            BBE_LSN_4X64_XP(vTmp, pV, 2 * sizeof(float32_t));
            V0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            BBE_LVN_2XF32_XP(Acc, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z0, V0, 2, 11);
            BBE_SVN_2XF32_XP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(Acc, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, Z1, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z1, V0, 2, 11);
            BBE_SVN_2XF32_XP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(Acc, pR, 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, Z2, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z2, V0, 2, 11);
            BBE_SVN_2XF32_XP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(Acc, pR, 2 * N0*sizeof(float32_t) - 3 * 2 * BBE_SIMD_WIDTH);
            BBE_MULMASN_2XF32(Acc, Z3, V0, 3, 4);
            BBE_MULMASN_2XF32(Acc, Z3, V0, 2, 11);
            BBE_SVN_2XF32_XP(Acc, pRw, 2 * N0*sizeof(float32_t) - 3 * 2 * BBE_SIMD_WIDTH);
        }

        pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
    }
}
#endif

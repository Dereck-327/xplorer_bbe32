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
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "qrnf_common.h"

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
    Z[N*L]
-------------------------------------------------------*/
// N=5..8
void qrnfUpdateR8(float32_t* restrict Z,
                float32_t* restrict R,
                const float32_t* restrict v,
                int SA, int M,int N, int N0, int L)
{
    int m, l;

    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pR0;
    const xb_vecN_2xf32 * restrict pR1;
    const xb_vecN_2xf32 * restrict pR_;
    const xtfloat       * restrict pV;
          xb_vecN_2xf32 * restrict pRw;

    xb_vecN_2xf32 Acc, Acc2, V0, R0, Z0;
    valign vR, vR1;

    int K = ((N + (BBE_SIMD_WIDTH / 2 - 1)) &~(BBE_SIMD_WIDTH / 2 - 1)) - N;
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N0 % 4 == 0);
    NASSERT(N >= 5 && N <= 8);
    /* this is trick: we may simply update bit wider matrix because
       we always updating upper triangle part of original matrix, so
       we have zeroes in K columns from the left side !
       */
    R -= K;
    N += K;
    NASSERT_ALIGN(R, sizeof(float32_t) * 4);

    pR_ = (const xb_vecN_2xf32 *)R;
    pV = (const xtfloat *)v;
    for (l = 0; l < L; l++)
    {
        pR = pR_;
        pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(N0, (uintptr_t)pR_);
        Acc = BBE_ZERON_2XF32();
        Acc2 = BBE_ZERON_2XF32();
        for (m = 0; m < M >> 1; m++)
        {
            BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
            V0 = BBE_REPN_2XF32(V0, 0);
            vR = BBE_LAN_2XF32_PP(pR);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            BBE_MULAN_2XF32(Acc, V0, R0);

            pR = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0 - 8, (uintptr_t)pR);

            BBE_LSN_2XF32_XP(V0, pV, sizeof(float32_t));
            V0 = BBE_REPN_2XF32(V0, 0);
            vR1 = BBE_LAN_2XF32_PP(pR1);
            BBE_LAN_2XF32_IP(R0, vR1, pR1);
            BBE_MULAN_2XF32(Acc2, V0, R0);

            pR1 = (const xb_vecN_2xf32 *)XT_ADDX4(2 * N0 - 8, (uintptr_t)pR1);
        }
        if (M & 1)
        {
            BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
            V0 = BBE_REPN_2XF32(V0, 0);
            vR = BBE_LAN_2XF32_PP(pR);
            BBE_LAN_2XF32_IP(R0, vR, pR);
            BBE_MULAN_2XF32(Acc, V0, R0);
        }
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Z0 = BBE_MULN_2XF32(Acc, 2.f);

        pV = (const xtfloat *)XT_ADDX4(-M, (uintptr_t)pV);
        pR0 = pR_;
        for (m = 0; m < M; m++)
        {
            BBE_LSN_2XF32_IP(V0, pV, sizeof(float32_t));
            V0 = BBE_REPN_2XF32(V0, 0);

            pR = pR0;
            pRw = (xb_vecN_2xf32 *)pR0;
            vR = BBE_LAN_2XF32_PP(pR);
            vR1 = BBE_ZALIGN();

            BBE_LAN_2XF32_IP(Acc, vR, pR);
            BBE_MULSN_2XF32(Acc, Z0, V0);
            BBE_SAN_2XF32_IP(Acc, vR1, pRw);

            BBE_SAN_2XF32POS_FP(vR1, pRw);

            pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(N0, (uintptr_t)pR0);
        }

        pR_ = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
    }
}
#endif

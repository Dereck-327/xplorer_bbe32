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

/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
void cholnfFwdrec16(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
{
    int l, n;
          xb_vecN_2xf32 * restrict pY;
          long long     * restrict pYw;
    const xb_vecN_2xf32 * restrict pR;
    const long long     * restrict pZ;
    const long long     * restrict pD;

    valign vy, vr;
    xb_vecN_2xf32 Acc, Acc1, Y0, R0, D0;
    xb_vecN_4x64 temp;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N <= 16);

    cholnfFwdrec12(y, R, D, Z, 12, L, SR, SD, SY, SZ);

    for (n = 12; n < N; n++)
    {
        pY = (xb_vecN_2xf32 *)(y);
        pYw = (long long *)(y + 2 * n);
        pR = (const xb_vecN_2xf32 *)(R + (n*(n + 1)));
        pZ = (const long long *)(Z + 2 * n);
        pD = (const long long *)(D + 2 * n);

        for (l = 0; l < L; l++)
        {
            vy = BBE_LAN_2XF32_PP(pY);
            vr = BBE_LAN_2XF32_PP(pR);
            
            BBE_LSN_4X64_XP(temp, pZ, SZ * sizeof(float32_t));
            Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));

            BBE_LAN_2XF32_IP(R0, vr, pR);
            BBE_LAN_2XF32_IP(Y0, vy, pY);
            Acc1 = BBE_MULMN_2XF32(R0, Y0, 3, 4);
            BBE_MULMASN_2XF32(Acc, R0, Y0, 1, 11);

            BBE_LAN_2XF32_IP(R0, vr, pR);
            BBE_LAN_2XF32_IP(Y0, vy, pY);
            BBE_MULMASN_2XF32(Acc1, R0, Y0, 3, 4);
            BBE_MULMASN_2XF32(Acc, R0, Y0, 1, 11);

            BBE_LAN_2XF32_IP(R0, vr, pR);
            BBE_LAN_2XF32_IP(Y0, vy, pY);
            BBE_MULMASN_2XF32(Acc1, R0, Y0, 3, 4);
            BBE_MULMASN_2XF32(Acc, R0, Y0, 1, 11);

            BBE_LAVN_2XF32_XP(R0, vr, pR, 2 * (n - 3 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(Y0, vy, pY, 2 * (n - 3 * (BBE_SIMD_WIDTH / 4)) * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc1, R0, Y0, 3, 4);
            BBE_MULMASN_2XF32(Acc, R0, Y0, 1, 11);

            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
            BBE_LSN_4X64_XP(temp, pD, SD * sizeof(float32_t));
            D0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(temp));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pYw, SY * sizeof(float32_t));

            pY = (xb_vecN_2xf32 *)XT_ADDX4(SY - 2 * n, (uintptr_t)pY);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - 2 * n, (uintptr_t)pR);
        }
    }
}
#endif

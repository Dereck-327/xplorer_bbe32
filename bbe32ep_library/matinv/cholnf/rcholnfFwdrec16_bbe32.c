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
    Cholesky decomposition, floating point real data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "rcholnf_common.h"

#if (HAVE_VFPU)
/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
void rcholnfFwdrec16(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
{
    int l, n;
          xb_vecN_2xf32 * restrict pY;
          xtfloat       * restrict pYw;
    const xb_vecN_2xf32 * restrict pR;
    const xtfloat       * restrict pZ;
    const xtfloat       * restrict pD;

    valign vy, vr;
    xb_vecN_2xf32 Acc, Y0, R0, D0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);
    NASSERT(N <= 16);

    for (n = 0; n < N; n++)
    {
        pY = (xb_vecN_2xf32 *)(y);
        pYw = (xtfloat *)(y + n);
        pR = (const xb_vecN_2xf32 *)(R + ((n*(n + 1)) >> 1));
        pZ = (const xtfloat *)(Z + n);
        pD = (const xtfloat *)(D + n);

        for (l = 0; l < L; l++)
        {
            vy = BBE_LAN_2XF32_PP(pY);
            vr = BBE_LAN_2XF32_PP(pR);
            
            BBE_LSN_2XF32_XP(Acc, pZ, SZ * sizeof(float32_t));

            BBE_LAVN_2XF32_XP(R0, vr, pR, n * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(Y0, vy, pY, n * sizeof(float32_t));
            BBE_MULSN_2XF32(Acc, R0, Y0);
            BBE_LAVN_2XF32_XP(R0, vr, pR, (n - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
            BBE_LAVN_2XF32_XP(Y0, vy, pY, (n - (BBE_SIMD_WIDTH / 2)) * sizeof(float32_t));
            BBE_MULSN_2XF32(Acc, R0, Y0);

            BBE_MULMASN_2XF32(Acc, Acc, 1.f, 0, 6);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_XP(Acc, pYw, SY * sizeof(float32_t));

            pY = (xb_vecN_2xf32 *)XT_ADDX4(SY - n, (uintptr_t)pY);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - n, (uintptr_t)pR);
        }
    }
}
#endif

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
    Cholesky decomposition, floating point real data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cholnf_common.h"

#if HAVE_VFPU

#define VECLEN (BBE_SIMD_WIDTH/2)

/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
void rcholnfFwdrec(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
{
#if 0
    int l;
    int n, m;
    const float32_t *pR;
    //NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);

    for (n = 0; n < N; n++)
    {
        for (l = 0; l < L; l++)
        {
            float32_t B_re;
            pR = R + ((n*(n + 1)) >> 1) + l*SR;
            // calculate A(:,n)'*B-Rn'*Y, 1xP
            B_re = Z[SZ*l + n];
            for (m = 0; m < n; m++)
            {
                float32_t r_re;
                float32_t y_re;
                r_re = pR[m];
                y_re = y[l*SY + m];
                B_re -= (y_re*r_re);
            }
            y[l*SY + m] = B_re*D[l*SD + n];
        }
    }
#endif // 0

    int l, n, m;
          xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pR;
    const xtfloat       * restrict pZ;
    const xtfloat       * restrict pD;

    valign vy, vr;
    xb_vecN_2xf32 Acc, Y0, R0, D0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    for (n = 0; n < N; n++)
    {
        pY = (xb_vecN_2xf32 *)(y);
        pR = (const xb_vecN_2xf32 *)(R + ((n*(n + 1)) >> 1));
        pZ = (const xtfloat *)(Z + n);
        pD = (const xtfloat *)(D + n);

        for (l = 0; l < L; l++)
        {
            vy = BBE_LAN_2XF32_PP(pY);
            vr = BBE_LAN_2XF32_PP(pR);
            
            BBE_LSN_2XF32_XP(Acc, pZ, SZ * sizeof(float32_t));

            for (m = n; m > 0; m -= VECLEN)
            {
                BBE_LAVN_2XF32_XP(R0, vr, pR, m * sizeof(float32_t));
                BBE_LAVN_2XF32_XP(Y0, vy, pY, m * sizeof(float32_t));
                BBE_MULSN_2XF32(Acc, R0, Y0);
            }
            BBE_MULMASN_2XF32(Acc, Acc, 1.f, 0, 6);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_4), Acc);
            Acc = BBE_ADDN_2XF32(BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_8), Acc);
            BBE_LSN_2XF32_XP(D0, pD, SD * sizeof(float32_t));
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SSN_2XF32_I(Acc, (xtfloat *)pY, 0);

            pY = (xb_vecN_2xf32 *)XT_ADDX4(SY - n, (uintptr_t)pY);
            pR = (const xb_vecN_2xf32 *)XT_ADDX4(SR - n, (uintptr_t)pR);
        }
    }
}
#endif

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

#define VECLEN (BBE_SIMD_WIDTH/4)

/*-------------------------------------------------
   update N-th diagonal element
   Input:
   Z[L][N+1][2]  convolutions in N-th column
   Input/output:
   y             pointer to the begining of column 
                 in matrix R[L][SR] (N+1 elements 
                 is written)
   Output:
   D[L][SD]      reciprocals of main diagonal 
                 (pointer to the N-th element
-------------------------------------------------*/
void cholnfDiagUpd(float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD)
{
#if 0
    float32_t B_re;
    int l, m;
    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    Z += 2 * N; // points to diagonal element
    for (l = 0; l < L; l++)
    {
        B_re = Z[2 * l*(N + 1) + 0];
        for (m = 0; m < N; m++)
        {
            float32_t r_re, r_im;
            r_re = y[2 * m + 0]; r_im = y[2 * m + 1];
            B_re -= (r_re*r_re) + (r_im*r_im);
        }
        D[0] = D[1] = 1.f / sqrtf(B_re);
        y[2 * N + 0] = D[0] * B_re;
        y[2 * N + 1] = 0.f;
        y += SR;
        D += SD;
    }
#endif // 0

    int l, m;
    const xtfloat       * restrict pZ;
          xtfloat       * restrict pD;
    const xb_vecN_2xf32 * restrict pY;

    valign vy;
    xb_vecN_2xf32 Acc, Y0, rsqrt;
    xtfloat out;

    NASSERT_ALIGN(Z, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    pZ = (const xtfloat *)(Z + 2 * N);
    pD = (xtfloat *)D;

    for (l = 0; l < L; l++)
    {
        BBE_LSN_2XF32_XP(Acc, pZ, 2 * (N + 1) * sizeof(float32_t));
        pY = (const xb_vecN_2xf32 *)(y);
        vy = BBE_LAN_2XF32_PP(pY);

        for (m = N; m > 0; m -= VECLEN)
        {
            BBE_LAVN_2XF32_XP(Y0, vy, pY, 2 * m * sizeof(float32_t));
            BBE_MULSN_2XF32(Acc, Y0, Y0);
        }
        out = BBE_RADDN_2XF32(Acc);
        Acc = BBE_MOVN_2XF32_FROMF32(out);
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);

        BBE_SSN_2XF32_I(rsqrt, pD, 4);
        BBE_SSN_2XF32_XP(rsqrt, pD, SD * sizeof(float32_t));
        //BBE_SSN_2XF32_I(0.f, (xtfloat *)pY, 4);
        BBE_SSN_2XF32_I(Acc, (xtfloat *)pY, 0);

        y += SR;
    }
}

#endif

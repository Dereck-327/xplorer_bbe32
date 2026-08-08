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
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "cqrsf_common.h"
#include <float.h>
#if HAVE_VFPU

/*-----------------------------------------------------------
    find Householder vectors (V and Fi), diagonal element D
    Input:
    x[]       pointer to the diagonal element of original matrix A given 
              in streaming order
    NL,SV,SD  strides for x, V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Householder vectors (M elements each)
-----------------------------------------------------------*/
// M==3
void cqrfsHousholder3(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX1;
    const xb_vecN_2xf32 * restrict pX2;
    const xb_vecN_2xf32 * restrict pZr;
          xb_vecN_2xf32 * restrict pZ;
          xb_vecN_2xf32 * restrict pD;
          xb_vecN_2xf32 * restrict pF;
          xb_vecN_2xf32 * restrict pV;

    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, X0, X1, F0, rsqrt;
    xb_vecN_2xf32 r0;
    const float32_t ALIGN(32) rep_tbl[8] = { 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f };

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);
    NASSERT(M == 3);

    // compute d=1/sqrt(x'x)
    pX = (const xb_vecN_2xf32 *)x;
    pX1 = (const xb_vecN_2xf32 *)(x + 2 * N*L);
    pX2 = (const xb_vecN_2xf32 *)(x + 4 * N*L);
    pZ = (xb_vecN_2xf32 *)pScr;
    pD = (xb_vecN_2xf32 *)D;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        Acc = FLT_MIN;
        Acc1 = FLT_MIN;

            BBE_LVN_2XF32_XP(X0, pX2, 2 * BBE_SIMD_WIDTH);
            X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
            Acc2 = BBE_MULN_2XF32(X0, X0);
            BBE_MULAN_2XF32(Acc, X1, X1);
            Acc3 = BBE_MULN_2XF32(X0, X0);
            BBE_MULAN_2XF32(Acc1, X1, X1);
            BBE_LVN_2XF32_XP(X0, pX1, 2 * BBE_SIMD_WIDTH);
            X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
            BBE_MULAN_2XF32(Acc2, X0, X0);
            BBE_MULAN_2XF32(Acc, X1, X1);
            BBE_MULAN_2XF32(Acc3, X0, X0);
            BBE_MULAN_2XF32(Acc1, X1, X1);
        BBE_LVN_2XF32_XP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        BBE_MULAN_2XF32(Acc2, X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);

        Acc1 = BBE_ADDN_2XF32(Acc1, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X(Acc1, pZ, 2 * L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
    }

    // compute fi=x0/sqrt(x0'x0)
    pX = (const xb_vecN_2xf32 *)x;
    pF = (xb_vecN_2xf32 *)Fi;
    r0 = BBE_LVN_2XF32_I((const xb_vecN_2xf32 *)rep_tbl, 0); /////////////
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        vboolN_2 tmp;

        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        Acc = BBE_NEGN_2XF32(Acc);
        Acc = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(r0), BBE_MOVNX16_FROMN_2XF32(Acc), BBE_MOVN_FROMN_2 (tmp))); /////////////
        BBE_SVN_2XF32_IP(Acc, pF, 2 * BBE_SIMD_WIDTH);
    }

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)x;
    pZr = (const xb_vecN_2xf32 *)pScr;
    pZ = (xb_vecN_2xf32 *)pScr;
    pF = (xb_vecN_2xf32 *)Fi;
    pV = (xb_vecN_2xf32 *)v;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_LVN_2XF32_X(pZr, 2 * L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(X1, pZr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        BBE_MULSN_2XF32(X0, X1, F0);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        BBE_MULAN_2XF32(Acc, X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SVN_2XF32_IP(rsqrt, pZ, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);
    }

    pX = (const xb_vecN_2xf32 *)(x + 2 * N*L);
    pZr = (const xb_vecN_2xf32 *)pScr;
    pV = (xb_vecN_2xf32 *)(v + 2 * L);
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IP(rsqrt, pZr, 2 * BBE_SIMD_WIDTH);
        
        X0 = BBE_LVN_2XF32_X(pX, 2 * N*L * sizeof(float32_t));
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_X(Acc, pV, 2 * L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

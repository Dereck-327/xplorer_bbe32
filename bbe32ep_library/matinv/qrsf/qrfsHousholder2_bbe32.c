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
  QR decomposition, floating point, real data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "qrsf_common.h"
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
// M==2
void  qrfsHousholder2(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX1;
    const xb_vecN_2xf32 * restrict pZr;
          xb_vecN_2xf32 * restrict pZ;
          xb_vecN_2xf32 * restrict pD;
          xb_vecN_2xf32 * restrict pF;
          xb_vecN_2xf32 * restrict pV;
          xb_vecN_2xf32 * restrict pV1;

    xb_vecN_2xf32 Acc, Acc1, X0, X1, F0, rsqrt;
    vboolN_2 tmp;

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(M == 2);

    // compute d=1/sqrt(x'x)
    pX = (const xb_vecN_2xf32 *)x;
    pZ = (xb_vecN_2xf32 *)pScr;
    pD = (xb_vecN_2xf32 *)D;
    pF = (xb_vecN_2xf32 *)Fi;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        X0 = BBE_LVN_2XF32_X(pX, N*L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(X0, X0);

        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(X0, X0);

        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SVN_2XF32_IP(X0, pF, 2 * BBE_SIMD_WIDTH);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_X(Acc1, pZ, L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
    }

    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pX = (const xb_vecN_2xf32 *)x;
    pX1 = (const xb_vecN_2xf32 *)(x + N*L);
    pZr = (const xb_vecN_2xf32 *)pScr;
    pV = (xb_vecN_2xf32 *)v;
    pV1 = (xb_vecN_2xf32 *)(v + L);
    pF = (xb_vecN_2xf32 *)Fi;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_LVN_2XF32_X(pZr, L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(X1, pZr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        BBE_MULSN_2XF32(X0, X1, F0);
        BBE_MULAN_2XF32(Acc, X0, X0);

        rsqrt = BBE_RSQRTN_2XF32(Acc);
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_IP(Acc1, pX1, 2 * BBE_SIMD_WIDTH);
        Acc1 = BBE_MULN_2XF32(Acc1, rsqrt);
        BBE_SVN_2XF32_IP(Acc1, pV1, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

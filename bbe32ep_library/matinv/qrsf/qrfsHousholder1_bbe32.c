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
/*
  special case: the same as above but for M==1
*/
 void qrfsHousholder1(           float32_t* restrict v,
                                 float32_t* restrict Fi,
                                 float32_t *restrict D,
                           const float32_t* restrict x, 
                           int SV, int M, int N, int L)
{
#if 0
    int l;
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(M == 1);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    (void)N, (void)M;
    for (l = 0; l < L; l++)
    {
        float32_t  x_re;
        x_re = x[l];
        Fi[l] = x_re > 0.f ? -1.f : 1.f;
        v[l] = 1.f;
        D[l] = 1.f / (fabsf(x_re));
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)D;
          xb_vecN_2xf32 * restrict pF = (      xb_vecN_2xf32 *)Fi;
          xb_vecN_2xf32 * restrict pV = (      xb_vecN_2xf32 *)v;

    xb_vecN_2xf32 X0, rsqrt;
    vboolN_2 tmp;

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(M == 1);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X0 = BBE_ABSN_2XF32(X0);
        rsqrt = BBE_RECIPN_2XF32(X0);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
    }
    pX = (const xb_vecN_2xf32 *)x;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SVN_2XF32_IP(X0, pF, 2 * BBE_SIMD_WIDTH);

        BBE_SVN_2XF32_IP(1.f, pV, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

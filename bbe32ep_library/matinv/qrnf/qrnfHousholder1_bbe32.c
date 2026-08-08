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
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)

    special case: M==1
-------------------------------------------------------*/
void qrnfHousholder1(float32_t* restrict v,
                     float32_t* restrict Fi,
                     float32_t *restrict D,
                     const float32_t* restrict x, 
                     int SD, int L)
{
#if 0
    int l;
    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++)
    {
        float32_t  x_re;
        x_re = x[l];
        Fi[l] = x_re > 0.f ? -1.f : 1.f;
        v[l] = 1.f;
        D[l*SD] = 1.f / (fabsf(x_re));
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xtfloat       * restrict pD = (      xtfloat       *)D;
          xb_vecN_2xf32 * restrict pF = (      xb_vecN_2xf32 *)Fi;
          xb_vecN_2xf32 * restrict pV = (      xb_vecN_2xf32 *)v;
          xb_vecN_2xf32 * restrict pZ = (      xb_vecN_2xf32 *)x;

    xb_vecN_2xf32 X0, rsqrt;
    valign vX, vF, vV;
    vboolN_2 tmp;

    vF = BBE_ZALIGN();
    vV = BBE_ZALIGN();
    for (l = 0; l < L / (BBE_SIMD_WIDTH / 2); l++)
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SAN_2XF32_IP(X0, vF, pF);

        BBE_SAN_2XF32_IP(1.f, vV, pV);
    }
    l = L % (BBE_SIMD_WIDTH / 2);
    if (l)
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        tmp = BBE_OGTN_2XF32(X0, BBE_ZERON_2XF32());
        X0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(-1.f), BBE_MOVNX16_FROMN_2XF32(1.f), BBE_MOVN_FROMN_2(tmp)));
        BBE_SAVN_2XF32_XP(X0, vF, pF, l * sizeof(float32_t));

        BBE_SAVN_2XF32_XP(1.f, vV, pV, l * sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vF, pF);
    BBE_SAN_2XF32POS_FP(vV, pV);

    pX = (const xb_vecN_2xf32 *)x;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X0 = BBE_ABSN_2XF32(X0);
        rsqrt = BBE_RECIPN_2XF32(X0);
        BBE_SVN_2XF32_IP(rsqrt, pZ, 2 * BBE_SIMD_WIDTH);
    }

    pX = (const xb_vecN_2xf32 *)x;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(X0, vX, pX, sizeof(float32_t));
        BBE_SSN_2XF32_XP(X0, pD, SD*sizeof(float32_t));
    }
}

#endif

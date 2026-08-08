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
void cqrnfHousholder1(       float32_t* restrict v,
                             float32_t* restrict Fi,
                             float32_t* restrict D,
                       const float32_t* restrict x, 
                       int SD, int L)
{
#if 0
    int l;
    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++)
    {
        float32_t f_re, f_im, x_re, x_im;
        float32_t d0, invd0;
        float32_t A_re, A_im;

        x_re = x[l * 2 + 0];
        x_im = x[l * 2 + 1];
        A_re = (x_re*x_re) + (x_im*x_im);
        d0 = 1.f / sqrtf(A_re);
        if (fabsf(A_re) <= FLT_MIN) d0 = 1.f;
        D[l*SD + 0] = D[l*SD + 1] = d0;
        // compute fi=x0/sqrt(x0'x0)
        f_re = Fi[2 * l + 0] = -x_re*d0;
        f_im = Fi[2 * l + 1] = -x_im*d0;
        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        invd0 = 1.f / d0;
        A_re = x_re;
        A_im = x_im;
        A_re -= (f_re*invd0);
        A_im -= (f_im*invd0);
        x_re = A_re;
        x_im = A_im;
        A_re = (x_re*x_re) + (x_im*x_im);
        d0 = 1.f / sqrtf(A_re);
        v[2 * l + 0] = x_re*d0;
        v[2 * l + 1] = x_im*d0;
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)x;
          long long     * restrict pD_= (      long long     *)D;
          xb_vecN_2xf32 * restrict pF = (      xb_vecN_2xf32 *)Fi;
          xb_vecN_2xf32 * restrict pV = (      xb_vecN_2xf32 *)v;
    const xb_vecN_2xf32 * restrict pVr= (const xb_vecN_2xf32 *)v;

    xb_vecN_2xf32 Acc, X0, X1, F0, rsqrt;
    vboolN_2 tmp;
    valign vV, vVr, vF, vX;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    vV = BBE_ZALIGN();
    vF = BBE_ZALIGN();
    for (l = L; l > 0; l -= (BBE_SIMD_WIDTH / 4))
    {
        // compute d=1/sqrt(x'x)
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        rsqrt = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(BBE_CONSTN_2XF32(1)), BBE_MOVNX16_FROMN_2XF32(rsqrt), BBE_MOVN_FROMN_2(tmp))); /////////////
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);

        // compute fi=x0/sqrt(x0'x0)
        F0 = BBE_MULN_2XF32(X0, rsqrt);
        F0 = BBE_NEGN_2XF32(F0);
        BBE_SAVN_2XF32_XP(F0, vF, pF, 2 * l*sizeof(float32_t));

        // prepare for update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_MULSN_2XF32(X0, F0, Acc);
        BBE_SAVN_2XF32_XP(X0, vV, pV, 2 * l*sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vF, pF);
    BBE_SAN_2XF32POS_FP(vV, pV);

    pX = (const xb_vecN_2xf32 *)x;
    for (l = 0; l < L; l++)
    {
        BBE_LAVN_2XF32_XP(Acc, vX, pX, 2 * sizeof(float32_t));
        BBE_SSN_4X64_XP(BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc)), pD_, SD*sizeof(float32_t));
    }

    pVr = (const xb_vecN_2xf32 *)v;
    pV = (xb_vecN_2xf32 *)v;
    vVr = BBE_LAN_2XF32_PP(pVr);
    vV = BBE_ZALIGN();
    for (l = L; l > 0; l -= (BBE_SIMD_WIDTH / 4))
    {
        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        BBE_LAN_2XF32_IP(X0, vVr, pVr);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        rsqrt = BBE_RSQRTN_2XF32(Acc);

        // normalization
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SAVN_2XF32_XP(Acc, vV, pV, 2 * l*sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vV, pV);
}
#endif

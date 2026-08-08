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
void cqrfsHousholder1(    float32_t* restrict v,
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
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);
    (void)N, (void)M;
    // compute d=1/sqrt(x'x)
    for (l = 0; l < L; l++)
    {
        float32_t  f_re, f_im, x_re, x_im;
        float32_t d0, invd0;
        float32_t A_re;

        x_re = x[l * 2 + 0];
        x_im = x[l * 2 + 1];
        A_re = (x_re*x_re) + (x_im*x_im);
        d0 = D[l * 2 + 0] = D[l * 2 + 1] = 1.f / sqrtf(A_re);
        // compute fi=x0/sqrt(x0'x0)
        if (fabsf(A_re) <= FLT_MIN) d0 = 1.f;
        f_re = -x_re*d0;
        f_im = -x_im*d0;
        Fi[l * 2 + 0] = f_re;
        Fi[l * 2 + 1] = f_im;
        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        invd0 = A_re*d0;
        x_re -= f_re*invd0;
        x_im -= f_im*invd0;
        A_re = (x_re*x_re) + (x_im*x_im);
        d0 = 1.f / sqrtf(A_re);
        // normalization
        v[l * 2 + 0] = x_re*d0;
        v[l * 2 + 1] = x_im*d0;
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)D;
          xb_vecN_2xf32 * restrict pF = (      xb_vecN_2xf32 *)Fi;
          xb_vecN_2xf32 * restrict pV = (      xb_vecN_2xf32 *)v;
    const xb_vecN_2xf32 * restrict pVr= (const xb_vecN_2xf32 *)v;

    xb_vecN_2xf32 Acc, X0, X1, F0, rsqrt;
    vboolN_2 tmp;

    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(M == 1);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        // compute d=1/sqrt(x'x)
        BBE_LVN_2XF32_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        tmp = BBE_OEQN_2XF32(Acc, BBE_ZERON_2XF32()); /////////////
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);
        rsqrt = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(BBE_CONSTN_2XF32(1)), BBE_MOVNX16_FROMN_2XF32(rsqrt), BBE_MOVN_FROMN_2(tmp))); /////////////

        // compute fi=x0/sqrt(x0'x0)
        F0 = BBE_MULN_2XF32(X0, rsqrt);
        F0 = BBE_NEGN_2XF32(F0);
        BBE_SVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        // prepare for update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_MULSN_2XF32(X0, F0, Acc);
        BBE_SVN_2XF32_IP(X0, pV, 2 * BBE_SIMD_WIDTH);
    }

    pVr = (const xb_vecN_2xf32 *)v;
    pV = (xb_vecN_2xf32 *)v;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
    {
        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        BBE_LVN_2XF32_IP(X0, pVr, 2 * BBE_SIMD_WIDTH);
        X1 = BBE_SHFLN_2XF32I(X0, BBE_SHFLI_SWAP_2);
        Acc = BBE_MULN_2XF32(X0, X0);
        BBE_MULAN_2XF32(Acc, X1, X1);
        rsqrt = BBE_RSQRTN_2XF32(Acc);
        
        // normalization
        Acc = BBE_MULN_2XF32(X0, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pV, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

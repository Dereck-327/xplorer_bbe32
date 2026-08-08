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
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "qrsf_common.h"

#if HAVE_VFPU

/*-----------------------------------------------------------
    partial update of R matrix - update begins from 
    v         pointer to m-th Householder vector 
    K         current number of column in the matrix R
    M,N       matrix size
    L         number of matrices
    Input/output:
    R         pointer to the beginning of column  
-----------------------------------------------------------*/
void qrfsUpdateR(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L)
{
#if 0
    int m, l, n;
    float32_t A_re, A_im;
    int _NL = N*L;
    float32_t v_re, b_re, z_re;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));

    for (n = 0; n < (N - K); n++)
    {
        for (l = 0; l < L; l++)
        {
            // compute v'B first
            A_re = A_im = 0;
            for (m = 0; m < M; m++)
            {
                v_re = v[m*L + l];
                b_re = R[m*_NL + l];
                A_re += (v_re*b_re);
            }
            z_re = 2.f*A_re;
            for (m = 0; m < M; m++)
            {
                v_re = v[m*L + l];
                A_re = R[m*_NL + l];
                A_re -= (z_re*v_re);
                R[m*_NL + l + 0] = A_re;
            }
        }
        R += L; // next column
    }
#endif // 0

    int m, l;

    const xb_vecN_2xf32 * restrict pRr;
          xb_vecN_2xf32 * restrict pR;
          xb_vecN_2xf32 * restrict pR0;
    const xb_vecN_2xf32 * restrict pV;
    const xb_vecN_2xf32 * restrict pV0;

    xb_vecN_2xf32 Acc, V0, R0, Z0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));

    WUR_CBEGIN((uintptr_t)v);
    WUR_CEND((uintptr_t)v + L * sizeof(float32_t));

    pR0 = (xb_vecN_2xf32 *)R;
    pV0 = (const xb_vecN_2xf32 *)v;
    for (l = 0; l < (N - K)*L; l += (BBE_SIMD_WIDTH / 2))
    {
        pV = pV0;
        pRr = pR0;
        Acc = BBE_ZERON_2XF32();
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(R0, pRr, N*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc, V0, R0);
        }
        pV = pV0;
        pR = pR0;
        pRr = pR0;
        Z0 = BBE_MULN_2XF32(Acc, 2.f);
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(R0, pRr, N*L * sizeof(float32_t));
            BBE_MULSN_2XF32(R0, Z0, V0);
            BBE_SVN_2XF32_XP(R0, pR, N*L * sizeof(float32_t));
        }
        BBE_LVN_2XF32_IC(V0, pV0);
        pR0 = (xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pR0);
    }
}
#endif

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
    rotate R[L][SA] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV] diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    Input/output:
    R[L][SA]  sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper trinagle matrix 
              are zeroed!
-------------------------------------------------------*/
void qrnfRotateR(float32_t* restrict R,const float32_t* restrict Fi,int N,int SA,int L)
{
#if 0
    int m, n, l;
    float32_t A_re;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % 4 == 0 && N > 0);
    NASSERT(L > 0);

    for (l = 0; l < L; l++)
    {
        for (m = 0; m < N; m++)
            for (n = 0; n < N; n++)
            {
                float32_t f_re, b_re;
                f_re = Fi[l + m*L];
                b_re = R[m*N + n];
                A_re = (f_re*b_re);
                R[m*N + n] = A_re;
            }
        R += SA;
    }
#endif // 0

    int m, n, l;
    int count;

    const xtfloat       * restrict pF;
    const xtfloat       * restrict pF_;
    const xb_vecN_2xf32 * restrict pR;
          xb_vecN_2xf32 * restrict pRw;
          xb_vecN_2xf32 * restrict pR_;

    xb_vecN_2xf32 Acc, F0, R0;
    valign vR, vRw;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % 4 == 0);
    NASSERT(L > 0);

    count = (N % (BBE_SIMD_WIDTH / 2))*sizeof(float32_t);
    pF_ = (const xtfloat *)Fi;
    pR_ = (xb_vecN_2xf32 *)R;
    for (l = 0; l < L; l++)
    {
        pR = pR_;
        pRw = pR_;
        pF = pF_;
        for (m = 0; m < N; m++)
        {
            BBE_LSN_2XF32_XP(F0, pF, L*sizeof(float32_t));
            F0 = BBE_REPN_2XF32(F0, 0);

            for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
            {
                BBE_LAN_2XF32_IP(R0, vR, pR);
                Acc = BBE_MULN_2XF32(F0, R0);
                BBE_SAN_2XF32_IP(Acc, vRw, pRw);
            }
            BBE_LAVN_2XF32_XP(R0, vR, pR, count);
            Acc = BBE_MULN_2XF32(F0, R0);
            BBE_SAVN_2XF32_XP(Acc, vRw, pRw, count);
        }
        pF_ = (const xtfloat *)XT_ADDX4(1, (uintptr_t)pF_);
        pR_ = (xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pR_);
    }
}
#endif


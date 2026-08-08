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

#if HAVE_VFPU

/*-----------------------------------------------------------
    rotate R[L][MxN] by diagonal matrix Fi'[L][N]
    Input:
    Fi[L][N]  diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    Input/output:
    R[L][MxN] sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper triangle matrix 
              remain unchanged
-----------------------------------------------------------*/
void qrfsRotateR(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L)
{
#if 0
    float32_t f_re, b_re;
    int m, n, l;
    float32_t A_re;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    for (m = 0; m < N; m++)
    {
        for (l = 0; l < L; l++)
        {
            f_re = Fi[l + m*L];
            for (n = 0; n < N; n++)
            {
                b_re = R[l + n*L];
                A_re = (f_re*b_re);
                R[l + n*L] = A_re;
            }
        }
        R += N*L;   // next row
    }
#endif // 0

    int m, l;

    const xb_vecN_2xf32 * restrict pRr;
          xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pF;

    xb_vecN_2xf32 Acc, F0, R0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    pRr = (const xb_vecN_2xf32 *)R;
    pR = (xb_vecN_2xf32 *)R;
    pF = (const xb_vecN_2xf32 *)Fi;
    for (m = 0; m < N; m++)
    {
        WUR_CBEGIN((uintptr_t)pF);
        WUR_CEND((uintptr_t)pF + L * sizeof(float32_t));

        for (l = 0; l < N*L; l += (BBE_SIMD_WIDTH / 2))
        {
            BBE_LVN_2XF32_IC(F0, pF);
            BBE_LVN_2XF32_XP(R0, pRr, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULN_2XF32(F0, R0);
            BBE_SVN_2XF32_XP(Acc, pR, 2 * BBE_SIMD_WIDTH);
        }
        pF = (const xb_vecN_2xf32 *)XT_ADDX4(L, (uintptr_t)pF);
    }
}
#endif

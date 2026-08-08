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
void cqrsfRotateBconj1(float32_t* B,const float32_t* Fi,int N,int L)
{
#if 0
    float32_t A_re, A_im, f_re, f_im, b_re, b_im;
    int n;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(N > 1 && L > 0 && (L % (BBE_SIMD_WIDTH / 4) == 0));
    for (n = 0; n < N*L; n++)
    {
        f_re = Fi[n * 2 + 0];
        f_im = Fi[n * 2 + 1];
        b_re = B[n * 2 + 0];
        b_im = B[n * 2 + 1];
        A_re = (f_re*b_re) + (f_im*b_im);
        A_im = (f_re*b_im) - (f_im*b_re);
        B[n * 2 + 0] = A_re;
        B[n * 2 + 1] = A_im;
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pF;

    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(N > 1 && L > 0 && (L % (BBE_SIMD_WIDTH / 4) == 0));

    pBr = (const xb_vecN_2xf32 *)B;
    pB = (xb_vecN_2xf32 *)B;
    pF = (const xb_vecN_2xf32 *)Fi;
    for (l = 0; l < N*L; l += (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_IP(B0, pBr, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
        BBE_SVN_2XF32_IP(Acc, pB, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

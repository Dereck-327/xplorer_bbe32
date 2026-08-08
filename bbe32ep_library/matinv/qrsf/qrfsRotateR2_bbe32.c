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
// N==2
void qrfsRotateR2(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pRr;
    const xb_vecN_2xf32 * restrict pRr0;
          xb_vecN_2xf32 * restrict pR;
          xb_vecN_2xf32 * restrict pR0;
    const xb_vecN_2xf32 * restrict pF;
    const xb_vecN_2xf32 * restrict pF0;

    xb_vecN_2xf32 Acc, F0, R0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(N == 2 && M >= 2);

    pRr = (const xb_vecN_2xf32 *)R;
    pRr0 = (const xb_vecN_2xf32 *)(R + 2 * L);
    pF = (const xb_vecN_2xf32 *)Fi;
    pF0 = (const xb_vecN_2xf32 *)(Fi + L);
    pR = (xb_vecN_2xf32 *)R;
    pR0 = (xb_vecN_2xf32 *)(R + 2 * L);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        BBE_LVN_2XF32_IP(F0, pF, 2 * BBE_SIMD_WIDTH);

        R0 = BBE_LVN_2XF32_X(pRr, L * sizeof(float32_t));
        Acc = BBE_MULN_2XF32(F0, R0);
        BBE_SVN_2XF32_X(Acc, pR, L * sizeof(float32_t));

        BBE_LVN_2XF32_IP(R0, pRr, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(F0, R0);
        BBE_SVN_2XF32_IP(Acc, pR, 2 * BBE_SIMD_WIDTH);


        BBE_LVN_2XF32_IP(F0, pF0, 2 * BBE_SIMD_WIDTH);

        R0 = BBE_LVN_2XF32_X(pRr0, L * sizeof(float32_t));
        Acc = BBE_MULN_2XF32(F0, R0);
        BBE_SVN_2XF32_X(Acc, pR0, L * sizeof(float32_t));

        BBE_LVN_2XF32_IP(R0, pRr0, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(F0, R0);
        BBE_SVN_2XF32_IP(Acc, pR0, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

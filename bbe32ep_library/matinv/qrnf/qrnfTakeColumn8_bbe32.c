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
    take column from sequence of block ordered matrices
    and put it to the linear array
    Input:
    A[L][SA]    matrices
    output:
    x [L*M]   contingious array
    xt[L*M]   the same but in the transposed form
-------------------------------------------------------*/
void qrnfTakeColumn8(float32_t* restrict x,float32_t* restrict xt,const float32_t* restrict A,int M,int N,int SA,int L)
{
#if 0
    int l, m;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);
    NASSERT(M >= 1 && M <= 8);
    for (l = 0; l < L; l++)
        for (m = 0; m < M; m++)
        {
            xt[m*L + l + 0] = x[l*M + m] = A[l*SA + m*N + 0];
        }
#endif // 0

    int l, m;

    const xtfloat       * restrict pA  = (const xtfloat       *)A;
          xtfloat       * restrict pXt = (      xtfloat       *)xt;
          xtfloat       * restrict pXt0= (      xtfloat       *)xt;
          xb_vecN_2xf32 * restrict pX  = (      xb_vecN_2xf32 *)x;

    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;
    valign vX;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(xt, 2 * BBE_SIMD_WIDTH);
    NASSERT(M >= 1 && M <= 8);

    for (l = 0; l < L - 1; l++)
    {
        BBE_LSN_2XF32_XP(t0, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t1, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t2, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t3, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t4, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t5, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t6, pA, N * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t7, pA, (SA - 7 * N) * sizeof(float32_t));
        
        t0 = BBE_SELN_2XF32I(t1, t0, BBE_SELI_INTERLEAVE_2_LO);
        t2 = BBE_SELN_2XF32I(t3, t2, BBE_SELI_INTERLEAVE_2_LO);
        t4 = BBE_SELN_2XF32I(t5, t4, BBE_SELI_INTERLEAVE_2_LO);
        t6 = BBE_SELN_2XF32I(t7, t6, BBE_SELI_INTERLEAVE_2_LO);
        t0 = BBE_SELN_2XF32I(t2, t0, BBE_SELI_INTERLEAVE_4_LO);
        t4 = BBE_SELN_2XF32I(t6, t4, BBE_SELI_INTERLEAVE_4_LO);
        t0 = BBE_SELN_2XF32I(t4, t0, BBE_SELI_EXTRACT_LO_HALVES);
        BBE_SAVN_2XF32_XP(t0, vX, pX, M * sizeof(float32_t));
    }
    for (m = 0; m < M; m++)
    {
        BBE_LSN_2XF32_XP(t0, pA, N * sizeof(float32_t));
        BBE_SAVN_2XF32_XP(t0, vX, pX, sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vX, pX);

    pX = (xb_vecN_2xf32 *)x;
    for (l = 0; l < L; l++)
    {
        pXt = pXt0;
        for (m = 0; m < M; m++)
        {
            BBE_LAVN_2XF32_XP(t0, vX, pX, sizeof(float32_t));
            BBE_SSN_2XF32_XP(t0, pXt, L * sizeof(float32_t));
        }
        pXt0 = (xtfloat *)XT_ADDX4(1, (uintptr_t)pXt0);
    }
}
#endif

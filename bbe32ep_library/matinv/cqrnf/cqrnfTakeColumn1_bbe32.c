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
#include "common.h"
#include "cqrnf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------
    take column from sequence of block ordered matrices
    and put it to the linear array
    M==1
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
void cqrnfTakeColumn1(float32_t* restrict x,const float32_t* restrict A,int SA,int L)
{
#if 0
    int l;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < L; l++)
    {
        x[l * 2 + 0] = A[l*SA + 0];
        x[l * 2 + 1] = A[l*SA + 1];
    }
#endif // 0

    int l;

    const long long     * restrict pA = (const long long     *)A;
          xb_vecN_2xf32 * restrict pX = (      xb_vecN_2xf32 *)x;

    xb_vecN_4x64 t0, t1, t2, t3;
    xb_vecN_2xf32 vTmp;
    valign vX;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);

    for (l = 0; l < L / (BBE_SIMD_WIDTH / 4); l++)
    {
        BBE_LSN_4X64_XP(t0, pA, SA * sizeof(float32_t));
        BBE_LSN_4X64_XP(t1, pA, SA * sizeof(float32_t));
        BBE_LSN_4X64_XP(t2, pA, SA * sizeof(float32_t));
        BBE_LSN_4X64_XP(t3, pA, SA * sizeof(float32_t));

        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        vTmp = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));
        BBE_SAN_2XF32_IP(vTmp, vX, pX);
    }
    for (l *= (BBE_SIMD_WIDTH / 4); l < L; l++)
    {
        BBE_LSN_4X64_XP(t0, pA, SA * sizeof(float32_t));

        BBE_SAVN_2XF32_XP(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(t0)), vX, pX, 2 * sizeof(float32_t));
    }
    BBE_SAN_2XF32POS_FP(vX, pX);
}

#endif

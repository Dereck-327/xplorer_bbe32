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
void cqrnfrotateR16(float32_t* restrict R,const float32_t* restrict Fi,int L)
{
    int m, l;
    const int N = 16;

    const long long     * restrict pF;
    const long long     * restrict pF_ = (const long long     *)Fi;
    const xb_vecN_2xf32 * restrict pR  = (const xb_vecN_2xf32 *)R;
          xb_vecN_2xf32 * restrict pRw = (      xb_vecN_2xf32 *)R;

    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 Acc, F0, R0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0);

    for (l = 0; l < L; l++)
    {
        pF = pF_;
        //BBE_LSN_4X64_XP(vTmp, pF, 2 * L*sizeof(float32_t));
        //F0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
        for (m = 0; m < N; m++)
        {
            BBE_LSN_4X64_XP(vTmp, pF, 2 * L*sizeof(float32_t));
            F0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULMN_2XF32(F0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, R0, 2, 11);
            BBE_SVN_2XF32_IP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULMN_2XF32(F0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, R0, 2, 11);
            BBE_SVN_2XF32_IP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULMN_2XF32(F0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, R0, 2, 11);
            BBE_SVN_2XF32_IP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULMN_2XF32(F0, R0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, R0, 2, 11);
            BBE_SVN_2XF32_IP(Acc, pRw, 2 * BBE_SIMD_WIDTH);

            //BBE_LSN_4X64_XP(vTmp, pF, 2 * L*sizeof(float32_t));
            //F0 = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(vTmp)), BBE_SHFLI_REP_0X4);
        }
        pF_ = (const long long *)XT_ADDX4(2, (uintptr_t)pF_);
    }
}
#endif

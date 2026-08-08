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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    Cholesky forward recursion, floating point complex data, 
    block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cholnf_common.h"

#if (HAVE_VFPU)

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    MxNxP=16x16x1
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void computeAB16x16x1f(float32_t* Z,const float32_t* A,const float32_t* B, int L)
{
#if 0
    float32_t B_re, B_im;
    int n, m, l;
    int SA = 2 * (16 * 16);
    int SB = 2 * (16 * 1);

    for (n = 0; n < 16; n++)
        for (l = 0; l < L; l++)
        {
            B_re = B_im = 0.f;
            for (m = 0; m < 16; m++)
            {
                float32_t a_re, a_im, b_re, b_im;
                a_re = A[l*SA + 2 * n + m * 16 * 2 + 0]; a_im = A[l*SA + 2 * n + m * 16 * 2 + 1];
                b_re = B[l*SB + m * 2 + 0]; b_im = B[l*SB + m * 2 + 1];
                B_re += (a_re*b_re) + (a_im*b_im);
                B_im += (a_re*b_im) - (a_im*b_re);
            }
            Z[2 * l * 16 + 2 * n + 0] = B_re;
            Z[2 * l * 16 + 2 * n + 1] = B_im;
        }
#endif // 0

    int l;
    const xb_vecN_2xf32 * restrict pA;
    const xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pZ;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, B0, temp;
    xb_vecN_2xf32 A0, A1, A2, A3;

    pA = (const xb_vecN_2xf32 *)(A);
    pB = (const xb_vecN_2xf32 *)(B);
    pZ = (xb_vecN_2xf32 *)(Z);
    for (l = 0; l < L; l++)
    {
        BBE_LVN_2XF32_XP(temp, pB, 2 * BBE_SIMD_WIDTH);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_0X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        Acc1 = BBE_MULMN_2XF32(A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        Acc2 = BBE_MULMN_2XF32(A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        Acc3 = BBE_MULMN_2XF32(A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_1X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_2X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_3X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        BBE_LVN_2XF32_XP(temp, pB, 2 * BBE_SIMD_WIDTH);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_0X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);
        
        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_1X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);
        
        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_2X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_3X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        BBE_LVN_2XF32_XP(temp, pB, 2 * BBE_SIMD_WIDTH);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_0X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_1X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_2X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_3X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        BBE_LVN_2XF32_XP(temp, pB, 2 * BBE_SIMD_WIDTH);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_0X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_1X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_2X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);

        B0 = BBE_SHFLN_2XF32I(temp, BBE_SHFLI_REP_3X4);
        BBE_LVN_2XF32_XP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A1, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A2, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(A3, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc, A0, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A1, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A2, B0, 2, 11);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A3, B0, 2, 11);


        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Acc1, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Acc2, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(Acc3, pZ, 2 * BBE_SIMD_WIDTH);
    }
}

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SA,SR,SD,SB,SY denote the number of data elements required to store
a matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next power
of two, otherwise it is matrix_size rounded up to the next multiple of the SIMD 
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SB = storage_size(M*P)
SY = storage_size(N*P)

Scratch size in bytes is defined by [r]cholfwd<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 M         Matrix dimension (number of rows in matrices A)
 N         Matrix dimension (number of columns and rows in 
           matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 A[L][SA]  Sequence of L complex matrices A
 B[L][SB]  Sequence of original right-side matrices B
 D[L][SD]  Reciprocal of main diagonal
Output:
 y[L][SY]   Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4 
4. M>=N
---------------------------------------------------------------------------*/
void cholfwd16x16x1nf(
            void *pScr,
            complex_float * restrict _y,
      const complex_float * restrict _R, 
      const complex_float * restrict _D,
      const complex_float * restrict _A, 
      const complex_float * restrict _B, 
            int L)
{
          float32_t* restrict y=(      float32_t*)_y;
    const float32_t* restrict R=(const float32_t*)_R;
    const float32_t* restrict D=(const float32_t*)_D;
    const float32_t* restrict A=(const float32_t*)_A;
    const float32_t* restrict B=(const float32_t*)_B;
    float32_t* Z=(float32_t* )pScr;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);

    // compute A'*B
    computeAB16x16x1f(Z,A,B, L);
    cholnfFwdrec16(y,R,D,Z,16,L,2*((16*(16+1))>>1),2*16,2*16,2*16);
}

size_t cholfwd16x16x1nf_getScratchSize(int M,int N, int P,int L) 
{
    size_t Z_size;
    NASSERT(M==16 && N==16 && P==1);
    L=XT_MAX(0,L);
    Z_size= (2*L*N*P*M*sizeof(float32_t));
    return (Z_size);
}

#else
DISCARD_FUN(void, cholfwd16x16x1nf,(
            void * pScr,
            complex_float * restrict _y,
      const complex_float * restrict _R, 
      const complex_float * restrict _D,
      const complex_float * restrict _A, 
      const complex_float * restrict _B, 
      int L ))

size_t cholfwd16x16x1nf_getScratchSize(int M,int N, int P,int L)
{
  (void)M; (void)N; (void)P; (void)L;
  return 0;
}

#endif

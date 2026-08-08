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
    Cholesky backward recursion, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
    */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cholnf_common.h"

#if (HAVE_VFPU)

static void cholnfTransformR16(float32_t* Rt, const float32_t* R, int N, int L)
{
    const int SR = N*(N + 1);
    int l;
    const long long * restrict pR;
    long long * restrict pRw = (long long *)Rt;
    const long long * restrict pR2 = (const long long *)R;
    xb_vecN_4x64 vTmp;

    __Pragma("loop_count min=1");
    for (l = 0; l < L; l++)
    {
        pR = (const long long *)XT_ADDX4(268, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(236, (uintptr_t)pR2);

        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(206, (uintptr_t)pR2);

        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(178, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(152, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(128, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(106, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR,                         0);
        BBE_SSN_4X64_IP(vTmp, pRw,    2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(86, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(68, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(52, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(38, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 6 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(26, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 5 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 6 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(16, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 4 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 5 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 6 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(8, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 3 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 4 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 5 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 6 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


        pR = (const long long *)XT_ADDX4(2, (uintptr_t)pR2);

        BBE_LSN_4X64_IP(vTmp, pR, 2 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 3 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 4 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 5 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 6 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 7 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 8 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 9 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 10 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 11 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 12 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 13 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 14 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_XP(vTmp, pR, 15 * 2 * sizeof(float32_t));
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        BBE_LSN_4X64_IP(vTmp, pR, 0);
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));




        pR2 = (const long long *)XT_ADDX4(SR, (uintptr_t)pR2);

        

    }
}

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SR,SD,SY,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SY = storage_size(N*P)
SX = storage_size(N*P)

Scratch size in bytes is defined by [r]cholbkw<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 N         Matrix dimension (number of columns and rows in matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal
 y[L][SY]  Sequence of intermediate decision matrices y
Output:         
 x[L][SX]  Sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4
---------------------------------------------------------------------------*/
void cholbkw16x1nf(
      void *pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
    int L)
{
          float32_t* restrict x=(      float32_t*)_x;
    const float32_t* restrict R=(const float32_t*)_R;
    const float32_t* restrict D=(const float32_t*)_D;
    const float32_t* restrict y=(const float32_t*)_y;
    float32_t* Rt = (float32_t*)pScr;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    cholnfTransformR16(Rt, R, 16, L);
    cholnfBkwnx1(x, Rt,D,y, 16,L);
}

/*
    return Scratch size
    Input:
    N:      matrix size
    P:      ignored
    L       number of matrices
    */
size_t cholbkw16x1nf_getScratchSize(int N, int P, int L)
{
    NASSERT(N==16 && P==1);
    L=XT_MAX(0,L);
    return (L*N*(N - 1)*sizeof(float32_t));
}

#else
DISCARD_FUN(void, cholbkw16x1nf,(
            void * pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
      int L ))

size_t cholbkw16x1nf_getScratchSize(int N, int P, int L)
{
  (void)N; (void)P; (void)L;
  return 0;
}

#endif

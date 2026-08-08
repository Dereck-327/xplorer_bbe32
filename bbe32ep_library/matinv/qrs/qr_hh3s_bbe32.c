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
  NatureDSP_Baseband library. QR-based matrix decomposition and inversion for streaming order
    helper routines for Housholder algorithm, len_x==3
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "qr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, qrHouseholder3,(  
    const int16_t * restrict A, /* input matrix         */
    int16_t * restrict  v ,     /* output streaming order       */     
    const int M, 
    const int N, 
    const int L
    ))

    DISCARD_FUN(void, qrUpdateR3,(  int16_t* restrict R,
    const int16_t * restrict V,     /*i  pointer to Householder vector, Q14 */
    int m,int M,int N, int L))
#else
/*
    find Housholder vector by column and transform that column
    len_x==3
    */
void qrHouseholder3(
    const int16_t * restrict A, /* input matrix         */
    int16_t * restrict  v,     /* output streaming order       */
    const int M,
    const int N,
    const int L
    )
{
    int p;
    const xb_vecNx16 *restrict x = (const xb_vecNx16 *)A;
    xb_vecNx16 *restrict pv = (xb_vecNx16 *)v;
    xb_vecNx16 x0, t, norm_x, _0x4000 = 0x4000, d0, _0 = 0, k_norm;
    xb_vecNx40 a, b;
    vboolN sign;
    vsaN d1, exp;
    int strideV = 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(M, N) - 3 + 1);

    NASSERT(M > 1);

    for (p = 0; p < L; p += BBE_SIMD_WIDTH)
    {
        BBE_LVNX16_XP(x0, x, 2 * L*N); a = BBE_MULNX16(x0, x0);
        BBE_LVNX16_XP(x0, x, 2 * L*N); BBE_MULANX16(a, x0, x0);
        BBE_LVNX16_XP(x0, x, -2 * L*N*(3 - 1));
        BBE_MULANX16(a, x0, x0);

        a = BBE_ADDNX40(a, a);
        d1 = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, d1);
        BBE_RSQRTLUNX40_0(b, t, x0, a);
        BBE_RSQRTLUNX40_1(b, t, x0, a);
        BBE_MULUUSNX16(b, x0, t);
        a = BBE_SRAINX40(b, 23);
        d0 = BBE_PACKLNX40(a);
        d1 = BBE_SUBSR1SAVSN(16 + 4, d1);
        BBE_LVNX16_XP(x0, x, 2 * L*N);
        a = BBE_MULUSRNX16(d0, x0, d1);
        norm_x = BBE_PACKVNX40(a, d1);
        sign = BBE_LTNX16(x0, _0);
        BBE_ADDSNX16F(norm_x, norm_x, _0x4000, sign);
        BBE_SUBSNX16T(norm_x, norm_x, _0x4000, sign);
        x0 = BBE_ABSSNX16(norm_x);
        x0 = BBE_ADDNX16(x0, x0);
        a = BBE_UNPKUNX16(x0);

        exp = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, exp);
        exp = BBE_SUBSR1SAVSN(35, exp);
        BBE_RSQRTLUNX40_0(b, t, x0, a);
        BBE_RSQRTLUNX40_1(b, t, x0, a);
        BBE_MULUUSNX16(b, x0, t);
        k_norm = BBE_PACKVNX40(b, exp);
        x0 = BBE_MULNX16PACKQ(k_norm, norm_x);
        BBE_SVNX16_IP(x0, pv, 2 * BBE_SIMD_WIDTH);

        BBE_LVNX16_XP(x0, x, 2 * L*N);
        a = BBE_MULUSRNX16(d0, x0, d1);
        x0 = BBE_PACKVNX40(a, d1);
        x0 = BBE_MULNX16PACKQ(k_norm, x0);
        BBE_SVNX16_IP(x0, pv, 2 * BBE_SIMD_WIDTH);

        BBE_LVNX16_XP(x0, x, 2 * BBE_SIMD_WIDTH - 2 * L*N*(3 - 1));
        a = BBE_MULUSRNX16(d0, x0, d1);
        x0 = BBE_PACKVNX40(a, d1);
        x0 = BBE_MULNX16PACKQ(k_norm, x0);
        BBE_SVNX16_XP(x0, pv, strideV);
    }
}

/*
    update columns of matrix R by Housholder vectors
    len_x==3
    */
void qrUpdateR3(int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int L)
{
    int strideV;
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r0, r1, r2, v0, v1, v2, vr, _0x4000 = 0x4000;
    int l, j, sizeV = (SIZE_OF_V(M, N));
    pv = (const xb_vecNx16*)(V);
    pr = (const xb_vecNx16*)(R);
    pw = (xb_vecNx16*)(R);
    for (j = l = 0; l < (N - m)*((L >> (LOG2_BBE_SIMD_WIDTH))); l++)
    {
        j = BBE_ADDMOD16U(j, ((L >> (LOG2_BBE_SIMD_WIDTH)) << 16) | 1);
        strideV = 2 * BBE_SIMD_WIDTH*(sizeV - 2);
        XT_MOVEQZ(strideV, 2 * BBE_SIMD_WIDTH*(sizeV - 2) - 2 * L*sizeV, j);
        BBE_LVNX16_XP(r0, pr, 2 * N*L);
        BBE_LVNX16_XP(r1, pr, 2 * N*L);
        BBE_LVNX16_XP(r2, pr, 2 * BBE_SIMD_WIDTH - 2 * 2 * N*L);
        BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(v2, pv, strideV);
        acc = BBE_MULRNX16(r0, v0, q14);
        BBE_MULANX16(acc, r1, v1);
        BBE_MULANX16(acc, r2, v2);
        vr = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r2, _0x4000, q14); BBE_MULSNX16(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q14);
        BBE_SVNX16_XP(r0, pw, 2 * N*L);
        BBE_SVNX16_XP(r1, pw, 2 * N*L);
        BBE_SVNX16_XP(r2, pw, 2 * BBE_SIMD_WIDTH - 2 * 2 * N*L);
    }
}
#endif

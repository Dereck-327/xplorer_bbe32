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
    helper routines for Housholder algorithm, len_x==13
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
#include "cqr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void,cqrUpdateR13,(  int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m,int M,int N, int L))
#else
/*
    update columns of matrix R by Housholder vectors
    len_x==13
    */
void cqrUpdateR13(int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int L)
{
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, ra, rb, rc, v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, va, vb, vc, vr, _0x4000 = 0x4000;
    int l, j, sizeV = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    int strideV;
    pr = (const xb_vecNx16*)(R);
    pw = (xb_vecNx16*)(R);
    pv = (const xb_vecNx16*)(V);
    for (j = l = 0; l < (N - m)*((L >> (LOG2_BBE_SIMD_WIDTH - 1))); l++)
    {
        j = BBE_ADDMOD16U(j, ((L >> (LOG2_BBE_SIMD_WIDTH - 1)) << 16) | 1);
        strideV = 2 * BBE_SIMD_WIDTH*(sizeV - 12);
        XT_MOVEQZ(strideV, 2 * BBE_SIMD_WIDTH*(sizeV - 12) - 4 * L*sizeV, j);
        BBE_LVNX16_XP(r0, pr, 4 * N*L);
        BBE_LVNX16_XP(r1, pr, 4 * N*L);
        BBE_LVNX16_XP(r2, pr, 4 * N*L);
        BBE_LVNX16_XP(r3, pr, 4 * N*L);
        BBE_LVNX16_XP(r4, pr, 4 * N*L);
        BBE_LVNX16_XP(r5, pr, 4 * N*L);
        BBE_LVNX16_XP(r6, pr, 4 * N*L);
        BBE_LVNX16_XP(r7, pr, 4 * N*L);
        BBE_LVNX16_XP(r8, pr, 4 * N*L);
        BBE_LVNX16_XP(r9, pr, 4 * N*L);
        BBE_LVNX16_XP(ra, pr, 4 * N*L);
        BBE_LVNX16_XP(rb, pr, 4 * N*L);
        BBE_LVNX16_XP(rc, pr, -12 * 4 * N*L);
        BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v4, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v5, pv, 2 * BBE_SIMD_WIDTH);
        v6 = BBE_LVNX16_I(pv, 0 * 2 * BBE_SIMD_WIDTH);
        v7 = BBE_LVNX16_I(pv, 1 * 2 * BBE_SIMD_WIDTH);
        v8 = BBE_LVNX16_I(pv, 2 * 2 * BBE_SIMD_WIDTH);
        v9 = BBE_LVNX16_I(pv, 3 * 2 * BBE_SIMD_WIDTH);
        va = BBE_LVNX16_I(pv, 4 * 2 * BBE_SIMD_WIDTH);
        vb = BBE_LVNX16_I(pv, 5 * 2 * BBE_SIMD_WIDTH);
        vc = BBE_LVNX16_I(pv, 6 * 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r0, v0, q14);
        BBE_MULANX16J(acc, r1, v1);
        BBE_MULANX16J(acc, r2, v2);
        BBE_MULANX16J(acc, r3, v3);
        BBE_MULANX16J(acc, r4, v4);
        BBE_MULANX16J(acc, r5, v5);
        BBE_MULANX16J(acc, r6, v6);
        BBE_MULANX16J(acc, r7, v7);
        BBE_MULANX16J(acc, r8, v8);
        BBE_MULANX16J(acc, r9, v9);
        BBE_MULANX16J(acc, ra, va);
        BBE_MULANX16J(acc, rb, vb);
        BBE_MULANX16J(acc, rc, vc);
        vr = BBE_PACKVNX40(acc, q14);
        BBE_LVNX16_XP(r0, pr, 4 * N*L);
        BBE_LVNX16_XP(r1, pr, 4 * N*L);
        BBE_LVNX16_XP(r2, pr, 4 * N*L);
        BBE_LVNX16_XP(r3, pr, 4 * N*L);
        BBE_LVNX16_XP(r4, pr, 4 * N*L);
        BBE_LVNX16_XP(r5, pr, 4 * N*L);
        BBE_LVNX16_XP(r6, pr, 4 * N*L);
        BBE_LVNX16_XP(r7, pr, 4 * N*L);
        BBE_LVNX16_XP(r8, pr, 4 * N*L);
        BBE_LVNX16_XP(r9, pr, 4 * N*L);
        BBE_LVNX16_XP(ra, pr, 4 * N*L);
        BBE_LVNX16_XP(rb, pr, 4 * N*L);
        BBE_LVNX16_XP(rc, pr, 2 * BBE_SIMD_WIDTH - 12 * 4 * N*L);
        BBE_LVNX16_IP(v6, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v7, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v8, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v9, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(va, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vb, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(vc, pv, strideV);
        acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r2, _0x4000, q14); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r3, _0x4000, q14); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r4, _0x4000, q14); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r5, _0x4000, q14); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r6, _0x4000, q14); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r7, _0x4000, q14); BBE_MULSNX16C(acc, v7, vr); r7 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r8, _0x4000, q14); BBE_MULSNX16C(acc, v8, vr); r8 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r9, _0x4000, q14); BBE_MULSNX16C(acc, v9, vr); r9 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(ra, _0x4000, q14); BBE_MULSNX16C(acc, va, vr); ra = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(rb, _0x4000, q14); BBE_MULSNX16C(acc, vb, vr); rb = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(rc, _0x4000, q14); BBE_MULSNX16C(acc, vc, vr); rc = BBE_PACKVNX40(acc, q14);
        BBE_SVNX16_XP(r0, pw, 4 * N*L);
        BBE_SVNX16_XP(r1, pw, 4 * N*L);
        BBE_SVNX16_XP(r2, pw, 4 * N*L);
        BBE_SVNX16_XP(r3, pw, 4 * N*L);
        BBE_SVNX16_XP(r4, pw, 4 * N*L);
        BBE_SVNX16_XP(r5, pw, 4 * N*L);
        BBE_SVNX16_XP(r6, pw, 4 * N*L);
        BBE_SVNX16_XP(r7, pw, 4 * N*L);
        BBE_SVNX16_XP(r8, pw, 4 * N*L);
        BBE_SVNX16_XP(r9, pw, 4 * N*L);
        BBE_SVNX16_XP(ra, pw, 4 * N*L);
        BBE_SVNX16_XP(rb, pw, 4 * N*L);
        BBE_SVNX16_XP(rc, pw, 2 * BBE_SIMD_WIDTH - 12 * 4 * N*L);
    }
}
#endif

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
    cqr_calc_qb7x7xPs
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

#define SIZE_OF_V7x7    SIZE_OF_V(7,7)

/*-------------------------------------------------------------------------
cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
cqr_build_rMxN/qr_build_rMxN).

Fixed point representation of output matrices is the same as for input
matrices.

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L][C]                             Matrices B (L matrices of size MxP)
V[C*CQR_SIZE_V(M,N,L)/sizeof(int16_t)]  L sets of Householder vectors
C                                       1 for real, 2 for complex data
Output:
B[M*P][L][C]                            Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/

int cqr_calc_qb7x7xps(int16_t *B, const int16_t *V, int P, int L)
{
    int p, l;
    vsaN q = BBE_MOVVSA32(14);
    xb_vecNx16 vr, _0x4000 = 0x4000;
    xb_vecNx16 v0, v1, v2, v3, v4, v5, v6;
    xb_vecNx16 r0, r1, r2, r3, r4, r5, r6;
    xb_vecNx16 fi0, fi1, fi2, fi3, fi4, fi5, fi6;
    xb_vecNx40 acc;
    xb_vecNx16 * restrict _pw;
    const xb_vecNx16 * restrict _pr;
    const xb_vecNx16 * restrict _pv;
    const xb_vecNx16 *  restrict pFi;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    for (p = 0; p < P; p++)
    {
        _pv = (const xb_vecNx16 *)V;
        _pr = (const xb_vecNx16 *)B;
        _pw = (xb_vecNx16 *)B;
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(v6, _pv, 28 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, _pr, 4 * L*P);
            BBE_LVNX16_XP(r1, _pr, 4 * L*P);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -6 * 4 * L*P);
            acc = BBE_MULRNX16J(r0, v0, q);
            BBE_MULANX16J(acc, r1, v1);
            BBE_MULANX16J(acc, r2, v2);
            BBE_MULANX16J(acc, r3, v3);
            BBE_MULANX16J(acc, r4, v4);
            BBE_MULANX16J(acc, r5, v5);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            BBE_LVNX16_XP(r0, _pr, 4 * L*P);
            BBE_LVNX16_XP(r1, _pr, 4 * L*P);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -6 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r0, _pw, 4 * L*P);
            BBE_SVNX16_XP(r1, _pw, 4 * L*P);
            BBE_SVNX16_XP(r2, _pw, 4 * L*P);
            BBE_SVNX16_XP(r3, _pw, 4 * L*P);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);
            BBE_SVNX16_XP(r6, _pw, -6 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        // 2-nd column
        __Pragma("no_reorder")
            _pv = (const xb_vecNx16 *)(V + 7 * BBE_SIMD_WIDTH);
        _pr = (const xb_vecNx16 *)(B + 1 * 2 * L*P);
        _pw = (xb_vecNx16 *)_pr;
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(v6, _pv, 29 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r1, _pr, 4 * L*P);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -5 * 4 * L*P);
            acc = BBE_MULRNX16J(r1, v1, q);
            BBE_MULANX16J(acc, r2, v2);
            BBE_MULANX16J(acc, r3, v3);
            BBE_MULANX16J(acc, r4, v4);
            BBE_MULANX16J(acc, r5, v5);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            BBE_LVNX16_XP(r1, _pr, 4 * L*P);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -5 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r1, _pw, 4 * L*P);
            BBE_SVNX16_XP(r2, _pw, 4 * L*P);
            BBE_SVNX16_XP(r3, _pw, 4 * L*P);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);
            BBE_SVNX16_XP(r6, _pw, -5 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        // 3-rd column
        __Pragma("no_reorder")
            _pv = (const xb_vecNx16 *)(V + 13 * BBE_SIMD_WIDTH);
        _pr = (const xb_vecNx16 *)(B + 2 * 2 * L*P);
        _pw = (xb_vecNx16 *)_pr;
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(v6, _pv, 30 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -4 * 4 * L*P);
            acc = BBE_MULRNX16J(r2, v2, q);
            BBE_MULANX16J(acc, r3, v3);
            BBE_MULANX16J(acc, r4, v4);
            BBE_MULANX16J(acc, r5, v5);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -4 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r2, _pw, 4 * L*P);
            BBE_SVNX16_XP(r3, _pw, 4 * L*P);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);
            BBE_SVNX16_XP(r6, _pw, -4 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        // 4-th column
        __Pragma("no_reorder")
            _pv = (const xb_vecNx16 *)(V + 18 * BBE_SIMD_WIDTH);
        _pr = (const xb_vecNx16 *)(B + 3 * 2 * L*P);
        _pw = (xb_vecNx16 *)_pr;
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(v6, _pv, 31 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -3 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16J(r3, v3, q);
            BBE_MULANX16J(acc, r4, v4);
            BBE_MULANX16J(acc, r5, v5);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r3, _pw, 4 * L*P);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);
            BBE_SVNX16_XP(r6, _pw, -3 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        // 5 & 6-th columns
        __Pragma("no_reorder")
            _pv = (const xb_vecNx16 *)(V + 22 * BBE_SIMD_WIDTH);
        _pr = (const xb_vecNx16 *)(B + 4 * 2 * L*P);
        _pw = (xb_vecNx16 *)_pr;
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(v6, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            BBE_LVNX16_XP(r6, _pr, -2 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16J(r4, v4, q);
            BBE_MULANX16J(acc, r5, v5);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);

            BBE_LVNX16_IP(v5, _pv, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(v6, _pv, 30 * 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16J(r5, v5, q);
            BBE_MULANX16J(acc, r6, v6);
            vr = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r5, _0x4000, q); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q);
            acc = BBE_MULRNX16(r6, _0x4000, q); BBE_MULSNX16C(acc, v6, vr); r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);
            BBE_SVNX16_XP(r6, _pw, -2 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        // final rotation
        __Pragma("no_reorder")
            _pr = (const xb_vecNx16 *)B;
        _pw = (xb_vecNx16 *)B;
        pFi = (const xb_vecNx16*)(V + BBE_SIMD_WIDTH*SIZE_OF_V7x7);
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r0, fi0, q);
            r0 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r0, _pw, 4 * L*P);

            BBE_LVNX16_IP(fi1, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r1, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r1, fi1, q);
            r1 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r1, _pw, 4 * L*P);

            BBE_LVNX16_IP(fi2, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r2, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r2, fi2, q);
            r2 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r2, _pw, 4 * L*P);

            BBE_LVNX16_IP(fi3, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r3, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r3, fi3, q);
            r3 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r3, _pw, 4 * L*P);

            BBE_LVNX16_IP(fi4, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r4, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r4, fi4, q);
            r4 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r4, _pw, 4 * L*P);

            BBE_LVNX16_IP(fi5, pFi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r5, _pr, 4 * L*P);
            acc = BBE_MULRNX16J(r5, fi5, q);
            r5 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r5, _pw, 4 * L*P);

            BBE_LVNX16_XP(fi6, pFi, (SIZE_OF_V7x7 + 1) * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r6, _pr, -6 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
            acc = BBE_MULRNX16J(r6, fi6, q);
            r6 = BBE_PACKVNX40(acc, q);
            BBE_SVNX16_XP(r6, _pw, -6 * 4 * L*P + 2 * BBE_SIMD_WIDTH);
        }
        B = (int16_t*)XT_ADD((uintptr_t)B, 4 * L);
    }

    return 0;
} /* cqr_calc_qb6x6xps() */

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
    cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#define SIZE_OF_V4x2  9
#define SIZE_OF_FI4x2 4

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
B[M*P]L]                   Matrices B (L matrices of size MxP)
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
Output:
B[M*P][L]                  Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
void cqr_calc_qb4x2x1s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    int l;
    vsaN q = BBE_MOVVSA32(14);
    xb_vecNx16 vr, _0x4000 = 0x4000;
    xb_vecNx16 v0, v1, v2, v3;
    xb_vecNx16 r0, r1, r2, r3;
    xb_vecNx16 fi0, fi1, fi2, fi3;
    xb_vecNx40 acc;
    xb_vecNx16 * restrict _pw;
    const xb_vecNx16 * restrict _pr;
    const xb_vecNx16 * restrict _pv;
    const xb_vecNx16 *  restrict pFi;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    _pv = (const xb_vecNx16 *)V;
    _pr = (const xb_vecNx16 *)B;
    _pw = (xb_vecNx16 *)B;
    pFi = (const xb_vecNx16*)(V + BBE_SIMD_WIDTH*SIZE_OF_V4x2);
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        // 1-st column
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(r0, _pr, 4 * L);
        BBE_LVNX16_XP(r1, _pr, 4 * L);
        BBE_LVNX16_XP(r2, _pr, 4 * L);
        BBE_LVNX16_XP(r3, _pr, -3 * 4 * L + 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r0, v0, q);
        BBE_MULANX16J(acc, r1, v1);
        BBE_MULANX16J(acc, r2, v2);
        BBE_MULANX16J(acc, r3, v3);
        vr = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
        // second column
        BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, _pv, (SIZE_OF_FI4x2 + SIZE_OF_V4x2 - 7 + 1) * 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r1, v1, q);
        BBE_MULANX16J(acc, r2, v2);
        BBE_MULANX16J(acc, r3, v3);
        vr = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q);
        // final rotation
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r0, fi0, q);
        r0 = BBE_PACKVNX40(acc, q);
        BBE_SVNX16_XP(r0, _pw, 4 * L);

        BBE_LVNX16_IP(fi1, pFi, 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r1, fi1, q);
        r1 = BBE_PACKVNX40(acc, q);
        BBE_SVNX16_XP(r1, _pw, 4 * L);

        BBE_LVNX16_IP(fi2, pFi, 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r2, fi2, q);
        r2 = BBE_PACKVNX40(acc, q);
        BBE_SVNX16_XP(r2, _pw, 4 * L);

        BBE_LVNX16_IP(fi3, pFi, (SIZE_OF_V4x2 + 1) * 2 * BBE_SIMD_WIDTH);
        acc = BBE_MULRNX16J(r3, fi3, q);
        r3 = BBE_PACKVNX40(acc, q);
        BBE_SVNX16_XP(r3, _pw, -3 * 4 * L + 2 * BBE_SIMD_WIDTH);
    }
} /* cqr_calc_qb4x2x1s() */

size_t cqr_calc_qb4x2x1s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return 0;
} /* cqr_calc_qb4x2x1s_getScratchSize() */

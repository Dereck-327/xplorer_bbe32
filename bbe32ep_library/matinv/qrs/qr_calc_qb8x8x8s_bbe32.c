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
/* Common utility declarations. */
#include "qr_common.h"

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
void  qr_calc_qb8x8x8s (void *pScr, int16_t *B, const int16_t *V , int L)
{
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L%BBE_SIMD_WIDTH == 0);
    (void)pScr;
    qr_calc_qb8x8xps(B, V, 8, L);
} /* qr_calc_qb8x8x8s() */

size_t  qr_calc_qb8x8x8s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return 0;
} /* qr_calc_qb8x8x8s_getScratchSize() */



/*
    calcqb function for 8x8xP matrices
*/
void qr_calc_qb8x8xps(int16_t *B, const int16_t *V, int P, int L)
{
    vsaN q = BBE_MOVVSA32(14);
    xb_vecNx16 vr, _0x4000 = 0x4000, r0, v0;
    xb_vecNx40 acc;
    int p, l;
    int stride, vstride;
    xb_vecNx16 * restrict _pw;
    const xb_vecNx16 * restrict _pr;
    const xb_vecNx16 * restrict _pv;

    NASSERT_ALIGN(B, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(V, BBE_SIMD_WIDTH * 2);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH) == 0);

    stride = -7 * 2 * L*P + 2 * BBE_SIMD_WIDTH;

    //m=0;
    _pr = (const xb_vecNx16*)(B);
    _pw = (xb_vecNx16*)_pr;
    _pv = (const xb_vecNx16*)V;
    for (p = l = 0; l < P*L; l += BBE_SIMD_WIDTH)
    {
        p = BBE_ADDMOD16U(p, (L << 16) | BBE_SIMD_WIDTH);
        vstride = (7 + 6 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH;
        XT_MOVEQZ(vstride, (7 + 6 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH - 35 * 2 * L, p);

        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 0 * 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 6 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, -7 * 2 * L*P); v0 = BBE_LVNX16_I(_pv, 7 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, stride); BBE_LVNX16_XP(v0, _pv, vstride); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, stride);
    }
    __Pragma("no_reorder")
        //m=1;
        B = (int16_t*)XT_ADD(2 * L*P, (uintptr_t)B);
    _pr = (const xb_vecNx16*)(B);
    _pw = (xb_vecNx16*)_pr;
    stride = XT_ADD(2 * L*P, stride);
    _pv = (const xb_vecNx16*)XT_ADD(8 * 2 * BBE_SIMD_WIDTH, (uintptr_t)V);
    for (p = l = 0; l < P*L; l += BBE_SIMD_WIDTH)
    {
        p = BBE_ADDMOD16U(p, (L << 16) | BBE_SIMD_WIDTH);
        vstride = (8 + 0 + 6 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH;
        XT_MOVEQZ(vstride, (8 + 0 + 6 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH - 35 * 2 * L, p);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 0 * 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, -6 * 2 * L*P); v0 = BBE_LVNX16_I(_pv, 6 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, stride); BBE_LVNX16_XP(v0, _pv, vstride); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, stride);
    }
    __Pragma("no_reorder")
        //m=2;
        B = (int16_t*)XT_ADD(2 * L*P, (uintptr_t)B);
    _pr = (const xb_vecNx16*)(B);
    _pw = (xb_vecNx16*)_pr;
    stride = XT_ADD(2 * L*P, stride);
    _pv = (const xb_vecNx16*)XT_ADD((8 + 7) * 2 * BBE_SIMD_WIDTH, (uintptr_t)V);
    for (p = l = 0; l < P*L; l += BBE_SIMD_WIDTH)
    {
        p = BBE_ADDMOD16U(p, (L << 16) | BBE_SIMD_WIDTH);
        vstride = (8 + 7 + 0 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH;
        XT_MOVEQZ(vstride, (8 + 7 + 0 + 5 + 4 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH - 35 * 2 * L, p);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 0 * 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, -5 * 2 * L*P); v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_LVNX16_XP(r0, _pr, stride); BBE_LVNX16_XP(v0, _pv, vstride); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, stride);
    }
    __Pragma("no_reorder")
        //m=3...4
        B = (int16_t*)XT_ADD(2 * L*P, (uintptr_t)B);
    _pr = (const xb_vecNx16*)(B);
    _pw = (xb_vecNx16*)_pr;
    stride = XT_ADD(2 * L*P, stride);
    _pv = (const xb_vecNx16*)XT_ADD((8 + 7 + 6) * 2 * BBE_SIMD_WIDTH, (uintptr_t)V);
    for (p = l = 0; l < P*L; l += BBE_SIMD_WIDTH)
    {
        p = BBE_ADDMOD16U(p, (L << 16) | BBE_SIMD_WIDTH);
        vstride = (8 + 7 + 6 + 0 + 0 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH;
        XT_MOVEQZ(vstride, (8 + 7 + 6 + 0 + 0 + 3 + 2 + 1) * 2 * BBE_SIMD_WIDTH - 35 * 2 * L, p);
        xb_vecNx16 r0, r1, r2, r3, r4, v0, v1, v2, v3, v4;
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, v0, q);
        BBE_LVNX16_XP(r1, _pr, 2 * L*P); BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r1, v1);
        BBE_LVNX16_XP(r2, _pr, 2 * L*P); BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r2, v2);
        BBE_LVNX16_XP(r3, _pr, 2 * L*P); BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r3, v3);
        BBE_LVNX16_XP(r4, _pr, stride); BBE_LVNX16_XP(v4, _pv, -4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r4, v4);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v4, _pv, 2 * BBE_SIMD_WIDTH);

        acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16(acc, vr, v1); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16(acc, vr, v2); r2 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16(acc, vr, v3); r3 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16(acc, vr, v4); r4 = BBE_PACKVNX40(acc, q);

        BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r1, v1, q);
        BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r2, v2);
        BBE_LVNX16_IP(v3, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r3, v3);
        BBE_LVNX16_XP(v4, _pv, vstride); BBE_MULANX16(acc, r4, v4);
        vr = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16(acc, vr, v1); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16(acc, vr, v2); r2 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r3, _0x4000, q); BBE_MULSNX16(acc, vr, v3); r3 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r4, _0x4000, q); BBE_MULSNX16(acc, vr, v4); r4 = BBE_PACKVNX40(acc, q);

        BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_SVNX16_XP(r1, _pw, 2 * L*P);
        BBE_SVNX16_XP(r2, _pw, 2 * L*P);
        BBE_SVNX16_XP(r3, _pw, 2 * L*P);
        BBE_SVNX16_XP(r4, _pw, stride);
    }
    __Pragma("no_reorder")
        //m=5..6
        B = (int16_t*)XT_ADDX2(2 * L*P, (uintptr_t)B);
    _pr = (const xb_vecNx16*)(B);
    _pw = (xb_vecNx16*)_pr;
    stride = XT_ADDX2(2 * L*P, stride);
    _pv = (const xb_vecNx16*)XT_ADD((8 + 7 + 6 + 5 + 4) * 2 * BBE_SIMD_WIDTH, (uintptr_t)V);
    for (p = l = 0; l < P*L; l += BBE_SIMD_WIDTH)
    {
        p = BBE_ADDMOD16U(p, (L << 16) | BBE_SIMD_WIDTH);
        vstride = (8 + 7 + 6 + 5 + 4 + 0 + 0 + 1) * 2 * BBE_SIMD_WIDTH;
        XT_MOVEQZ(vstride, (8 + 7 + 6 + 5 + 4 + 0 + 0 + 1) * 2 * BBE_SIMD_WIDTH - 35 * 2 * L, p);
        xb_vecNx16 r0, r1, r2, v0, v1, v2;
        BBE_LVNX16_XP(r0, _pr, 2 * L*P); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, v0, q);
        BBE_LVNX16_XP(r1, _pr, 2 * L*P); BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r1, v1);
        BBE_LVNX16_XP(r2, _pr, stride); BBE_LVNX16_IP(v2, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16(acc, r2, v2);
        vr = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16(acc, vr, v0); r0 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16(acc, vr, v1); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16(acc, vr, v2); r2 = BBE_PACKVNX40(acc, q);

        BBE_LVNX16_IP(v1, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r1, v1, q);
        BBE_LVNX16_XP(v2, _pv, vstride); BBE_MULANX16(acc, r2, v2);
        vr = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r1, _0x4000, q); BBE_MULSNX16(acc, vr, v1); r1 = BBE_PACKVNX40(acc, q);
        acc = BBE_MULRNX16(r2, _0x4000, q); BBE_MULSNX16(acc, vr, v2); r2 = BBE_PACKVNX40(acc, q);

        BBE_SVNX16_XP(r0, _pw, 2 * L*P);
        BBE_SVNX16_XP(r1, _pw, 2 * L*P);
        BBE_SVNX16_XP(r2, _pw, stride);
    }
}

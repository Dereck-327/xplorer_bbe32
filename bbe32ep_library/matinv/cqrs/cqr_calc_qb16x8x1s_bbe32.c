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
#include "cqr_common.h"

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
void cqr_calc_qb16x8x1s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    //int M=16,N=8;
    vsaN q = BBE_MOVVSA32(14);
    xb_vecNx16 vr, _0x4000 = 0x4000;
    int l, sizeV, strideV, strideR, vinc;
    int voff;  // sum of arith progression (M,M-1...)
    xb_vecNx16 * restrict _pw;
    const xb_vecNx16 * restrict _pr;
    const xb_vecNx16 * restrict _pv;
    const xb_vecNx16 *  restrict pFi;
    xb_vecNx16 v0, r0;
    xb_vecNx40 acc;

    NASSERT_ALIGN(B, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(V, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    sizeV = (SIZE_OF_V(16, 8) + SIZE_OF_FI(16)) * 2 * BBE_SIMD_WIDTH;
    vinc = 16 * 2 * BBE_SIMD_WIDTH;
    voff = 0;
    strideV = -vinc;
    strideR = -16 * 4 * L;
    _pr = (const xb_vecNx16*)B;
    _pw = (xb_vecNx16*)_pr;

    //m=0
    strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    _pv = (const xb_vecNx16 *)V;
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=1
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=2
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=3
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")

        //m=4
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=5
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=6
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        //m=7
        strideV = XT_ADD(2 * BBE_SIMD_WIDTH, strideV);
    strideR = XT_ADD(4 * L, strideR);
    voff += vinc;
    _pv = (const xb_vecNx16 *)XT_ADD(voff, (uintptr_t)V);
    vinc = XT_SUB(vinc, 2 * BBE_SIMD_WIDTH);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16J(r0, v0, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, 4 * L);  BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc, r0, v0);
        BBE_LVNX16_XP(r0, _pr, strideR); BBE_LVNX16_XP(v0, _pv, strideV); BBE_MULANX16J(acc, r0, v0);
        vr = BBE_PACKVNX40(acc, q);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, 4 * L); BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(r0, _pr, strideR + 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(v0, _pv, strideV + sizeV); acc = BBE_MULRNX16(r0, _0x4000, q); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR + 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
        // final rotation
        _pr = (const xb_vecNx16 *)B;
    _pw = (xb_vecNx16 *)B;
    pFi = (const xb_vecNx16*)(V + BBE_SIMD_WIDTH*SIZE_OF_V(16, 8));
    strideR = -15 * 4 * L + 2 * BBE_SIMD_WIDTH;
    strideV = (SIZE_OF_V(16, 8) + 1) * 2 * BBE_SIMD_WIDTH;
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        xb_vecNx16 fi0, r0;
        xb_vecNx40 acc;
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_IP(fi0, pFi, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(r0, _pr, 4 * L); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, 4 * L);
        BBE_LVNX16_XP(fi0, pFi, strideV); BBE_LVNX16_XP(r0, _pr, strideR); acc = BBE_MULRNX16J(r0, fi0, q); r0 = BBE_PACKVNX40(acc, q); BBE_SVNX16_XP(r0, _pw, strideR);
    }
} /* cqr_calc_qb16x8x1s() */

size_t cqr_calc_qb16x8x1s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return 0;
} /* cqr_calc_qb16x8x1s_getScratchSize() */

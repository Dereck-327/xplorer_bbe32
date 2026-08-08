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

#define CURRENT_M 4
#define CURRENT_P 4

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
void cqr_calc_qb4x4x4s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    const int16_t * startV = V;
    int16_t * startB = B;
    xb_vecNx16 * restrict pv__;
    const xb_vecNx16 * restrict pRrd__;
    xb_vecNx16 * restrict pRwr__;
    int l;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);
    // rotate 2-nd column
    V = startV;
    B = startB;
    {
        int  j, rstep, vstep;
        vsaN _14 = BBE_MOVVSA32(14);
        xb_vecNx16 _0x4000 = BBE_MOVPINT16(16);
        const int v_back_step = -(CURRENT_M - 0 - 1) * 2 * BBE_SIMD_WIDTH;
        const int v_supr_step = v_back_step + 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
        const int R_last_step = 4 * L - CURRENT_P*L * 4 * (CURRENT_M - 0 - 1);
        const int R_supr_step = R_last_step - 16 * L + 2 * BBE_SIMD_WIDTH;
        xb_vecNx40 acc0;
        xb_vecNx16  v0, vR0;
        xb_vecNx16  ar[4];

        pv__ = (xb_vecNx16 *)V;
        pRrd__ = (xb_vecNx16 *)(B);
        pRwr__ = (xb_vecNx16 *)(B);
        for (j = 1, l = 0; l < L; l += (BBE_SIMD_WIDTH / 4), j ^= 1)
        {
            // note: inner loop is unrolled twice implicitely and rstep,vstep control right indexing
            rstep = R_last_step;
            vstep = v_back_step;
            XT_MOVEQZ(rstep, R_supr_step, j);
            XT_MOVEQZ(vstep, v_supr_step, j);
            /* Calculate vB = v__' * B(m__:end, :) */
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P * 4 * L);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[1], pRrd__, CURRENT_P * 4 * L);
            BBE_MULANX16J(acc0, ar[1], v0);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[2], pRrd__, CURRENT_P * 4 * L);
            BBE_MULANX16J(acc0, ar[2], v0);
            BBE_LVNX16_XP(v0, pv__, v_back_step);    BBE_LVNX16_XP(ar[3], pRrd__, R_last_step);
            BBE_MULANX16J(acc0, ar[3], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[2], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[2] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[2], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_XP(v0, pv__, v_back_step);
            acc0 = BBE_MULRNX16(ar[3], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[3] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[3], pRwr__, R_last_step);      /* go to next column of R__*/

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P * 4 * L);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[1], pRrd__, CURRENT_P * 4 * L);
            BBE_MULANX16J(acc0, ar[1], v0);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[2], pRrd__, CURRENT_P * 4 * L);
            BBE_MULANX16J(acc0, ar[2], v0);
            BBE_LVNX16_XP(v0, pv__, v_back_step);    BBE_LVNX16_XP(ar[3], pRrd__, rstep);
            BBE_MULANX16J(acc0, ar[3], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[2], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[2] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[2], pRwr__, CURRENT_P * 4 * L);
            BBE_LVNX16_XP(v0, pv__, vstep);
            acc0 = BBE_MULRNX16(ar[3], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[3] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[3], pRwr__, rstep);      /* go to next column of R__*/
        }
    }

    // rotate 3-nd column
    V = startV + BBE_SIMD_WIDTH*(CURRENT_M - 0);
    B = startB + CURRENT_P*L * 2 * 1;
    {
        vsaN _14 = BBE_MOVVSA32(14);
        xb_vecNx16 _0x4000 = BBE_MOVPINT16(16);
        int  j, rstep, vstep;
        const int v_back_step = -(CURRENT_M - 1 - 1)*BBE_SIMD_WIDTH*(int)sizeof(int16_t);
        const int R_last_step = -CURRENT_P*L * 4 * (CURRENT_M - 1 - 1) + L * 4;
        const int R_supr_step = R_last_step - 16 * L + 2 * BBE_SIMD_WIDTH;
        const int v_supr_step = v_back_step + 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
        xb_vecNx40 acc0;
        xb_vecNx16  v0, vR0;
        xb_vecNx16  ar[3];
        pRrd__ = (xb_vecNx16 *)(B);
        pRwr__ = (xb_vecNx16 *)(B);
        pv__ = (xb_vecNx16 *)V;
        for (j = 1, l = 0; l < L; l += (BBE_SIMD_WIDTH / 4), j ^= 1)
        {
            // note: inner loop is unrolled twice implicitely and rstep,vstep control right indexing
            rstep = R_last_step;
            vstep = v_back_step;
            XT_MOVEQZ(rstep, R_supr_step, j);
            XT_MOVEQZ(vstep, v_supr_step, j);

            /* Calculate vB = v__' * B(m__:end, :) */
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P*L * 4);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[1], pRrd__, CURRENT_P*L * 4);
            BBE_MULANX16J(acc0, ar[1], v0);
            BBE_LVNX16_XP(v0, pv__, v_back_step);     BBE_LVNX16_XP(ar[2], pRrd__, R_last_step);
            BBE_MULANX16J(acc0, ar[2], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_XP(v0, pv__, v_back_step);
            acc0 = BBE_MULRNX16(ar[2], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[2] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[2], pRwr__, R_last_step);       /* go to next column of R__*/

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P*L * 4);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[1], pRrd__, CURRENT_P*L * 4);
            BBE_MULANX16J(acc0, ar[1], v0);
            BBE_LVNX16_XP(v0, pv__, v_back_step);     BBE_LVNX16_XP(ar[2], pRrd__, rstep);
            BBE_MULANX16J(acc0, ar[2], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_XP(v0, pv__, vstep);
            acc0 = BBE_MULRNX16(ar[2], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[2] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[2], pRwr__, rstep);       /* go to next column of R__*/
        }
    }

    // rotate 4-th column
    V = startV + BBE_SIMD_WIDTH*(2 * CURRENT_M - 1);
    B = startB + CURRENT_P*L * 4;
    {
        vsaN _14 = BBE_MOVVSA32(14);
        xb_vecNx16 _0x4000 = BBE_MOVPINT16(16);
        int  j, rstep, vstep;
        const int v_back_step = -(CURRENT_M - 2 - 1)*BBE_SIMD_WIDTH * 2;
        const int R_last_step = -CURRENT_P*L * 4 * (CURRENT_M - 2 - 1) + L * 4;
        const int R_supr_step = R_last_step - 16 * L + 2 * BBE_SIMD_WIDTH;
        const int v_supr_step = v_back_step + 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
        xb_vecNx40 acc0;
        xb_vecNx16  v0, vR0;
        xb_vecNx16  ar[2];
        pRrd__ = (xb_vecNx16 *)(B);
        pRwr__ = (xb_vecNx16 *)(B);
        pv__ = (xb_vecNx16 *)V;
        for (j = 1, l = 0; l < L; l += (BBE_SIMD_WIDTH / 4), j ^= 1)
        {
            // note: inner loop is unrolled twice implicitely and rstep,vstep control right indexing
            rstep = R_last_step;
            vstep = v_back_step;
            XT_MOVEQZ(rstep, R_supr_step, j);
            XT_MOVEQZ(vstep, v_supr_step, j);
            /* Calculate vB = v__' * B(m__:end, :) */
            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P*L * 4);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_XP(v0, pv__, v_back_step);     BBE_LVNX16_XP(ar[1], pRrd__, R_last_step);
            BBE_MULANX16J(acc0, ar[1], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_XP(v0, pv__, v_back_step);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, R_last_step);/* go to next column of R__*/

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[0], pRrd__, CURRENT_P*L * 4);
            acc0 = BBE_MULRNX16J(ar[0], v0, _14);
            BBE_LVNX16_XP(v0, pv__, v_back_step);     BBE_LVNX16_XP(ar[1], pRrd__, rstep);
            BBE_MULANX16J(acc0, ar[1], v0);
            vR0 = BBE_PACKVNX40(acc0, _14);

            BBE_LVNX16_IP(v0, pv__, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16(ar[0], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[0] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[0], pRwr__, CURRENT_P*L * 4);
            BBE_LVNX16_XP(v0, pv__, vstep);
            acc0 = BBE_MULRNX16(ar[1], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[1] = BBE_PACKVNX40(acc0, _14);
            BBE_SVNX16_XP(ar[1], pRwr__, rstep);/* go to next column of R__*/
        }
    }


#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    // final rotation
    {
        vsaN _14 = BBE_MOVVSA32(14);

        int binc = 2 * BBE_SIMD_WIDTH - 15 * 4 * L;

        V = startV + BBE_SIMD_WIDTH*(3 * CURRENT_M - 3);
        pRrd__ = (const xb_vecNx16 *)(startB);
        pRwr__ = (xb_vecNx16 *)(startB);
        pv__ = (xb_vecNx16*)V;
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            xb_vecNx40 acc0;
            xb_vecNx16  fi, r;
            BBE_LVNX16_IP(fi, pv__, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);

            BBE_LVNX16_IP(fi, pv__, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);

            BBE_LVNX16_IP(fi, pv__, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);

            BBE_LVNX16_IP(fi, pv__, -3 * 2 * BBE_SIMD_WIDTH + 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M)));
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, 4 * L); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, 4 * L);
            BBE_LVNX16_XP(r, pRrd__, binc); acc0 = BBE_MULRNX16J(r, fi, _14); r = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_XP(r, pRwr__, binc);
        }
    }
} /* cqr_calc_qb4x4x4s() */

size_t cqr_calc_qb4x4x4s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return 0;
} /* cqr_calc_qb4x4x4s_getScratchSize() */

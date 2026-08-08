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
void  qr_calc_qbmxnxps (void *pScr, int16_t *B, const int16_t *V , int M, int N, int P, int L)
{
    vsaN q = BBE_MOVVSA32(14);
    xb_vecNx16 vr, _0x4000 = 0x4000, r0, v0;
    xb_vecNx40 acc;
    int Ncolumns = M > N ? N : N - 1;  // update N columns for overdetermined linear equation system
    int m, p, l;
    xb_vecNx16 * restrict _pw;
    const xb_vecNx16 * restrict _pr;
    const xb_vecNx16 * restrict _pv;

    NASSERT_ALIGN(B, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(V, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH) == 0);

    for (p = 0; p < P; p++)
    {
        int16_t *b = B;
        _pv = (const xb_vecNx16*)V;
        for (l = 0; l < L; l += BBE_SIMD_WIDTH)
        {
            int i;
            for (m = 0; m < Ncolumns; m++)
            {
                _pr = (const xb_vecNx16*)(b);
                _pw = (xb_vecNx16*)_pr;

                BBE_LVNX16_XP(r0, _pr, 2 * L*P);
                BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);

                acc = BBE_MULRNX16(r0, v0, q);
                for (i = 0; i < M - m - 1; i++)
                {
                    BBE_LVNX16_XP(r0, _pr, 2 * L*P);
                    BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);
                    BBE_MULANX16(acc, r0, v0);
                }
                vr = BBE_PACKVNX40(acc, q);

                _pr = (const xb_vecNx16*)_pw;
                _pv = (const xb_vecNx16*)XT_ADD(-(M - m) * 2 * BBE_SIMD_WIDTH, (uintptr_t)_pv);
                for (i = 0; i < M - m; i++)
                {
                    BBE_LVNX16_XP(r0, _pr, 2 * L*P);
                    BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);
                    acc = BBE_MULRNX16(r0, _0x4000, q);
                    BBE_MULSNX16(acc, vr, v0);
                    r0 = BBE_PACKVNX40(acc, q);
                    BBE_SVNX16_XP(r0, _pw, 2 * L*P);
                }
                b = (int16_t*)XT_ADD((uintptr_t)b, 2 * L*P);
            }
            b += -Ncolumns * 1 * L*P + BBE_SIMD_WIDTH;
        }
        B += L;
    }
} /* qr_calc_qbmxnxps() */

size_t  qr_calc_qbmxnxps_getScratchSize (int M, int P, int L)
{
    (void)M;
    (void)P;
    (void)L;
    return 0;
} /* qr_calc_qbmxnxps_getScratchSize() */

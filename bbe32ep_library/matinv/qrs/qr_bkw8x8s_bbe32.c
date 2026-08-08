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
    cqr_bkwNxPs/qr_bkwNxPs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if HAVE_VSAMATH && HAVE_DIV && 1

static void Solve8x8_real(
    const int16_t *R,
    int16_t *X,
    int L,
    const xb_vecNx16 * restrict invDiag
    )
{
    int k;
    const xb_vecNx16 * restrict  pR0;
    const xb_vecNx16 * restrict pR00;
    xb_vecNx16 * restrict pX;
    const xb_vecNx16 * restrict pB;

    const int R_last_step0 = -2 * 8 * L;
    const int R_last_step_inc = 2 * L;
    const int X_last_step = (8 * 8 - 1) * 2 * L;
    int R_last_step;
    pR00 = (const xb_vecNx16 *)(R + (8 * 8 - 8 - 1)*L);
    pX = (xb_vecNx16 *)(X + (8 * 8 - 1)*L);
    pB = (const xb_vecNx16 *)pX;
    vsaN q11 = BBE_MOVVSA32(11);
    xb_vecNx16 dp;
    xb_vecNx16 r0, ex, f0, b;
    xb_vecNx40 acc0;
    xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;

    for (k = 8 - 1; k >= 0; k--)
    {
        pR0 = pR00;
        R_last_step = R_last_step0;
        //m=7;
        {
            BBE_LVNX16_XP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(b, pB, -2 * L * 8);
            acc0 = BBE_MULRNX16(f0, b, ex);
            x7 = BBE_PACKVNX40(acc0, ex);
        }
        //m=6;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L * 8);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x6 = BBE_PACKVNX40(acc0, ex);
        }
        //m=5;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x6, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L * 8);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x5 = BBE_PACKVNX40(acc0, ex);
        }
        //m=4;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, 7 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x6, r0);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x5, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -5 * 2 * L * 8 + X_last_step);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x4 = BBE_PACKVNX40(acc0, ex);
        }
        BBE_SVNX16_XP(x7, pX, -2 * L * 8);
        BBE_SVNX16_XP(x6, pX, -2 * L * 8);
        BBE_SVNX16_XP(x5, pX, -2 * L * 8);
        BBE_SVNX16_XP(x4, pX, -5 * 2 * L * 8 + X_last_step);
    }
    BBE_LVNX16_XP(ex, invDiag, -8 * 2 * BBE_SIMD_WIDTH);
    __Pragma("no_reorder")
        // m=3...2
        pR00 = (const xb_vecNx16 *)(R + (8 * 8 - 8 - 1)*L);
    pB = (xb_vecNx16 *)(X + (8 * 8 - 1)*L);
    pX = (xb_vecNx16 *)XT_ADD(-4 * 2 * L * 8, (uintptr_t)pB);
    pR00 = (const xb_vecNx16*)XT_ADD(3 * R_last_step0 + 3 * R_last_step_inc - 3 * 2 * L, (uintptr_t)pR00);
    for (k = 8 - 1; k >= 0; k--)
    {
        pR0 = pR00;
        R_last_step = R_last_step0 + 3 * R_last_step_inc;

        BBE_LVNX16_XP(x7, pB, -2 * L * 8);
        BBE_LVNX16_XP(x6, pB, -2 * L * 8);
        BBE_LVNX16_XP(x5, pB, -2 * L * 8);
        BBE_LVNX16_XP(x4, pB, -2 * L * 8);
        //m=3;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, -2 * BBE_SIMD_WIDTH);
            acc0 = 0;
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x6, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x5, r0);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x4, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L * 8);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x3 = BBE_PACKVNX40(acc0, ex);
        }
        //m=2;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, 3 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x6, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x5, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x4, r0);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x3, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -3 * 2 * L * 8 + X_last_step);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x2 = BBE_PACKVNX40(acc0, ex);
        }
        BBE_SVNX16_XP(x3, pX, -2 * L * 8);
        BBE_SVNX16_XP(x2, pX, -7 * 2 * L * 8 + X_last_step);
    }
    BBE_LVNX16_XP(ex, invDiag, -4 * 2 * BBE_SIMD_WIDTH);
    __Pragma("no_reorder")

        // m=1...0
        pR00 = (const xb_vecNx16 *)(R + (8 * 8 - 8 - 1)*L);
    pB = (xb_vecNx16 *)(X + (8 * 8 - 1)*L);
    pX = (xb_vecNx16 *)XT_ADD(-6 * 2 * L * 8, (uintptr_t)pB);
    pR00 = (const xb_vecNx16*)XT_ADD(5 * R_last_step0 + 10 * R_last_step_inc - 10 * 2 * L, (uintptr_t)pR00);

    for (k = 8 - 1; k >= 0; k--)
    {
        pR0 = pR00;
        R_last_step = R_last_step0 + 5 * R_last_step_inc;
        BBE_LVNX16_XP(x7, pB, -2 * L * 8);
        BBE_LVNX16_XP(x6, pB, -2 * L * 8);
        BBE_LVNX16_XP(x5, pB, -2 * L * 8);
        BBE_LVNX16_XP(x4, pB, -2 * L * 8);
        BBE_LVNX16_XP(x3, pB, -2 * L * 8);
        BBE_LVNX16_XP(x2, pB, -2 * L * 8);
        //m=1;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x6, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x5, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x4, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x3, r0);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x2, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L * 8);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x1 = BBE_PACKVNX40(acc0, ex);
        }
        //m=0;
        {
            BBE_LVNX16_IP(ex, invDiag, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, invDiag, 3 * 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            acc0 = BBE_MULRNX16(x7, r0, q11);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x6, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x5, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x4, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x3, r0);
            BBE_LVNX16_XP(r0, pR0, -2 * L);
            BBE_MULANX16(acc0, x2, r0);
            BBE_LVNX16_XP(r0, pR0, R_last_step);
            BBE_MULANX16(acc0, x1, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L * 8 + X_last_step);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x0 = BBE_PACKVNX40(acc0, ex);
        }
        BBE_SVNX16_XP(x1, pX, -2 * L * 8);
        BBE_SVNX16_XP(x0, pX, -7 * 2 * L * 8 + X_last_step);
    }
}

/*-------------------------------------------------------------------------
cqr_bkwNxPs/qr_bkwNxPs

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Fixed-point representation for output data is a function of fixed-point
format of input data: FPP(X) = FPP(QB)-FPP(R)+10, where FPP(x) stands for
the Fixed-Point Position of data item x.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/

void  qr_bkw8x8s (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L)
{
    int i, k;
    xb_vecNx16 *invDiagR;
    const xb_vecNx16 *restrict pR = (const xb_vecNx16*)R;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH == 0);
    for (i = 0; i < L; i += BBE_SIMD_WIDTH)
    {
        xb_vecNx16 r0, f0, e0, _19 = 19;
        vsaN nsa;
        xb_vecNx40 _0x20000000 = BBE_MOVWA32(0x20000000);
        invDiagR = (xb_vecNx16*)pScr;
        pR = (const xb_vecNx16*)R;
        for (k = 0; k < 8; k++)
        {
            BBE_LVNX16_XP(r0, pR, (8 + 1) * 2 * L);
            nsa = BBE_NSANX16(r0);
            e0 = BBE_MOVVSV(nsa, 0);
            e0 = BBE_SUBNX16(_19, e0);
            r0 = BBE_SLANX16(r0, nsa);
            f0 = BBE_QUONX32(_0x20000000, r0);
            BBE_SVNX16_IP(f0, invDiagR, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(e0, invDiagR, 2 * BBE_SIMD_WIDTH);
        }
        Solve8x8_real(R, B, L, invDiagR - 1);
        R += BBE_SIMD_WIDTH;
        B += BBE_SIMD_WIDTH;
    }
} /* qr_bkw8x8s() */

size_t  qr_bkw8x8s_getScratchSize (int N, int P, int L)
{
    (void)P;
    (void)L;
    return 4 * N*BBE_SIMD_WIDTH;
} /* qr_bkw8x8s_getScratchSize() */
#else
DISCARD_FUN(void,qr_bkw8x8s,(void* pScr, int16_t* restrict B, const int16_t* restrict R, int L))
size_t  qr_bkw8x8s_getScratchSize (int N, int P, int L) { (void)N;(void)P;(void)L; return 0; }
#endif

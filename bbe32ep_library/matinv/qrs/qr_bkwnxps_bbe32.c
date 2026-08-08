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

inline_ void SolveMxPx8_real(
    const int16_t *R,
    int16_t *X,
    int N,
    int P,
    int L,
    const xb_vecNx16 *invDiag_end
    )
{
    int k, m, i;
    const xb_vecNx16 * restrict  pR0;
    const xb_vecNx16 * restrict pR00;
    xb_vecNx16 * pX;
    const xb_vecNx16 * restrict pB;
    const xb_vecNx16 * restrict pf;

    const int R_last_step0 = -2 * N*L;
    const int R_last_step_inc = 2 * L;
    const int X_last_step = (N*P - 1) * 2 * L;
    int R_last_step;
    pR00 = (const xb_vecNx16 *)(R + (N*N - N - 1)*L);
    pX = (xb_vecNx16 *)(X + (N*P - 1)*L);
    pB = (const xb_vecNx16 *)pX;
    vsaN q11 = BBE_MOVVSA32(11);

    NASSERT(N <= 16);
    for (k = P - 1; k >= 0; k--)
    {
        xb_vecNx16 dp;
        xb_vecNx16 r0, ex, f0, b, x;
        xb_vecNx40 acc0;
        xb_vecNx16 Xreg[16];

        pR0 = pR00;
        pf = (const xb_vecNx16  *)invDiag_end;
        R_last_step = R_last_step0;
        // separate case: m==M-1
        m = N - 1;
        {
            BBE_LVNX16_IP(ex, pf, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, pf, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(b, pB, -2 * L*P);
            acc0 = BBE_MULRNX16(f0, b, ex);
            x = BBE_PACKVNX40(acc0, ex);
            Xreg[m] = x;
            BBE_SVNX16_XP(x, pX, -2 * L*P);
        }
        for (m--; m >= 0; m--)
        {
            BBE_LVNX16_IP(ex, pf, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f0, pf, -2 * BBE_SIMD_WIDTH);
            acc0 = 0;
            for (i = N - 1; i > m + 1; i--)
            {
                BBE_LVNX16_XP(r0, pR0, -2 * L);
                BBE_MULANX16(acc0, Xreg[i], r0);
            }
            if (m != N - 1)
            {
                BBE_LVNX16_XP(r0, pR0, R_last_step);
                BBE_MULANX16(acc0, Xreg[i], r0);
                R_last_step += R_last_step_inc;
            }
            acc0 = BBE_RNDSADJNX40(acc0, q11);
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -2 * L*P);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x = BBE_PACKVNX40(acc0, ex);
            Xreg[m] = x;
            BBE_SVNX16_XP(x, pX, -2 * L*P);
        }
        pX = (xb_vecNx16*)XT_ADD(X_last_step, (uintptr_t)pX);
        pB = (xb_vecNx16*)XT_ADD(X_last_step, (uintptr_t)pB);
    } //for(k=P-1; k>=0; k--)
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

void  qr_bkwnxps (void* pScr, int16_t* restrict B, const int16_t* restrict R, int N, int P, int L)
{
    int i, k;
    xb_vecNx16 *invDiagR;
    const xb_vecNx16 *restrict pR = (const xb_vecNx16*)R;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT((L%BBE_SIMD_WIDTH) == 0 && L > 0);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);

    for (i = 0; i < L; i += BBE_SIMD_WIDTH)
    {
        xb_vecNx16 r0, f0, e0, _19 = 19;
        vsaN nsa;
        xb_vecNx40 _0x20000000 = BBE_MOVWA32(0x20000000);
        invDiagR = (xb_vecNx16*)pScr;
        pR = (const xb_vecNx16*)R;
        for (k = 0; k < N; k++)
        {
            BBE_LVNX16_XP(r0, pR, (N + 1) * 2 * L);
            nsa = BBE_NSANX16(r0);
            e0 = BBE_MOVVSV(nsa, 0);
            e0 = BBE_SUBNX16(_19, e0);
            r0 = BBE_SLANX16(r0, nsa);
            f0 = BBE_QUONX32(_0x20000000, r0);
            BBE_SVNX16_IP(f0, invDiagR, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(e0, invDiagR, 2 * BBE_SIMD_WIDTH);
        }
        SolveMxPx8_real(R, B, N, P, L, invDiagR - 1);
        R += BBE_SIMD_WIDTH;
        B += BBE_SIMD_WIDTH;
    }
} /* qr_bkwnxps() */

size_t  qr_bkwnxps_getScratchSize (int N, int P, int L)
{
    (void)P;
    (void)L;
    return 4 * N*BBE_SIMD_WIDTH;
} /* qr_bkwnxps_getScratchSize() */
#else
DISCARD_FUN(void,qr_bkwnxps,(void* pScr, int16_t* restrict B, const int16_t* restrict R, int N, int P, int L))
size_t  qr_bkwnxps_getScratchSize (int N, int P, int L) { (void)N; (void)P; (void)L; return 0; }
#endif

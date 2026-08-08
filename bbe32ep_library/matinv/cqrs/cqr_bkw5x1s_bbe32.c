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
/* Common utility declarations. */
#include "cqr_common.h"

#if HAVE_VSAMATH && 1

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

void cqr_bkw5x1s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    int l, i, k;
    xb_vecNx16 *restrict invDiagR;
    const xb_vecNx16 *restrict pR;
    const xb_vecNx16 * restrict pR00;
    xb_vecNx16 * restrict pX;
    const xb_vecNx16 * restrict pB;
    xb_vecNx16 dp, r0, ex, f0, b;
    xb_vecNx40 acc0;
    xb_vecNx16 x0, x1, x2, x3, x4;
    vsaN q11 = BBE_MOVVSA32(11);
    int R_last_step0, R_last_step_inc, R_last_step;

    NASSERT_ALIGN(B, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    NASSERT((L % (BBE_SIMD_WIDTH / 2)) == 0 && L > 0);
    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    // compute reciprocal of the main diagonal
    invDiagR = (xb_vecNx16*)pScr;
    {
        int stride;
        xb_vecNx16 r0, f0, e0, _19 = 19;
        vsaN nsa;
        xb_vecNx40 _0x20000000 = BBE_MOVWA32(0x20000000);
        pR = (const xb_vecNx16*)(R);
        for (i = k = 0; k < 5 * (L / (BBE_SIMD_WIDTH / 2)); k++)
        {
            i = BBE_ADDMOD16U(i, (5 << 16) | 1);
            stride = (5 + 1) * 4 * L;
            XT_MOVEQZ(stride, -(5 - 1)*(5 + 1) * 4 * L + 2 * BBE_SIMD_WIDTH, i);
            BBE_LVNX16_XP(r0, pR, stride);
            nsa = BBE_NSANX16(r0);
            e0 = BBE_MOVVSV(nsa, 0);
            e0 = BBE_SUBNX16(_19, e0);
            r0 = BBE_SLANX16(r0, nsa);
            QUO8X32(f0, _0x20000000, r0);
            f0 = BBE_SELNX16I(e0, f0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
            BBE_SVNX16_IP(f0, invDiagR, 2 * BBE_SIMD_WIDTH);
        }
    }
    // back sustitution
    R_last_step0 = -4 * 5 * L;
    R_last_step_inc = 4 * L;

    invDiagR = (xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH * 5, (uintptr_t)pScr);
    invDiagR--;
    pX = (xb_vecNx16 *)(B + (5 - 1)*L * 2);
    pB = (const xb_vecNx16 *)pX;
    pR00 = (const xb_vecNx16 *)(R + (5 * 5 - 5 - 1)*L * 2);
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        pR = pR00;
        //m=4;
        {
            BBE_LVNX16_IP(f0, invDiagR, -2 * BBE_SIMD_WIDTH);
            ex = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP(b, pB, -4 * L);
            acc0 = BBE_MULRNX16(f0, b, ex);
            x4 = BBE_PACKVNX40(acc0, ex);
        }
        //m=3;
        {
            BBE_LVNX16_IP(f0, invDiagR, 6 * 2 * BBE_SIMD_WIDTH);
            ex = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP(r0, pR, -18 * 4 * L);
            acc0 = BBE_MULRNX16C(x4, r0, q11);
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, 4 * L * 1 + 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x3 = BBE_PACKVNX40(acc0, ex);
        }

        BBE_SVNX16_XP(x4, pX, -4 * L);
        BBE_SVNX16_XP(x3, pX, 4 * L * 1 + 2 * BBE_SIMD_WIDTH);
        pR00++;
    }
    __Pragma("no_reorder")
        invDiagR = (xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH * 3, (uintptr_t)pScr);
    invDiagR--;
    pB = (const xb_vecNx16 *)(B + (5 - 1)*L * 2);
    pX = (xb_vecNx16 *)pB;
    pR00 = (const xb_vecNx16 *)(R + (5 * 5 - 5 - 1)*L * 2);
    pR00 = (const xb_vecNx16 *)XT_ADD(R_last_step0, (uintptr_t)pR00);
    pX = (xb_vecNx16 *)XT_ADD(-4 * L * 2, (uintptr_t)pX);
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        pR = pR00;
        R_last_step = -4 * 4 * L;
        BBE_LVNX16_XP(x4, pB, -4 * L);
        BBE_LVNX16_XP(x3, pB, -4 * L);
        //m=2;
        {
            BBE_LVNX16_IP(f0, invDiagR, -2 * BBE_SIMD_WIDTH);
            ex = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            acc0 = BBE_MULRNX16C(x4, r0, q11);
            BBE_LVNX16_XP(r0, pR, R_last_step);
            BBE_MULANX16C(acc0, x3, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -4 * L);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x2 = BBE_PACKVNX40(acc0, ex);
        }
        //m=1;
        {
            BBE_LVNX16_IP(f0, invDiagR, -2 * BBE_SIMD_WIDTH);
            ex = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            acc0 = BBE_MULRNX16C(x4, r0, q11);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            BBE_MULANX16C(acc0, x3, r0);
            BBE_LVNX16_XP(r0, pR, R_last_step);
            BBE_MULANX16C(acc0, x2, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, -4 * L);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x1 = BBE_PACKVNX40(acc0, ex);
        }
        //m=0;
        {
            BBE_LVNX16_IP(f0, invDiagR, 7 * 2 * BBE_SIMD_WIDTH);
            ex = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            acc0 = BBE_MULRNX16C(x4, r0, q11);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            BBE_MULANX16C(acc0, x3, r0);
            BBE_LVNX16_XP(r0, pR, -4 * L);
            BBE_MULANX16C(acc0, x2, r0);
            BBE_LVNX16_XP(r0, pR, R_last_step);
            BBE_MULANX16C(acc0, x1, r0);
            R_last_step += R_last_step_inc;
            dp = BBE_PACKVNX40(acc0, q11);
            BBE_LVNX16_XP(b, pB, 4 * L * 4 + 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            ex = BBE_ADDSAVSN(-1, ex);
            x0 = BBE_PACKVNX40(acc0, ex);
        }
        BBE_SVNX16_XP(x2, pX, -4 * L);
        BBE_SVNX16_XP(x1, pX, -4 * L);
        BBE_SVNX16_XP(x0, pX, 4 * L * 2 + 2 * BBE_SIMD_WIDTH);
        pR00++;
    }
} /* cqr_bkw5x1s() */

size_t cqr_bkw5x1s_getScratchSize (int N, int P, int L)
{
    (void)P;
    return 4 * N*L;
} /* cqr_bkw5x1s_getScratchSize() */
#else
DISCARD_FUN(void,cqr_bkw5x1s,(void* pScr, complex_fract16* restrict B,const complex_fract16* restrict R, int L))
size_t cqr_bkw5x1s_getScratchSize (int N, int P, int L) { (void)N;(void)P;(void)L; return 0; }
#endif

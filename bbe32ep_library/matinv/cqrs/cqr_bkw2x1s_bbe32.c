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

void cqr_bkw2x1s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    int l;
    const xb_vecNx16 *restrict pR;
    xb_vecNx16 *restrict pX;
    const xb_vecNx16 *restrict pB;

    xb_vecNx16 dp, r10, r00, r11, ex1, ex0, f1, f0, b, _19 = 19;
    xb_vecNx40 acc0;
    xb_vecNx16 x0, x1;
    vsaN q11 = BBE_MOVVSA32(11), nsa;
    xb_vecNx40 _0x20000000 = BBE_MOVWA32(0x20000000);

    NASSERT_ALIGN(B, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    NASSERT((L % (BBE_SIMD_WIDTH / 2)) == 0 );
    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    (void)pScr;
    if (L<=0) return;
    // compute reciprocal of the main diagonal
    // back sustitution
    pX = (xb_vecNx16 *)(B);
    pB = (const xb_vecNx16 *)pX;
    pR = (const xb_vecNx16 *)(R);

    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        BBE_LVNX16_IP(r00, pR, 2 * BBE_SIMD_WIDTH);
        r10 = BBE_LVNX16_X(pR, 4 * L - 2 * BBE_SIMD_WIDTH);
        r11 = BBE_LVNX16_X(pR, 3 * 4 * L - 2 * BBE_SIMD_WIDTH);

        nsa = BBE_NSANX16(r00);
        r00 = BBE_SLANX16(r00, nsa);
        ex0 = BBE_MOVVSV(nsa, 0);

        nsa = BBE_NSANX16(r11);
        r11 = BBE_SLANX16(r11, nsa);
        ex1 = BBE_MOVVSV(nsa, 0);

        ex0 = BBE_SUBNX16(_19, ex0);
        ex1 = BBE_SUBNX16(_19, ex1);
        QUO8X32(f1, _0x20000000, r11);
        QUO8X32(f0, _0x20000000, r00);
        ex1 = BBE_SHFLNX16I(ex1, BBE_SHFLI_DUPLICATE_1_EVEN);
        ex0 = BBE_SHFLNX16I(ex0, BBE_SHFLI_DUPLICATE_1_EVEN);
        f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DUPLICATE_1_EVEN);
        f1 = BBE_SHFLNX16I(f1, BBE_SHFLI_DUPLICATE_1_EVEN);
        b = BBE_LVNX16_X(pB, 4 * L);
        acc0 = BBE_MULRNX16(f1, b, ex1);
        x1 = BBE_PACKVNX40(acc0, ex1);

        acc0 = BBE_MULRNX16C(x1, r10, q11);
        dp = BBE_PACKVNX40(acc0, q11);
        BBE_LVNX16_IP(b, pB, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        ex0 = BBE_ADDSAVSN(-1, ex0);
        x0 = BBE_PACKVNX40(acc0, ex0);

        BBE_SVNX16_X(x1, pX, 4 * L);
        BBE_SVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH); 
    }
} /* cqr_bkw2x1s() */
#else
DISCARD_FUN(void,cqr_bkw2x1s,(void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L))
#endif

size_t cqr_bkw2x1s_getScratchSize (int N, int P, int L)
{
    (void)P; (void)N; (void)L;
    return 0;
} /* cqr_bkw2x1s_getScratchSize() */

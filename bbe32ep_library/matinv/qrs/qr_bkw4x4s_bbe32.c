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

#define UPDATEB(b,f0,shift)          \
{                                    \
    xb_vecNx40 acc0;                 \
    acc0 = BBE_MULRNX16(f0, b,shift);\
    BBE_MULSNX16(acc0,f0,dp);        \
    b = BBE_PACKVNX40(acc0, shift);  \
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

void  qr_bkw4x4s (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L)
{
    const int q = 10;
    int i;
    int16_t *QB = B;

    const int16_t *start_R = R;
    xb_vecNx16 * restrict  pB0;
    xb_vecNx16 * restrict  pB1;
    xb_vecNx16 * restrict  pB2;
    xb_vecNx16 * restrict  pB3;

    const xb_vecNx16 * restrict  pr12;
    const xb_vecNx16 * restrict  pr23;
    const xb_vecNx16 * restrict  pr34;
    const xb_vecNx16 * restrict  pr44;
    int _0x2000000 = 0x2000000;
    xb_vecNx40 z0;//  = BBE_MOVWA32(0x2000000);
    xb_vecNx16 z4 = BBE_MOVVINT16(4);
    xb_vecNx16 c19 = BBE_MOVVA16(30 - 4 - q + 3);
    xb_vecNx16 nsa;
    vsaN nsa_vsa;
#define SHRB(b) b=BBE_SRAINX16(b,1);

#ifdef COMPILER_XTENSA
#pragma ymemory( pr12 )
#pragma ymemory( pr23 )
#pragma ymemory( pr34 )
#endif

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT((L%BBE_SIMD_WIDTH) == 0 && L > 0);
    (void)pScr;

    pr44 = (xb_vecNx16 *)(R + 3 * (4 + 1)*L * 1);
    R = start_R;
    {

        const int start_m = 4 - 1;

        pr34 = (xb_vecNx16 *)(R + (start_m - 1) * 4 * L * 1/*row = start_m-1*/ + start_m*L * 1 /*column = start_m*/);
        pr23 = pr34 - (1 * 4 + 1)*L * 2 / sizeof(*pr34);
        pr12 = pr34 - (2 * 4 + 2)*L * 2 / sizeof(*pr34);

        pB0 = (xb_vecNx16 *)(QB + start_m * 4 * L * 1 + 0 * L * 1);
        pB1 = (xb_vecNx16 *)(QB + start_m * 4 * L * 1 + 1 * L * 1);
        pB2 = (xb_vecNx16 *)(QB + start_m * 4 * L * 1 + 2 * L * 1);
        pB3 = (xb_vecNx16 *)(QB + start_m * 4 * L * 1 + 3 * L * 1);

        for (i = 0; i < L; i += BBE_SIMD_WIDTH)
        {
            vsaN _11 = BBE_MOVVSA32(11);
            xb_vecNx16 ex;
            xb_vecNx16 dp, b, b3;
            xb_vecNx40 acc0;
            xb_vecNx40 acc00;
            xb_vecNx16 inv_r, inv_r00;
            vsaN       inv_re, inv_r11e;
            xb_vecNx16 r34, r14, r13, r12, r24, r23;

            /**************cols 0...3 row 3 *********************/
            BBE_LVNX16_XP(inv_r, pr44, -(4 + 1)*L * 2); //
            nsa_vsa = BBE_NSANX16(inv_r);
            nsa = BBE_MOVVVS(nsa_vsa);
            ex = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa, z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r = BBE_SLANX16(inv_r, nsa_vsa);
            z0 = BBE_MOVWA32(_0x2000000);
            XT_MOVEQZ(_0x2000000, _0x2000000, _0x2000000);
            inv_r = BBE_QUONX32(z0, inv_r);
            inv_re = BBE_MOVVSV(ex, 0);
            b = BBE_LVNX16_I(pB0, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB0, -4 * L * 2);
            b = BBE_LVNX16_I(pB1, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB1, -4 * L * 2);
            b = BBE_LVNX16_I(pB2, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB2, -4 * L * 2);
            b = BBE_LVNX16_I(pB3, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB3, -4 * L * 2);
            /************** cols 0...3 row 2 *********************/
            BBE_LVNX16_XP(inv_r, pr44, -(4 + 1)*L * 2);
            nsa_vsa = BBE_NSANX16(inv_r);
            nsa = BBE_MOVVVS(nsa_vsa);
            ex = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa, z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r = BBE_SLANX16(inv_r, nsa_vsa);
            z0 = BBE_MOVWA32(_0x2000000);
            XT_MOVEQZ(_0x2000000, _0x2000000, _0x2000000);
            inv_r = BBE_QUONX32(z0, inv_r);
            inv_re = BBE_MOVVSV(ex, 0);
            inv_re = BBE_ADDSAVSN(-1, inv_re);
            BBE_LVNX16_XP(r34, pr34, -4 * L * 2);
            BBE_LVNX16_XP(r24, pr34, -4 * L * 2);
            BBE_LVNX16_XP(r14, pr34, 2 * BBE_SIMD_WIDTH + 2 * 4 * L * 2);
            b3 = BBE_LVNX16_X(pB0, 4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r34, _11);
            dp = BBE_PACKVNX40(acc0, _11);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_re);
            BBE_SVNX16_IP(b, pB0, 0);
            b3 = BBE_LVNX16_X(pB1, 4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r34, _11);
            dp = BBE_PACKVNX40(acc0, _11);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_re);
            BBE_SVNX16_IP(b, pB1, 0);
            b3 = BBE_LVNX16_X(pB2, 4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r34, _11);
            dp = BBE_PACKVNX40(acc0, _11);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_re);
            BBE_SVNX16_IP(b, pB2, 0);
            b3 = BBE_LVNX16_X(pB3, 4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r34, _11);
            dp = BBE_PACKVNX40(acc0, _11);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_re);
            BBE_SVNX16_IP(b, pB3, 0);

            /************** col 0 row 1 *********************/
            BBE_LVNX16_XP(inv_r, pr44, -(4 + 1)*L * 2);
            nsa_vsa = BBE_NSANX16(inv_r);
            nsa = BBE_MOVVVS(nsa_vsa);
            ex = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa, z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r = BBE_SLANX16(inv_r, nsa_vsa);
            z0 = BBE_MOVWA32(_0x2000000);
            XT_MOVEQZ(_0x2000000, _0x2000000, _0x2000000);
            inv_r = BBE_QUONX32(z0, inv_r);

            inv_r11e = BBE_MOVVSV(ex, 0);
            inv_r11e = BBE_ADDSAVSN(-1, inv_r11e);
            BBE_LVNX16_XP(inv_r00, pr44, 2 * BBE_SIMD_WIDTH + (4 - 1)*(4 + 1)*L * 2);
            nsa_vsa = BBE_NSANX16(inv_r00);
            nsa = BBE_MOVVVS(nsa_vsa);
            ex = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa, z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r00 = BBE_SLANX16(inv_r00, nsa_vsa);
            inv_r00 = BBE_QUONX32(z0, inv_r00);

            inv_re = BBE_MOVVSV(ex, 0);
            inv_re = BBE_ADDSAVSN(-1, inv_re);
            BBE_LVNX16_XP(r23, pr23, -4 * L * 2);
            BBE_LVNX16_XP(r13, pr23, 2 * BBE_SIMD_WIDTH + 4 * L * 2);
            BBE_LVNX16_IP(r12, pr12, 2 * BBE_SIMD_WIDTH);

            b3 = BBE_LVNX16_X(pB0, 4 * L * 2);
            BBE_LVNX16_XP(b, pB0, -4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r24, _11);
            BBE_MULANX16(acc0, b, r23);
            dp = BBE_PACKVNX40(acc0, _11);
            acc00 = BBE_MULRNX16(b3, r14, _11);
            BBE_MULANX16(acc00, b, r13);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_r11e);
            BBE_SVNX16_XP(b, pB0, -4 * L * 2);
            /************** col 0 row 0 *********************/
            BBE_MULANX16(acc00, b, r12);
            dp = BBE_PACKVNX40(acc00, _11);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b, inv_r00, inv_re);
            BBE_SVNX16_XP(b, pB0, 2 * BBE_SIMD_WIDTH + (4 - 1) * 4 * L * 2);
            /************** col 1 row 1 *********************/
            b3 = BBE_LVNX16_X(pB1, 4 * L * 2);
            BBE_LVNX16_XP(b, pB1, -4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r24, _11);
            BBE_MULANX16(acc0, b, r23);
            dp = BBE_PACKVNX40(acc0, _11);
            acc00 = BBE_MULRNX16(b3, r14, _11);
            BBE_MULANX16(acc00, b, r13);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_r11e);
            BBE_SVNX16_XP(b, pB1, -4 * L * 2);
            /************** col 1 row 0 *********************/
            BBE_MULANX16(acc00, b, r12);
            dp = BBE_PACKVNX40(acc00, _11);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b, inv_r00, inv_re);
            BBE_SVNX16_XP(b, pB1, 2 * BBE_SIMD_WIDTH + (4 - 1) * 4 * L * 2);
            /************** col 2 row 1 *********************/
            b3 = BBE_LVNX16_X(pB2, 4 * L * 2);
            BBE_LVNX16_XP(b, pB2, -4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r24, _11);
            BBE_MULANX16(acc0, b, r23);
            dp = BBE_PACKVNX40(acc0, _11);
            acc00 = BBE_MULRNX16(b3, r14, _11);
            BBE_MULANX16(acc00, b, r13);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_r11e);
            BBE_SVNX16_XP(b, pB2, -4 * L * 2);
            /************** col 2 row 0 *********************/
            BBE_MULANX16(acc00, b, r12);
            dp = BBE_PACKVNX40(acc00, _11);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b, inv_r00, inv_re);
            BBE_SVNX16_XP(b, pB2, 2 * BBE_SIMD_WIDTH + (4 - 1) * 4 * L * 2);
            /************** col 3 row 1 *********************/
            b3 = BBE_LVNX16_X(pB3, 4 * L * 2);
            BBE_LVNX16_XP(b, pB3, -4 * L * 2);
            acc0 = BBE_MULRNX16(b3, r24, _11);
            BBE_MULANX16(acc0, b, r23);
            dp = BBE_PACKVNX40(acc0, _11);
            acc00 = BBE_MULRNX16(b3, r14, _11);
            BBE_MULANX16(acc00, b, r13);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b, inv_r, inv_r11e);
            BBE_SVNX16_XP(b, pB3, -4 * L * 2);
            /************** col 3 row 0 *********************/
            BBE_MULANX16(acc00, b, r12);
            dp = BBE_PACKVNX40(acc00, _11);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b, inv_r00, inv_re);
            BBE_SVNX16_XP(b, pB3, 2 * BBE_SIMD_WIDTH + (4 - 1) * 4 * L * 2);
        }
    }
} /* qr_bkw4x4s() */
#else
DISCARD_FUN(void,qr_bkw4x4s,(void* pScr, int16_t* restrict B, const int16_t* restrict R, int L))
#endif

size_t  qr_bkw4x4s_getScratchSize (int N, int P, int L)
{
    (void)N;(void)P;(void)L;
    return 0;
} /* qr_bkw4x4s_getScratchSize() */

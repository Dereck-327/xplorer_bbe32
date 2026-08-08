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
#include "qr_common.h"

#if !(HAVE_VSAMATH && HAVE_DIV && 1)
DISCARD_FUN(void, qr_bkw8x1s, (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L))
#else

#define CURRENT_M 8
#define CURRENT_P 1

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

void  qr_bkw8x1s (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L)
{
    xb_vecNx16 * restrict pb0 = (xb_vecNx16 *)(B);
    xb_vecNx16 * restrict pb4 = (xb_vecNx16 *)(B + 4 * (1 * L));
    const xb_vecNx16 * restrict pR27 = (const xb_vecNx16 *)(R + (2 * CURRENT_M + 7)*L);
    const xb_vecNx16 * restrict pR67 = (const xb_vecNx16 *)(R + (6 * CURRENT_M + 7)*L);
    const xb_vecNx16 * restrict pr77 = (const xb_vecNx16 *)(R + (7 * CURRENT_M + 7)*L);

#ifdef COMPILER_XTENSA
#pragma ymemory( pr77 )
#pragma ymemory( pR27 ) 
#pragma ymemory( pR67 )
#endif

    int i;

    xb_vecNx16
        r01, r02, r03, r04, r05, r06, r07,
        r12, r13, r14, r15, r16, r17,
        r23, r24, r25, r26, r27,
        r34, r35, r36, r37,
        r45, r46, r47,
        r56, r57,
        r67;
    const int o1 = 1 * 1 * L*sizeof(int16_t);
    const int o2 = 2 * 1 * L*sizeof(int16_t);
    const int o3 = 3 * 1 * L*sizeof(int16_t);
    const int stepR = -(1 * L*(CURRENT_M + 1))*sizeof(int16_t);
    const int last_step_pR67 = sizeof(*pR67) - stepR * 6;
    const int last_step_pR27 = sizeof(*pR27) - stepR * 2;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH == 0);
    (void)pScr;

    for (i = 0; i < L; i += BBE_SIMD_WIDTH)
    {
        const xb_vecNx40 z0 = BBE_MOVWA32(0x20000000);
        xb_vecNx16 x0, f0;
        vsaN e;

        xb_vecNx16 dp;
        xb_vecNx16  b;
        xb_vecNx40 acc0;
        xb_vecNx16 b7, b6, b5, b4, b3, b2, b1, b0;

        const int q = 10;
        vsaN sh11 = BBE_MOVVSA32(q + 1);
        vsaN ex_vsaN;

        //m=7;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        dp = BBE_ZERONX16();

        b = BBE_LVNX16_X(pb4, o3);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b7 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b7, pb4, o3);
        //}
        //m=6;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        BBE_LVNX16_XP(r67, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r67, sh11);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_X(pb4, o2);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b6 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b6, pb4, o2);
        //}
        //m=5;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        r57 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r56, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r57, sh11);
        BBE_MULANX16(acc0, b6, r56);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_X(pb4, o1);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b5 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b5, pb4, o1);
        //}
        //m=4;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        r47 = BBE_LVNX16_X(pR67, o2);
        r46 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r45, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r47, sh11);
        BBE_MULANX16(acc0, b6, r46);
        BBE_MULANX16(acc0, b5, r45);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_I(pb4, 0);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b4 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_IP(b4, pb4, sizeof(*pb4));
        //}
        //m=3;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        r37 = BBE_LVNX16_X(pR67, o3);
        r36 = BBE_LVNX16_X(pR67, o2);
        r35 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r34, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r37, sh11);
        BBE_MULANX16(acc0, b6, r36);
        BBE_MULANX16(acc0, b5, r35);
        BBE_MULANX16(acc0, b4, r34);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_X(pb0, o3);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b3 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b3, pb0, o3);
        //}
        //m=2;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        BBE_LVNX16_XP(r27, pR27, stepR);
        r26 = BBE_LVNX16_X(pR67, o3);
        r25 = BBE_LVNX16_X(pR67, o2);
        r24 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r23, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r27, sh11);
        BBE_MULANX16(acc0, b6, r26);
        BBE_MULANX16(acc0, b5, r25);
        BBE_MULANX16(acc0, b4, r24);
        BBE_MULANX16(acc0, b3, r23);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_X(pb0, o2);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b2 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b2, pb0, o2);
        //}
        //m=1;
        //{
        BBE_LVNX16_XP(x0, pr77, -(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        r17 = BBE_LVNX16_X(pR27, o1);
        BBE_LVNX16_XP(r16, pR27, stepR);
        r15 = BBE_LVNX16_X(pR67, o3);
        r14 = BBE_LVNX16_X(pR67, o2);
        r13 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r12, pR67, stepR);
        acc0 = BBE_MULRNX16(b7, r17, sh11);
        BBE_MULANX16(acc0, b6, r16);
        BBE_MULANX16(acc0, b5, r15);
        BBE_MULANX16(acc0, b4, r14);
        BBE_MULANX16(acc0, b3, r13);
        BBE_MULANX16(acc0, b2, r12);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_X(pb0, o1);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b1 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_X(b1, pb0, o1);
        //}
        //m=0;
        //{
        BBE_LVNX16_XP(x0, pr77, 2 * BBE_SIMD_WIDTH + (CURRENT_M - 1)*(CURRENT_M + 1)*L*sizeof(int16_t));
        e = BBE_NSANX16(x0);
        x0 = BBE_SLANX16(x0, e);
        f0 = BBE_QUONX32(z0, x0);
        ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

        r07 = BBE_LVNX16_X(pR27, o2);
        r06 = BBE_LVNX16_X(pR27, o1);
        BBE_LVNX16_XP(r05, pR27, last_step_pR27);
        r04 = BBE_LVNX16_X(pR67, o3);
        r03 = BBE_LVNX16_X(pR67, o2);
        r02 = BBE_LVNX16_X(pR67, o1);
        BBE_LVNX16_XP(r01, pR67, last_step_pR67);
        acc0 = BBE_MULRNX16(b7, r07, sh11);
        BBE_MULANX16(acc0, b6, r06);
        BBE_MULANX16(acc0, b5, r05);
        BBE_MULANX16(acc0, b4, r04);
        BBE_MULANX16(acc0, b3, r03);
        BBE_MULANX16(acc0, b2, r02);
        BBE_MULANX16(acc0, b1, r01);

        dp = BBE_PACKVNX40(acc0, sh11);
        b = BBE_LVNX16_I(pb0, 0);
        acc0 = BBE_MULNX16(f0, b);
        acc0 = BBE_SRAINX40(acc0, 1);
        BBE_MULSNX16(acc0, f0, dp);
        acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
        b0 = BBE_PACKVNX40(acc0, ex_vsaN);
        BBE_SVNX16_IP(b0, pb0, sizeof(*pb0));
        //}
    }
    //return qr_bkwnxps(pScr, B, R, 8, 1, L); 
} /* qr_bkw8x1s() */
#endif

size_t  qr_bkw8x1s_getScratchSize (int N, int P, int L)
{
    (void)N;(void)P;(void)L;
    return 0;
} /* qr_bkw8x1s_getScratchSize() */

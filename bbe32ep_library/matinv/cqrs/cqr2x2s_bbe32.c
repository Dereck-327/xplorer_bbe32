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
    These functions apply QR decomposition procedure to the sequence of complex
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if (HAVE_VSAMATH && HAVE_NSAENX40 && 1)

#define CURRENT_M 2

inline_ void chouseholder_2x16(xb_vecNx16 *pr0,
                               xb_vecNx16 *pr1,
                               xb_vecNx16 *pv0,
                               xb_vecNx16 *pv1,
                               xb_vecNx16 *pfi)
{
    xb_vecNx16 tmp0, tmp1;
    vsaN c_vec, v_exp;
    xb_vecNx40 acc0;
    xb_vecNx16 k_norm, mant;
    xb_vecNx16 norm_x, Fi;
    xb_vecNx16  x1, x2;

    x1 = *pr0;

    acc0 = BBE_MAGINX16C(x1, x1);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    acc0 = BBE_SRAINX40(acc0, 24);
    tmp0 = BBE_PACKLNX40(acc0);
    tmp0 = BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
    c_vec = BBE_SUBSR1SAVSN(18 + 1, c_vec);
    acc0 = BBE_MULUSRNX16(tmp0, x1, c_vec);
    norm_x = BBE_PACKVNX40(acc0, c_vec);

    Fi = norm_x;
    *pfi = norm_x;
    /* Build Householder's vector   v = x/sqrt(x'*x) + e1*x(1)/abs(x(1) */
    x2 = *pr1;

    acc0 = BBE_MAGINX16C(x1, x1);
    BBE_MAGIANX16C(acc0, x2, x2);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    acc0 = BBE_SRAINX40(acc0, 24);
    tmp0 = BBE_PACKLNX40(acc0);
    mant = BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
    v_exp = BBE_SUBSR1SAVSN(18 + 1, c_vec);
    acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
    norm_x = BBE_PACKVNX40(acc0, v_exp);

    x1 = *pr0;
    norm_x = BBE_ADDNX16(norm_x, Fi);
    acc0 = BBE_MULNX16J(norm_x, Fi);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    c_vec = BBE_SUBSR1SAVSN(28, c_vec);
    k_norm = BBE_PACKVNX40(acc0, c_vec);
    k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);

    acc0 = BBE_MULNX16(norm_x, k_norm);
    norm_x = BBE_PACKQNX40(acc0);

    *pv0 = norm_x;
    x1 = *pr1;
    acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
    norm_x = BBE_PACKVNX40(acc0, v_exp);
    norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
    *pv1 = norm_x;
}

/*-------------------------------------------------------------------------
These functions apply QR decomposition procedure to the sequence of complex 
matrices written in a streaming order. The transformation is done in-place 
so the result replaces the original input.
Rotation matrix Q is calculated in Q15 fixed point representation. Fixed 
point representation of upper-diagonal matrix R is the same as of input. 

Functions return nonzero if overflow is detected 

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L][C]   input matrices
C              1 for real, 2 for complex data
Output:
Q[M*M][L][C]   output rotation matrices (L matrices of size MxM)
R[M*N][L][C]   output upper triangular matrices (L matrices of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Matrix sizes M,N,L must be greater than 1
4. Scratch memory must be aligned on 32-byte boundary. Its size  (in bytes)
   is defined by cqr2x2s_getScratchSize, qr2x2s_getScratchSize
---------------------------------------------------------------------------*/
void cqr2x2s (void* pScr, complex_fract16 * restrict _Q, complex_fract16 * restrict _R, int L)
{
    int16_t * restrict Q=(int16_t *)_Q;
    int16_t * restrict R=(int16_t *)_R;
    int l;

    const int q = 14;
    xb_vecNx16  vR0, vR1, ar[CURRENT_M*CURRENT_M];
    xb_vecNx40 acc0;

    const xb_vecNx16 * restrict pRRr;
    xb_vecNx16 * restrict pRRw;
    xb_vecNx16 * restrict pV = (xb_vecNx16*)pScr;

#ifdef COMPILER_XTENSA

#pragma ymemory( pRRr )
#pragma ymemory( pRRw )

#endif

    xb_vecNx16 * restrict pQ00 = (xb_vecNx16 *)(Q + 0 * L * 2);
    xb_vecNx16 * restrict pQ01 = (xb_vecNx16 *)(Q + 1 * L * 2);
    xb_vecNx16 * restrict pQ10 = (xb_vecNx16 *)(Q + 2 * L * 2);
    xb_vecNx16 * restrict pQ11 = (xb_vecNx16 *)(Q + 3 * L * 2);
    xb_vecNx16 one = BBE_MOVVA16C(1 << q);

    NASSERT_ALIGN(pScr, (BBE_SIMD_WIDTH)* 2);
    NASSERT_ALIGN(Q, (BBE_SIMD_WIDTH)* 2);
    NASSERT_ALIGN(R, (BBE_SIMD_WIDTH)* 2);
    NASSERT(L > 0 && L % ((BBE_SIMD_WIDTH) / 2) == 0);

    pRRr = (xb_vecNx16 *)(R);
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        xb_vecNx16 av[CURRENT_M];
        xb_vecNx16 afi[CURRENT_M];

        ar[1] = BBE_LVNX16_X(pRRr, 2 * 4 * L);
        BBE_LVNX16_IP(ar[0], pRRr, 2 * BBE_SIMD_WIDTH);

        chouseholder_2x16(&ar[0], &ar[1],
            &av[0], &av[1],
            &afi[0]);

        afi[0] = BBE_SLSINX16(afi[0], 1);  // Convert to Q15
        BBE_SVNX16_IP(av[0], pV, sizeof(*pV));
        BBE_SVNX16_IP(av[1], pV, sizeof(*pV));
        BBE_SVNX16_IP(afi[0], pV, sizeof(*pV));
    }

    pRRr = (const xb_vecNx16 *)R;
    pRRw = (xb_vecNx16 *)R;
    pV = (xb_vecNx16*)pScr;

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        xb_vecNx16 _0x4000 = BBE_MOVPINT16(16);
        vsaN _14 = BBE_MOVVSA32(14);
        xb_vecNx16 av[CURRENT_M];
        xb_vecNx16 afi[CURRENT_M];
        xb_vecNx16 tmp0, tmp1;

        BBE_LVNX16_XP(ar[0], pRRr, 4 * L);
        BBE_LVNX16_XP(ar[2], pRRr, 4 * L);
        BBE_LVNX16_XP(ar[1], pRRr, 4 * L);
        BBE_LVNX16_XP(ar[3], pRRr, 2 * BBE_SIMD_WIDTH - 3 * 4 * L);

        BBE_LVNX16_IP(av[0], pV, sizeof(*pV));
        BBE_LVNX16_IP(av[1], pV, sizeof(*pV));
        BBE_LVNX16_IP(afi[0], pV, sizeof(*pV));

        acc0 = BBE_MULRNX16J(ar[0], av[0], _14);
        BBE_MULANX16J(acc0, ar[1], av[1]);
        vR0 = BBE_PACKVNX40(acc0, _14);
        acc0 = BBE_MULRNX16(ar[0], _0x4000, _14);
        BBE_MULSNX16C(acc0, vR0, av[0]);
        ar[0] = BBE_PACKVNX40(acc0, _14);
        ar[1] = 0;
        acc0 = BBE_MULRNX16J(ar[2], av[0], _14);
        BBE_MULANX16J(acc0, ar[3], av[1]);
        vR1 = BBE_PACKVNX40(acc0, _14);
        acc0 = BBE_MULRNX16(ar[2], _0x4000, _14);
        BBE_MULSNX16C(acc0, vR1, av[0]);
        ar[2] = BBE_PACKVNX40(acc0, _14);
        acc0 = BBE_MULRNX16(ar[3], _0x4000, _14);
        BBE_MULSNX16C(acc0, vR1, av[1]);
        ar[3] = BBE_PACKVNX40(acc0, _14);

        {
            vsaN c_vec;
            acc0 = BBE_MAGINX16C(ar[3], ar[3]);
            acc0 = BBE_ADDNX40(acc0, acc0);
            c_vec = BBE_NSAENX40(acc0);
            acc0 = BBE_SLLNX40(acc0, c_vec);
            BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
            BBE_MULUUSNX16(acc0, tmp1, tmp0);
            acc0 = BBE_SRAINX40(acc0, 24);
            tmp0 = BBE_PACKLNX40(acc0);
            tmp0 = BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
            c_vec = BBE_SUBSR1SAVSN(18, c_vec);
            acc0 = BBE_MULUSNX16(tmp0, ar[3]);
            afi[1] = BBE_PACKVNX40(acc0, c_vec);
        }

        ar[0] = BBE_MULNX16JPACKQ(ar[0], afi[0]);
        ar[2] = BBE_MULNX16JPACKQ(ar[2], afi[0]);
        ar[3] = BBE_MULNX16JPACKQ(ar[3], afi[1]);
        BBE_SVNX16_XP(ar[0], pRRw, 4 * L);
        BBE_SVNX16_XP(ar[2], pRRw, 4 * L);
        BBE_SVNX16_XP(ar[1], pRRw, 4 * L);
        BBE_SVNX16_XP(ar[3], pRRw, 2 * BBE_SIMD_WIDTH - 3 * 4 * L);

        {
            vsaN _14 = BBE_MOVVSA32(14);
            xb_vecNx16 q00, q01, q11, q10;

            acc0 = BBE_MULRNX16J(av[0], av[0], _14); //q00
            q00 = BBE_PACKVNX40(acc0, _14);
            q00 = BBE_SUBSNX16(one, q00);

            acc0 = BBE_MULRNX16J(av[0], av[1], _14); //q01
            q01 = BBE_PACKVNX40(acc0, _14);
            q01 = BBE_NEGSNX16(q01);
            q10 = BBE_CONJNX16C(q01);

            acc0 = BBE_MULRNX16J(av[1], av[1], _14); //q11
            q11 = BBE_PACKVNX40(acc0, _14);
            q11 = BBE_SUBSNX16(one, q11);

            acc0 = BBE_MULRNX16C(q00, afi[0], _14); //q00
            q00 = BBE_PACKVNX40(acc0, _14); // result Q.15

            acc0 = BBE_MULRNX16C(q10, afi[0], _14); //q01
            q10 = BBE_PACKVNX40(acc0, _14);

            acc0 = BBE_MULRNX16C(q01, afi[1], _14); //q10
            q01 = BBE_PACKVNX40(acc0, _14);

            acc0 = BBE_MULRNX16C(q11, afi[1], _14); //q11
            q11 = BBE_PACKVNX40(acc0, _14);

            BBE_SVNX16_IP(q00, pQ00, sizeof(*pQ00));
            BBE_SVNX16_IP(q01, pQ01, sizeof(*pQ01));
            BBE_SVNX16_IP(q10, pQ10, sizeof(*pQ10));
            BBE_SVNX16_IP(q11, pQ11, sizeof(*pQ11));
        }
    }
}

size_t cqr2x2s_getScratchSize (int M, int N,int L)
{
    (void)M; (void)N;
    return 4 * 2 * (L)*sizeof(int16_t);
} /* cqr2x2s_getScratchSize() */
#else
DISCARD_FUN(void, cqr2x2s, (void* pScr,
    complex_fract16 * restrict Q, 
    complex_fract16 * restrict R, 
    int L))
size_t cqr2x2s_getScratchSize (int M, int N,int L){(void)M;(void)N;(void)L;return 0;}
#endif

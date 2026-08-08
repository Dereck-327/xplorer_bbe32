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
    cqr_build_rMxNs
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

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, cqr_build_r4x4s, ( void* pScr,complex_fract16 * restrict V, complex_fract16 * restrict R, int L))
#else

#define CURRENT_M 4
//#include "qr_common.h"

#define BBE_MULANX16CPACKQ(_result, _vR0, _v0)          \
{                                                       \
    xb_vecNx16 tmp = BBE_MULNX16CPACKQ(_vR0, _v0);      \
    _result = BBE_ADDNX16(tmp, _result);                \
}

#define DEF_CHOUSEHOLDER_FN_BODY                                            \
        xb_vecNx40 acc0;                                                    \
        xb_vecNx16 tmp0,tmp1;                                               \
        vsaN v_exp,c_vec;                                                   \
        int i = 0 ;                                                         \
        xb_vecNx16 k_norm, mant;                                            \
        xb_vecNx16 norm_x, Fi;                                              \
        xb_vecNx16  x1, x2;                                                 \
                                                                            \
        BBE_LVNX16_XP(x1, px, stride_bytes);                                \
        acc0 = BBE_MAGINX16C(x1, x1);                                       \
        acc0=BBE_ADDNX40(acc0,acc0);                                        \
        c_vec=BBE_NSAENX40(acc0);                                           \
        acc0=BBE_SLLNX40(acc0,c_vec);                                       \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                           \
        BBE_MULUUSNX16( acc0, tmp1,  tmp0);                                 \
        acc0=BBE_SRAINX40(acc0,24);                                         \
        tmp0=BBE_PACKLNX40(acc0);                                           \
        tmp0= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);              \
        c_vec = BBE_SUBSR1SAVSN(18+1,c_vec);                                \
        acc0 = BBE_MULUSRNX16(tmp0,x1, c_vec);                              \
        norm_x = BBE_PACKVNX40(acc0, c_vec);                                \
        Fi = norm_x;                                                        \
        BBE_SVNX16_XP(norm_x, pfi, V_Lstep*2*BBE_SIMD_WIDTH);               \
                                                                            \
      /* Build Householder's vector   v = x/sqrt(x'*x) + e1*x(1)/abs(x(1)*/ \
        acc0 = BBE_MAGINX16C( x1, x1);                                      \
        if ((len_x&1)==0)                                                   \
        {                                                                   \
            BBE_LVNX16_XP(x2, px, stride_bytes);                            \
            BBE_MAGIANX16C(acc0, x2, x2);                                   \
        }                                                                   \
        for(i=2-(len_x&1); i<len_x; i+=2)                                   \
        {                                                                   \
          BBE_LVNX16_XP(x1, px, stride_bytes);                              \
          BBE_LVNX16_XP(x2, px, stride_bytes);                              \
          BBE_MAGIANX16C(acc0, x1, x1);                                     \
          BBE_MAGIANX16C(acc0, x2, x2);                                     \
        }                                                                   \
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px);                     \
        BBE_LVNX16_XP(x1, px, stride_bytes);                                \
        acc0=BBE_ADDNX40(acc0,acc0);                                        \
        c_vec=BBE_NSAENX40(acc0);                                           \
        acc0=BBE_SLLNX40(acc0,c_vec);                                       \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                           \
        BBE_MULUUSNX16( acc0, tmp1,  tmp0);                                 \
        acc0=BBE_SRAINX40(acc0,24);                                         \
        tmp0=BBE_PACKLNX40(acc0);                                           \
        mant= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);              \
        v_exp = BBE_SUBSR1SAVSN(18+1,c_vec);                                \
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);                              \
        norm_x = BBE_PACKVNX40(acc0, v_exp);                                \
        norm_x = BBE_ADDNX16(norm_x, Fi);                                   \
                                                                            \
        acc0 = BBE_MULNX16J( norm_x, Fi);                                   \
        acc0=BBE_ADDNX40(acc0,acc0);                                        \
        c_vec=BBE_NSAENX40(acc0);                                           \
        acc0=BBE_SLLNX40(acc0,c_vec);                                       \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                           \
        BBE_MULUUSNX16(acc0, tmp1,  tmp0);                                  \
        c_vec = BBE_SUBSR1SAVSN(28,c_vec);                                  \
        k_norm=BBE_PACKVNX40(acc0,c_vec);                                   \
        k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);         \
                                                                            \
        acc0 = BBE_MULNX16( norm_x, k_norm);                                \
        norm_x = BBE_PACKQNX40(acc0);                                       \
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));

inline_  void chouseholder_M4xL_m0(
    const int16_t * A,        /* input matrix              */
    int16_t         * v,      /* output streaming order */
    int16_t         * Fi_out, /* output phase rotator      */
    const int L
    )
{
    const int m = 0;
    int l;
    const int stride_bytes = L * 4 * CURRENT_M;
    const xb_vecNx16  * restrict px = (xb_vecNx16*)(A + L * 2 * (CURRENT_M + 1)*m);
    xb_vecNx16  * restrict pv = (xb_vecNx16*)(v);
    xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out;
    const int V_Lstep = (SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
    const int len_x = CURRENT_M - m;

    NASSERT_ALIGN(px, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(v, (BBE_SIMD_WIDTH * 2));
    NASSERT(len_x > 1);
    NASSERT(CURRENT_M > 1);
    NASSERT(m < CURRENT_M - 1);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        DEF_CHOUSEHOLDER_FN_BODY;
        {
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
        }
        pv += (V_Lstep - (CURRENT_M - m));
        px -= (stride_bytes*(CURRENT_M - m)) / (2 * BBE_SIMD_WIDTH) - 1;
    } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_  void chouseholder_M4xL_m1(
    const int16_t * A,        /* input matrix              */
    int16_t         * v,      /* output streaming order */
    int16_t         * Fi_out, /* output phase rotator      */
    const int L
    )
{
    const int m = 1;
    int l;
    const int stride_bytes = L * 4 * CURRENT_M;
    const xb_vecNx16  * restrict px = (xb_vecNx16*)(A + L * 2 * (CURRENT_M + 1)*m);
    xb_vecNx16  * restrict pv = (xb_vecNx16*)(v);
    xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out;
    const int V_Lstep = (SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
    const int len_x = CURRENT_M - m;

    NASSERT_ALIGN(px, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(v, (BBE_SIMD_WIDTH * 2));
    NASSERT(len_x > 1);
    NASSERT(CURRENT_M > 1);
    NASSERT(m < CURRENT_M - 1);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        DEF_CHOUSEHOLDER_FN_BODY;
        {
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
        }
        pv += (V_Lstep - (CURRENT_M - m));
        px -= (stride_bytes*(CURRENT_M - m)) / (2 * BBE_SIMD_WIDTH) - 1;
    } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_  void chouseholder_M4xL_m2(
    const int16_t * A,        /* input matrix              */
    int16_t         * v,      /* output streaming order */
    int16_t         * Fi_out, /* output phase rotator      */
    const int L
    )
{
    const int m = 2;
    int l;
    const int stride_bytes = L * 4 * CURRENT_M;
    const xb_vecNx16  * restrict px = (xb_vecNx16*)(A + L * 2 * (CURRENT_M + 1)*m);
    xb_vecNx16  * restrict pv = (xb_vecNx16*)(v);
    xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out;
    const int V_Lstep = (SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
    const int len_x = CURRENT_M - m;

    NASSERT_ALIGN(px, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(v, (BBE_SIMD_WIDTH * 2));
    NASSERT(len_x > 1);
    NASSERT(CURRENT_M > 1);
    NASSERT(m < CURRENT_M - 1);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        DEF_CHOUSEHOLDER_FN_BODY;
        {
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2 * BBE_SIMD_WIDTH));
        }
        pv += (V_Lstep - (CURRENT_M - m));
        px -= (stride_bytes*(CURRENT_M - m)) / (2 * BBE_SIMD_WIDTH) - 1;
    } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

/*-----------------------------------------------------------------------
[c]qr_build_rMxNs

QR decomposition of MxN complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Fixed point representation of output matrices R is the same as for
input matrices, but Householder vectors V are always Q14.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]                  matrices A (L matrices of size MxN)
Output:
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
R[M*N][L]                  upper triangular matrices (L matrices
                           of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
void cqr_build_r4x4s  ( void* pScr, complex_fract16 * restrict _V, complex_fract16 * restrict _R, int L)
{
    int16_t * restrict V=(int16_t *)_V;
    int16_t * restrict R=(int16_t *)_R;
    int l;
    int m = 0;
    int V_Lstep = (SIZE_OF_V(CURRENT_M, CURRENT_M) + SIZE_OF_FI(CURRENT_M));
    //static const int v_offset[3] = { 0, 4, 4+3}; 

    int Fi_Offset = SIZE_OF_V(CURRENT_M, CURRENT_M);
    xb_vecNx16  vR0;
    xb_vecNx40 acc0;
    xb_vecNx16 * restrict pV;
    xb_vecNx16 * restrict pcol0;
    xb_vecNx16 * restrict pcol1;
    xb_vecNx16 * restrict pcol2;
    xb_vecNx16 * restrict pcol3;
    xb_vecNx16 * restrict pfi;
    xb_vecNx16 * restrict pcol01;
    xb_vecNx16 * restrict pcol11;

#ifdef COMPILER_XTENSA
#pragma ymemory( pcol0 )
#pragma ymemory( pcol1 )
#pragma ymemory( pcol2 )
#pragma ymemory( pcol3 )
#pragma ymemory( pcol01)
#pragma ymemory( pcol11)
#endif


    const int col_step = CURRENT_M*L * 4;
    int v_last_step = (V_Lstep*BBE_SIMD_WIDTH - BBE_SIMD_WIDTH * 3) * 2;
    int col_back_step = -3 * (CURRENT_M*L * 4);
    int col_next_step = col_back_step + 2 * BBE_SIMD_WIDTH;
    xb_vecNx16 av[CURRENT_M];
    xb_vecNx16 ar[CURRENT_M];

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    m = 0;
    {
        pcol0 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * m);
        pcol1 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 1));
        pcol2 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 2));
        pcol3 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 3));
        chouseholder_M4xL_m0(R,
            V + BBE_SIMD_WIDTH * 0,
            V + BBE_SIMD_WIDTH*(Fi_Offset + m),
            L);

        pV = (xb_vecNx16*)(V + BBE_SIMD_WIDTH * 0);
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            vsaN _13 = BBE_MOVVSA32(13);
            xb_vecNx16 tmp;
            // ********* Process column 0 ***********
            BBE_LVNX16_IP(av[0], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(av[1], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(av[2], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(av[3], pV, v_last_step);

            BBE_LVNX16_XP(ar[0], pcol0, col_step);
            BBE_LVNX16_XP(ar[1], pcol0, col_step);
            BBE_LVNX16_XP(ar[2], pcol0, col_step);
            BBE_LVNX16_XP(ar[3], pcol0, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);
            BBE_MULANX16J(acc0, ar[3], av[3]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            tmp = BBE_MULNX16CPACKQ(vR0, av[0]);
            ar[0] = BBE_SUBNX16(ar[0], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[1]);
            ar[1] = BBE_SUBNX16(ar[1], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[2]);
            ar[2] = BBE_SUBNX16(ar[2], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[3]);
            ar[3] = BBE_SUBNX16(ar[3], tmp);

            BBE_SVNX16_XP(ar[0], pcol0, col_step);
            BBE_SVNX16_XP(ar[1], pcol0, col_step);
            BBE_SVNX16_XP(ar[2], pcol0, col_step);
            BBE_SVNX16_XP(ar[3], pcol0, col_next_step);

            // ********* Process column 1 ***********

            BBE_LVNX16_XP(ar[0], pcol1, col_step);
            BBE_LVNX16_XP(ar[1], pcol1, col_step);
            BBE_LVNX16_XP(ar[2], pcol1, col_step);
            BBE_LVNX16_XP(ar[3], pcol1, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);
            BBE_MULANX16J(acc0, ar[3], av[3]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            tmp = BBE_MULNX16CPACKQ(vR0, av[0]);
            ar[0] = BBE_SUBNX16(ar[0], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[1]);
            ar[1] = BBE_SUBNX16(ar[1], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[2]);
            ar[2] = BBE_SUBNX16(ar[2], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[3]);
            ar[3] = BBE_SUBNX16(ar[3], tmp);

            BBE_SVNX16_XP(ar[0], pcol1, col_step);
            BBE_SVNX16_XP(ar[1], pcol1, col_step);
            BBE_SVNX16_XP(ar[2], pcol1, col_step);
            BBE_SVNX16_XP(ar[3], pcol1, col_next_step);
        }

        pV = (xb_vecNx16*)(V + BBE_SIMD_WIDTH * 0);
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            vsaN _13 = BBE_MOVVSA32(13);
            xb_vecNx16 tmp;
            // ********* Process column 2 ***********
            BBE_LVNX16_IP(av[0], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(av[1], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(av[2], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(av[3], pV, v_last_step);

            BBE_LVNX16_XP(ar[0], pcol2, col_step);
            BBE_LVNX16_XP(ar[1], pcol2, col_step);
            BBE_LVNX16_XP(ar[2], pcol2, col_step);
            BBE_LVNX16_XP(ar[3], pcol2, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);
            BBE_MULANX16J(acc0, ar[3], av[3]);

            vR0 = BBE_PACKVNX40(acc0, _13);

            tmp = BBE_MULNX16CPACKQ(vR0, av[0]);
            ar[0] = BBE_SUBNX16(ar[0], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[1]);
            ar[1] = BBE_SUBNX16(ar[1], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[2]);
            ar[2] = BBE_SUBNX16(ar[2], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[3]);
            ar[3] = BBE_SUBNX16(ar[3], tmp);

            BBE_SVNX16_XP(ar[0], pcol2, col_step);
            BBE_SVNX16_XP(ar[1], pcol2, col_step);
            BBE_SVNX16_XP(ar[2], pcol2, col_step);
            BBE_SVNX16_XP(ar[3], pcol2, col_next_step);
            // ********* Process column 4 ***********

            BBE_LVNX16_XP(ar[0], pcol3, col_step);
            BBE_LVNX16_XP(ar[1], pcol3, col_step);
            BBE_LVNX16_XP(ar[2], pcol3, col_step);
            BBE_LVNX16_XP(ar[3], pcol3, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);
            BBE_MULANX16J(acc0, ar[3], av[3]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            tmp = BBE_MULNX16CPACKQ(vR0, av[0]);
            ar[0] = BBE_SUBNX16(ar[0], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[1]);
            ar[1] = BBE_SUBNX16(ar[1], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[2]);
            ar[2] = BBE_SUBNX16(ar[2], tmp);
            tmp = BBE_MULNX16CPACKQ(vR0, av[3]);
            ar[3] = BBE_SUBNX16(ar[3], tmp);

            BBE_SVNX16_XP(ar[0], pcol3, col_step);
            BBE_SVNX16_XP(ar[1], pcol3, col_step);
            BBE_SVNX16_XP(ar[2], pcol3, col_step);
            BBE_SVNX16_XP(ar[3], pcol3, col_next_step);
        }
    }
    //print_mtx(startR, 4, 4, 0, L, 0, 3, 1, "R(m=0)"); 
    m = 1;
    col_back_step += col_step;
    col_next_step += col_step;
    v_last_step += 2 * BBE_SIMD_WIDTH;
    {


        pcol0 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * m);
        pcol1 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 1));
        pcol2 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 2));
        pV = (xb_vecNx16*)(V + BBE_SIMD_WIDTH * 4);

        chouseholder_M4xL_m1(R,
            V + BBE_SIMD_WIDTH * 4,
            V + BBE_SIMD_WIDTH*(Fi_Offset + m),
            L);
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            vsaN _13 = BBE_MOVVSA32(13);
            // ********* Process column 0 ***********
            BBE_LVNX16_IP(av[0], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(av[1], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(av[2], pV, v_last_step);

            BBE_LVNX16_XP(ar[0], pcol0, col_step);
            BBE_LVNX16_XP(ar[1], pcol0, col_step);
            BBE_LVNX16_XP(ar[2], pcol0, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            vR0 = BBE_NEGNX16(vR0);

            BBE_MULANX16CPACKQ(ar[0], vR0, av[0]);
            BBE_MULANX16CPACKQ(ar[1], vR0, av[1]);
            BBE_MULANX16CPACKQ(ar[2], vR0, av[2]);

            BBE_SVNX16_XP(ar[0], pcol0, col_step);
            BBE_SVNX16_XP(ar[1], pcol0, col_step);
            BBE_SVNX16_XP(ar[2], pcol0, col_next_step);

            // ********* Process column 1 ***********

            BBE_LVNX16_XP(ar[0], pcol1, col_step);
            BBE_LVNX16_XP(ar[1], pcol1, col_step);
            BBE_LVNX16_XP(ar[2], pcol1, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            vR0 = BBE_NEGNX16(vR0);

            BBE_MULANX16CPACKQ(ar[0], vR0, av[0]);
            BBE_MULANX16CPACKQ(ar[1], vR0, av[1]);
            BBE_MULANX16CPACKQ(ar[2], vR0, av[2]);

            BBE_SVNX16_XP(ar[0], pcol1, col_step);
            BBE_SVNX16_XP(ar[1], pcol1, col_step);
            BBE_SVNX16_XP(ar[2], pcol1, col_next_step);

            // ********* Process column 2 ***********

            BBE_LVNX16_XP(ar[0], pcol2, col_step);
            BBE_LVNX16_XP(ar[1], pcol2, col_step);
            BBE_LVNX16_XP(ar[2], pcol2, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);
            BBE_MULANX16J(acc0, ar[2], av[2]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            vR0 = BBE_NEGNX16(vR0);

            BBE_MULANX16CPACKQ(ar[0], vR0, av[0]);
            BBE_MULANX16CPACKQ(ar[1], vR0, av[1]);
            BBE_MULANX16CPACKQ(ar[2], vR0, av[2]);

            BBE_SVNX16_XP(ar[0], pcol2, col_step);
            BBE_SVNX16_XP(ar[1], pcol2, col_step);
            BBE_SVNX16_XP(ar[2], pcol2, col_next_step);
        }
    }
    //print_mtx(startR, 4, 4, 0, L, 0, 3, 1, "R(m=1)"); 
    m = 2;
    col_back_step += col_step;
    col_next_step += col_step;
    v_last_step += 2 * BBE_SIMD_WIDTH;

    {
        pcol0 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * m);
        pcol1 = (xb_vecNx16 *)(R + CURRENT_M*L * 2 * m + L * 2 * (m + 1));
        pcol01 = pcol0 + col_step / sizeof(*pcol01);
        pcol11 = pcol1 + col_step / sizeof(*pcol01);

        pV = (xb_vecNx16*)(V + BBE_SIMD_WIDTH * 7);
        chouseholder_M4xL_m2(R,
            V + BBE_SIMD_WIDTH * 7,
            V + BBE_SIMD_WIDTH*(Fi_Offset + m),
            L);
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            vsaN _13 = BBE_MOVVSA32(13);
            // ********* Process column 0 ***********
            BBE_LVNX16_IP(av[0], pV, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(av[1], pV, v_last_step);

            ar[0] = BBE_LVNX16_I(pcol0, 0);
            ar[1] = BBE_LVNX16_I(pcol01, 0);
            ar[2] = BBE_LVNX16_I(pcol1, 0);
            ar[3] = BBE_LVNX16_I(pcol11, 0);

            acc0 = BBE_MULRNX16J(ar[0], av[0], _13);
            BBE_MULANX16J(acc0, ar[1], av[1]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            vR0 = BBE_NEGNX16(vR0);

            BBE_MULANX16CPACKQ(ar[0], vR0, av[0]);
            BBE_MULANX16CPACKQ(ar[1], vR0, av[1]);

            BBE_SVNX16_IP(ar[0], pcol0, (2 * BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(ar[1], pcol01, (2 * BBE_SIMD_WIDTH));

            // ********* Process column 1 ***********
            acc0 = BBE_MULRNX16J(ar[2], av[0], _13);
            BBE_MULANX16J(acc0, ar[3], av[1]);

            vR0 = BBE_PACKVNX40(acc0, _13);
            vR0 = BBE_NEGNX16(vR0);

            BBE_MULANX16CPACKQ(ar[2], vR0, av[0]);
            BBE_MULANX16CPACKQ(ar[3], vR0, av[1]);

            BBE_SVNX16_IP(ar[2], pcol1, (2 * BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(ar[3], pcol11, (2 * BBE_SIMD_WIDTH));
        }
    }
    //print_mtx(startR, 4, 4, 0, L, 0, 3, 1, "R(m=2)"); 
    m = 3;

    {
        xb_vecNx16  x0, _fi;
        xb_vecNx16 *  pa;
        xb_vecNx16 *  pfi = (xb_vecNx16 *)(V + BBE_SIMD_WIDTH*(Fi_Offset + m));
        pa = (xb_vecNx16 *)(R + m*L * 2 + m*CURRENT_M*L * 2);
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            xb_vecNx16 tmp0, tmp1;
            vsaN c_vec, v_exp;
            xb_vecNx40 acc0;
            xb_vecNx16  mant;
            BBE_LVNX16_IP(x0, pa, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MAGINX16C(x0, x0);
            acc0 = BBE_ADDNX40(acc0, acc0);
            c_vec = BBE_NSAENX40(acc0);
            acc0 = BBE_SLLNX40(acc0, c_vec);
            BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
            BBE_MULUUSNX16(acc0, tmp1, tmp0);
            acc0 = BBE_SRAINX40(acc0, 24);
            tmp0 = BBE_PACKLNX40(acc0);
            mant = BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
            v_exp = BBE_SUBSR1SAVSN(18 + 1, c_vec);
            acc0 = BBE_MULUSRNX16(mant, x0, v_exp);
            _fi = BBE_PACKVNX40(acc0, v_exp);
            BBE_SVNX16_XP(_fi, pfi, V_Lstep * 2 * BBE_SIMD_WIDTH);
        }
    }
    //print_mtx(startR, 4, 4, 0, L, 0, 3, 1, "R(m=3)"); 
    {
        pfi = (xb_vecNx16 *)(V + BBE_SIMD_WIDTH*(Fi_Offset + 0));
        pcol0 = (xb_vecNx16 *)(R + L * 2 * 0);
        pcol1 = (xb_vecNx16 *)(R + L * 2 * 1);
        pcol2 = (xb_vecNx16 *)(R + L * 2 * 2);
        pcol3 = (xb_vecNx16 *)(R + L * 2 * 3);
        col_back_step = -3 * (CURRENT_M*L * 4);
        col_next_step = col_back_step + 2 * BBE_SIMD_WIDTH;


        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            vsaN _14 = BBE_MOVVSA32(14);
            xb_vecNx16 fi0, fi1, fi2, fi3;
            BBE_LVNX16_IP(fi0, pfi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(fi1, pfi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(fi2, pfi, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(fi3, pfi, V_Lstep * 2 * BBE_SIMD_WIDTH - 3 * 2 * BBE_SIMD_WIDTH);

            BBE_LVNX16_XP(ar[0], pcol0, col_step);
            BBE_LVNX16_XP(ar[1], pcol0, col_step);
            BBE_LVNX16_XP(ar[2], pcol0, col_step);
            BBE_LVNX16_XP(ar[3], pcol0, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], fi0, _14); ar[0] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[1], fi1, _14); ar[1] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[2], fi2, _14); ar[2] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[3], fi3, _14); ar[3] = BBE_PACKVNX40(acc0, _14);

            BBE_SVNX16_XP(ar[0], pcol0, col_step);
            BBE_SVNX16_XP(ar[1], pcol0, col_step);
            BBE_SVNX16_XP(ar[2], pcol0, col_step);
            BBE_SVNX16_XP(ar[3], pcol0, col_next_step);

            BBE_LVNX16_XP(ar[0], pcol1, col_step);
            BBE_LVNX16_XP(ar[1], pcol1, col_step);
            BBE_LVNX16_XP(ar[2], pcol1, col_step);
            BBE_LVNX16_XP(ar[3], pcol1, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], fi0, _14); ar[0] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[1], fi1, _14); ar[1] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[2], fi2, _14); ar[2] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[3], fi3, _14); ar[3] = BBE_PACKVNX40(acc0, _14);

            BBE_SVNX16_XP(ar[0], pcol1, col_step);
            BBE_SVNX16_XP(ar[1], pcol1, col_step);
            BBE_SVNX16_XP(ar[2], pcol1, col_step);
            BBE_SVNX16_XP(ar[3], pcol1, col_next_step);

            BBE_LVNX16_XP(ar[0], pcol2, col_step);
            BBE_LVNX16_XP(ar[1], pcol2, col_step);
            BBE_LVNX16_XP(ar[2], pcol2, col_step);
            BBE_LVNX16_XP(ar[3], pcol2, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], fi0, _14); ar[0] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[1], fi1, _14); ar[1] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[2], fi2, _14); ar[2] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[3], fi3, _14); ar[3] = BBE_PACKVNX40(acc0, _14);

            BBE_SVNX16_XP(ar[0], pcol2, col_step);
            BBE_SVNX16_XP(ar[1], pcol2, col_step);
            BBE_SVNX16_XP(ar[2], pcol2, col_step);
            BBE_SVNX16_XP(ar[3], pcol2, col_next_step);

            BBE_LVNX16_XP(ar[0], pcol3, col_step);
            BBE_LVNX16_XP(ar[1], pcol3, col_step);
            BBE_LVNX16_XP(ar[2], pcol3, col_step);
            BBE_LVNX16_XP(ar[3], pcol3, col_back_step);

            acc0 = BBE_MULRNX16J(ar[0], fi0, _14); ar[0] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[1], fi1, _14); ar[1] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[2], fi2, _14); ar[2] = BBE_PACKVNX40(acc0, _14);
            acc0 = BBE_MULRNX16J(ar[3], fi3, _14); ar[3] = BBE_PACKVNX40(acc0, _14);

            BBE_SVNX16_XP(ar[0], pcol3, col_step);
            BBE_SVNX16_XP(ar[1], pcol3, col_step);
            BBE_SVNX16_XP(ar[2], pcol3, col_step);
            BBE_SVNX16_XP(ar[3], pcol3, col_next_step);
        }
    }
    //print_mtx(startR, 4, 4, 0, L, 0, 3, 1, "R(final)"); 
} /* cqr_build_r4x4s() */
#endif

size_t cqr_build_r4x4s_getScratchSize (int M, int N,int L)
{
    (void)M; (void)N; (void)L;
    return 0;
} /* cqr_build_r4x4s_getScratchSize() */

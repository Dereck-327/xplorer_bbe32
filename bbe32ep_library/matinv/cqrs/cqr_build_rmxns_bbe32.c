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

#if (HAVE_VSAMATH && HAVE_NSAENX40 && 1)
/*
    find Housholder vector by column and transform that column
    */
static void cqrHouseholder(
    void* pScr,          /* scratch (4*L bytes)        */
    const int16_t *   A, /* input matrix               */
    int16_t *    v,     /* output (housholder vector) */
    int16_t*    Fi_out,  /* output phase rotator       */
    const int M,         /* number of rows             */
    const int N,         /* number of columns          */
    const int len_x,     /* length of given column     */
    const int L          /* number of matrices         */
    )
{
    int i, l;
    xb_vecNx16* restrict _pwrs;   // write pointer to the scratch data
    const xb_vecNx16* restrict _prds;   // read pointer to the scratch data
    const xb_vecNx16* restrict px;
    xb_vecNx16* restrict pFi;
    xb_vecNx16 *restrict pv;
    const xb_vecNx16 *px_v16;
    const int V_Lstep = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    xb_vecNx16 x, t, mant, k_norm, norm_x, fi;
    vsaN exp;
    xb_vecNx40 a;

    NASSERT_ALIGN(v, BBE_SIMD_WIDTH * 2);
    NASSERT(len_x > 1);
    NASSERT(M > 1);
    NASSERT(len_x > 1);

    //-----------------------------
    // compute diagonal element
    //-----------------------------
    pFi = (xb_vecNx16*)Fi_out;
    px = (xb_vecNx16*)A;
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        BBE_LVNX16_IP(norm_x, px, 2 * BBE_SIMD_WIDTH);
        a = BBE_MULNX16J(norm_x, norm_x);
        a = BBE_ADDNX40(a, a);
        exp = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, exp);
        BBE_RSQRTLUNX40_0(a, x, t, a);
        BBE_MULUUSNX16(a, t, x);
        a = BBE_SRAINX40(a, 23);
        mant = BBE_PACKLNX40(a);
        exp = BBE_SUBSR1SAVSN(16 + 4, exp);
        x = BBE_MOVVVS(exp);
        mant = BBE_SHFLNX16I(mant, BBE_SHFLI_DUPLICATE_1_EVEN);
        x = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_EVEN);
        exp = BBE_MOVVSV(x, 0);
        a = BBE_MULUSRNX16(mant, norm_x, exp);
        norm_x = BBE_PACKVNX40(a, exp);
        BBE_SVNX16_XP(norm_x, pFi, V_Lstep * 2 * BBE_SIMD_WIDTH);
    }
    //-----------------------------
    // compute rotation vector
    //-----------------------------
    _pwrs = (xb_vecNx16*)pScr;
    px = (xb_vecNx16*)A;
    if (len_x == 2)
    {
        px_v16 = (const xb_vecNx16*)px;
        for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
        {
            BBE_LVNX16_XP(x, px_v16, 4 * N*L);
            a = BBE_MULNX16J(x, x);
            BBE_LVNX16_XP(x, px_v16, 2 * BBE_SIMD_WIDTH - 4 * N*L);
            BBE_MULANX16J(a, x, x);
            a = BBE_ADDNX40(a, a);
            exp = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, exp);
            BBE_RSQRTLUNX40_0(a, x, t, a);
            BBE_MULUUSNX16(a, t, x);
            a = BBE_SRAINX40(a, 23);
            mant = BBE_PACKLNX40(a);
            exp = BBE_SUBSR1SAVSN(16 + 4, exp);
            x = BBE_MOVVVS(exp);
            x = BBE_SELNX16I(x, mant, BBE_SELI_INTERLEAVE_1_EVEN);
            BBE_SVNX16_IP(x, _pwrs, 2 * BBE_SIMD_WIDTH);
        }
    }
    else
    {
        for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
        {
            px_v16 = (const xb_vecNx16*)px;
            a = 0;
            __Pragma("loop_count min=3")
                for (i = 0; i < len_x; i++)
                {
                    BBE_LVNX16_XP(x, px_v16, 4 * N*L);
                    BBE_MULANX16J(a, x, x);
                }
            a = BBE_ADDNX40(a, a);
            exp = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, exp);
            BBE_RSQRTLUNX40_0(a, x, t, a);
            BBE_MULUUSNX16(a, t, x);
            a = BBE_SRAINX40(a, 23);
            mant = BBE_PACKLNX40(a);
            exp = BBE_SUBSR1SAVSN(16 + 4, exp);
            x = BBE_MOVVVS(exp);
            x = BBE_SELNX16I(x, mant, BBE_SELI_INTERLEAVE_1_EVEN);
            BBE_SVNX16_IP(x, _pwrs, 2 * BBE_SIMD_WIDTH);

            px += (BBE_SIMD_WIDTH * 2) / sizeof(*px);
        }
    }
    __Pragma("no_reorder")
        //-------------------------------------
        // affine projection of original column
        //-------------------------------------
        _prds = (const xb_vecNx16*)pScr;
    pFi = (xb_vecNx16*)Fi_out;
    pv = (xb_vecNx16*)(v);
    px = (xb_vecNx16*)A;
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        int i;
        BBE_LVNX16_XP(fi, pFi, V_Lstep * 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x, _prds, 2 * BBE_SIMD_WIDTH);
        mant = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_EVEN);
        x = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_ODD);
        exp = BBE_MOVVSV(x, 0);

        px_v16 = (const xb_vecNx16*)px;
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        norm_x = BBE_PACKVNX40(a, exp);
        norm_x = BBE_ADDNX16(norm_x, fi);
        {
            vsaN c_vec;
            a = BBE_MULNX16J(norm_x, fi);
            a = BBE_ADDNX40(a, a);
            c_vec = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, c_vec);
            c_vec = BBE_SUBSR1SAVSN(28, c_vec);
            BBE_RSQRTLUNX40_0(a, t, x, a);
            BBE_MULUUSNX16(a, x, t);
            k_norm = BBE_PACKVNX40(a, c_vec);
            k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);
        }
        norm_x = BBE_MULNX16PACKQ(k_norm, norm_x);
        BBE_SVNX16_IP(norm_x, pv, BBE_SIMD_WIDTH * 2);
        for (i = 1; i < len_x; i++)
        {
            BBE_LVNX16_XP(x, px_v16, 4 * N*L);
            a = BBE_MULUSRNX16(mant, x, exp);
            x = BBE_PACKVNX40(a, exp);
            x = BBE_MULNX16PACKQ(x, k_norm);
            BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        }
        pv += (V_Lstep - len_x)*(BBE_SIMD_WIDTH * 2) / sizeof(*pv);
        px += (BBE_SIMD_WIDTH * 2) / sizeof(*px);
    }
}

/*
    update columns of matrix R by Housholder vectors
    */
static void cqrUpdateR(int16_t* restrict R,
    const int16_t * restrict V,             /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int len_x, int L)
{
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r, v, vr, _0x4000 = 0x4000;
    int i, l, j, sizeV = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    NASSERT(len_x > 1);
    if (len_x == 2)
    {
        for (j = 0; j < N - m; j++)
        {
            pv = (xb_vecNx16*)(V);
            pr = (xb_vecNx16*)(R);
            pw = (xb_vecNx16*)(R);
            for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
            {
                xb_vecNx16 v0, v1, r0, r1;
                r1 = BBE_LVNX16_X(pr, 4 * N*L);
                BBE_LVNX16_IP(r0, pr, 2 * BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(v1, pv, 2 * BBE_SIMD_WIDTH*(sizeV - 1));
                acc = BBE_MULRNX16J(r0, v0, q14);
                BBE_MULANX16J(acc, r1, v1);
                vr = BBE_PACKVNX40(acc, q14);
                acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
                acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
                BBE_SVNX16_X(r1, pw, 4 * N*L);
                BBE_SVNX16_IP(r0, pw, 2 * BBE_SIMD_WIDTH);
#if 0
                //pr=(xb_vecNx16*)(R + 2*l);
                pw=(xb_vecNx16*)pr;
                BBE_LVNX16_XP(r, pr, 4*N*L); 
                BBE_LVNX16_IP(v, pv, 2*BBE_SIMD_WIDTH); 
                acc = BBE_MULRNX16J(r, v, q14); 
                for(i=1; i<len_x; i++)
                {
                    BBE_LVNX16_XP(r, pr, 4*N*L); 
                    BBE_LVNX16_IP(v,pv, 2*BBE_SIMD_WIDTH); 
                    BBE_MULANX16J(acc, r, v); 
                }
                vr=BBE_PACKVNX40(acc,q14);
                pv=(const xb_vecNx16*)XT_ADD(-2*BBE_SIMD_WIDTH*(len_x),(uintptr_t)pv);
                pr=(const xb_vecNx16*)pw;
                for(i=0; i<len_x; i++)
                {
                    BBE_LVNX16_XP(r,pr, 4*N*L); 
                    BBE_LVNX16_IP(v,pv, 2*BBE_SIMD_WIDTH); 
                    acc=BBE_MULRNX16(r,_0x4000,q14);
                    BBE_MULSNX16C(acc,v,vr);
                    r=BBE_PACKVNX40(acc,q14);
                    BBE_SVNX16_XP(r, pw, 4*N*L);
                }
                pv=(const xb_vecNx16*)XT_ADD(2*BBE_SIMD_WIDTH*(sizeV-len_x),(uintptr_t)pv);
                pr=(const xb_vecNx16*)XT_ADD(2*BBE_SIMD_WIDTH-len_x*4*N*L  ,(uintptr_t)pr);
#endif
            }
            R = (int16_t*)XT_ADDX4(L, (uintptr_t)R);
        }
    }
    else
    {
        for (j = 0; j < N - m; j++)
        {
            pv = (xb_vecNx16*)(V);
            pr = (xb_vecNx16*)(R);
            for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
            {
                //pr=(xb_vecNx16*)(R + 2*l);
                pw = (xb_vecNx16*)pr;
                BBE_LVNX16_XP(r, pr, 4 * N*L);
                BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                acc = BBE_MULRNX16J(r, v, q14);
                __Pragma("loop_count min=2")
                    for (i = 1; i < len_x; i++)
                    {
                        BBE_LVNX16_XP(r, pr, 4 * N*L);
                        BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                        BBE_MULANX16J(acc, r, v);
                    }
                vr = BBE_PACKVNX40(acc, q14);
                pv = (const xb_vecNx16*)XT_ADD(-2 * BBE_SIMD_WIDTH*(len_x), (uintptr_t)pv);
                pr = (const xb_vecNx16*)pw;
                //__Pragma("loop_count min=3") : can not be used here due to the compiler bug
                for (i = 0; i < len_x; i++)
                {
                    BBE_LVNX16_XP(r, pr, 4 * N*L);
                    BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                    acc = BBE_MULRNX16(r, _0x4000, q14);
                    BBE_MULSNX16C(acc, v, vr);
                    r = BBE_PACKVNX40(acc, q14);
                    BBE_SVNX16_XP(r, pw, 4 * N*L);
                }
                pv = (const xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH*(sizeV - len_x), (uintptr_t)pv);
                pr = (const xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH - len_x * 4 * N*L, (uintptr_t)pr);
            }
            R = (int16_t*)XT_ADDX4(L, (uintptr_t)R);
        }
    }
}


// final rotation of matrices R(M,N)xL
static void cqrRotateR(int16_t* R, const int16_t* pFi, int M, int N, int L)
{
    vsaN q14 = BBE_MOVVSA32(14);
    int k, i, j;
    xb_vecNx16  fi, r;
    xb_vecNx40  acc0;
    const xb_vecNx16* restrict _pFi;
    const xb_vecNx16* restrict _pr;
    xb_vecNx16* restrict _pw;

    _pr = (const xb_vecNx16*)(R);
    _pw = (xb_vecNx16*)(R);
    _pFi = (const xb_vecNx16*)pFi;
    NASSERT(N > 1);
    for (k = i = 0; i < M*(L >> (LOG2_BBE_SIMD_WIDTH - 1)); i++)
    {
        int strideFi, strideR;
        k = BBE_ADDMOD16U(k, (M << 16) | 1);
        strideFi = 2 * BBE_SIMD_WIDTH;
        strideR = 4 * L;
        XT_MOVEQZ(strideFi, 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(M, N) + 1), k);
        XT_MOVEQZ(strideR, -(M*N - 1) * 4 * L + 2 * BBE_SIMD_WIDTH, k);
        BBE_LVNX16_XP(fi, _pFi, strideFi);
        for (j = 0; j < N - 1; j++)
        {
            BBE_LVNX16_XP(r, _pr, 4 * L);
            acc0 = BBE_MULRNX16J(r, fi, q14);
            r = BBE_PACKVNX40(acc0, q14);
            BBE_SVNX16_XP(r, _pw, 4 * L);
        }
        BBE_LVNX16_XP(r, _pr, strideR);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, strideR);
    }
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
void cqr_build_rmxns  ( void* pScr, complex_fract16 * restrict _V, complex_fract16 * restrict _R, int M, int N,int L)
{
    int16_t * restrict V=(int16_t *)_V;
    int16_t * restrict R=(int16_t *)_R;
    int Ncolumns = M > N ? N : N - 1;  // update N columns for overdetermined linear equation system
    int m;
    int16_t * restrict pV;
    int16_t * restrict pFi;
    int16_t * restrict pR;

    NASSERT(M <= MAX_M);
    NASSERT(M >= N);
    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(V, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    pFi = V + BBE_SIMD_WIDTH* SIZE_OF_V(M, N);
    pR = R;
    pV = V;
    for (m = 0; m < Ncolumns; m++)
    {
        cqrHouseholder(pScr, pR, pV, pFi, M, N, (M - m), L);
        cqrUpdateR(pR, pV, m, M, N, (M - m), L);
        pFi += BBE_SIMD_WIDTH;
        pV += BBE_SIMD_WIDTH*(M - m);
        pR += (N + 1)*L * 2;
    }
    // for square matrices we have to calculate last rotation explicitely
    if (M == N)
    {
        cqrComputeLastRot(V + BBE_SIMD_WIDTH* SIZE_OF_V(M, N), R, M, L);
    }
    // final rotation
    pFi = V + BBE_SIMD_WIDTH* SIZE_OF_V(M, N);
    cqrRotateR(R, pFi, M, N, L);
}

size_t cqr_build_rmxns_getScratchSize (int M, int N,int L)
{
    (void)M; (void)N;
    return 4 * L;
} /* cqr_build_rmxns_getScratchSize() */
#else
DISCARD_FUN(void, cqr_build_rmxns, ( void* pScr,
    complex_fract16 * restrict V, 
    complex_fract16 * restrict R, 
    int M, 
    int N,
    int L))
size_t cqr_build_rmxns_getScratchSize (int M, int N,int L) { (void)M;(void)N;(void)L; return 0; }
#endif

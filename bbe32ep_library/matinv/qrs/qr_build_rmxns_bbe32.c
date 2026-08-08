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
#include "qr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, qr_build_rmxns,( void* pScr, int16_t * restrict V, int16_t * restrict R, int M,int N,int L))
#else

/*
  Build Householder's 16 vectors for column A(m:end,m)
  The A is input complex matrix MxM

  Algorithm:

  x = A(m:end,m);

  inv_norm2_x =  1/sqrt(x'*x);
  if abs(x(1))==0
  fi = 1;
  else
  fi = x(1)/abs(x(1);
  end

  alpha = - 1 * fi;
  v = x;
  v1 = x(1)*inv_norm2_x - alpha;
  k_norm = sqrt(2)/sqrt( 2 * real(v1 *conj(Fi(m)))) ;

  v(1) = v1*k_norm;
  v(2:end) = x(2:end)*inv_norm2_x*k_norm;

  */
static void qrHouseholder(
    const int16_t * restrict A, /* input matrix         */
    int16_t * restrict  v,     /* output streaming order       */
    const int M,
    const int N,
    const int len_x,
    const int L
    )
{
    int p;
    const xb_vecNx16 *restrict x = (const xb_vecNx16 *)A;
    xb_vecNx16 *restrict pv = (xb_vecNx16 *)v;
    xb_vecNx16 x0, t, norm_x, _0x4000 = 0x4000, d0, _0 = 0, k_norm;
    xb_vecNx40 a, b;
    vboolN sign;
    vsaN d1, exp;
    int i;
    int size_V = SIZE_OF_V(M, N);

    NASSERT(len_x > 1);
    NASSERT(M > 1);
    NASSERT(len_x > 1);
    if (len_x == 2)   // separate case
    {
        for (p = 0; p < L; p += BBE_SIMD_WIDTH)
        {
            BBE_LVNX16_XP(x0, x, 2 * L*N);
            a = BBE_MULNX16(x0, x0);
            BBE_LVNX16_XP(x0, x, -2 * L*N);
            BBE_MULANX16(a, x0, x0);

            a = BBE_ADDNX40(a, a);
            d1 = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, d1);
            BBE_RSQRTLUNX40_0(b, t, x0, a);
            BBE_RSQRTLUNX40_1(b, t, x0, a);
            BBE_MULUUSNX16(b, x0, t);
            a = BBE_SRAINX40(b, 23);
            d0 = BBE_PACKLNX40(a);
            d1 = BBE_SUBSR1SAVSN(16 + 4, d1);
            BBE_LVNX16_XP(x0, x, 2 * L*N);
            a = BBE_MULUSRNX16(d0, x0, d1);
            norm_x = BBE_PACKVNX40(a, d1);
            sign = BBE_LTNX16(x0, _0);
            BBE_ADDSNX16F(norm_x, norm_x, _0x4000, sign);
            BBE_SUBSNX16T(norm_x, norm_x, _0x4000, sign);
            x0 = BBE_ABSSNX16(norm_x);
            x0 = BBE_ADDNX16(x0, x0);
            a = BBE_UNPKUNX16(x0);

            exp = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, exp);
            exp = BBE_SUBSR1SAVSN(35, exp);
            BBE_RSQRTLUNX40_0(b, t, x0, a);
            BBE_RSQRTLUNX40_1(b, t, x0, a);
            BBE_MULUUSNX16(b, x0, t);
            k_norm = BBE_PACKVNX40(b, exp);
            x0 = BBE_MULNX16PACKQ(k_norm, norm_x);
            BBE_SVNX16_IP(x0, pv, 2 * BBE_SIMD_WIDTH);

            BBE_LVNX16_XP(x0, x, 2 * BBE_SIMD_WIDTH - 2 * L*N);
            a = BBE_MULUSRNX16(d0, x0, d1);
            x0 = BBE_PACKVNX40(a, d1);
            x0 = BBE_MULNX16PACKQ(k_norm, x0);
            BBE_SVNX16_XP(x0, pv, 2 * BBE_SIMD_WIDTH*(size_V - 1));
        }
    }
    else
    {
        for (p = 0; p < L; p += BBE_SIMD_WIDTH)
        {
            x0 = BBE_LVNX16_I(x, 0);
            a = 0;
            __Pragma("loop_count min=2")
                for (i = 0; i < len_x - 1; i++)
                {
                    BBE_LVNX16_XP(x0, x, 2 * L*N);
                    BBE_MULANX16(a, x0, x0);
                }
            BBE_LVNX16_XP(x0, x, -2 * L*N*(len_x - 1));
            BBE_MULANX16(a, x0, x0);

            a = BBE_ADDNX40(a, a);
            d1 = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, d1);
            BBE_RSQRTLUNX40_0(b, t, x0, a);
            BBE_RSQRTLUNX40_1(b, t, x0, a);
            BBE_MULUUSNX16(b, x0, t);
            a = BBE_SRAINX40(b, 23);
            d0 = BBE_PACKLNX40(a);
            d1 = BBE_SUBSR1SAVSN(16 + 4, d1);
            BBE_LVNX16_XP(x0, x, 2 * L*N);
            a = BBE_MULUSRNX16(d0, x0, d1);
            norm_x = BBE_PACKVNX40(a, d1);
            sign = BBE_LTNX16(x0, _0);
            BBE_ADDSNX16F(norm_x, norm_x, _0x4000, sign);
            BBE_SUBSNX16T(norm_x, norm_x, _0x4000, sign);
            x0 = BBE_ABSSNX16(norm_x);
            x0 = BBE_ADDNX16(x0, x0);
            a = BBE_UNPKUNX16(x0);

            exp = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, exp);
            exp = BBE_SUBSR1SAVSN(35, exp);
            BBE_RSQRTLUNX40_0(b, t, x0, a);
            BBE_RSQRTLUNX40_1(b, t, x0, a);
            BBE_MULUUSNX16(b, x0, t);
            k_norm = BBE_PACKVNX40(b, exp);
            x0 = BBE_MULNX16PACKQ(k_norm, norm_x);
            BBE_SVNX16_IP(x0, pv, 2 * BBE_SIMD_WIDTH);

            __Pragma("loop_count min=1")
                for (i = 1; i < len_x - 1; i++)
                {
                    BBE_LVNX16_XP(x0, x, 2 * L*N);
                    a = BBE_MULUSRNX16(d0, x0, d1);
                    x0 = BBE_PACKVNX40(a, d1);
                    x0 = BBE_MULNX16PACKQ(k_norm, x0);
                    BBE_SVNX16_IP(x0, pv, 2 * BBE_SIMD_WIDTH);
                }
            BBE_LVNX16_XP(x0, x, 2 * BBE_SIMD_WIDTH - 2 * L*N*(len_x - 1));
            a = BBE_MULUSRNX16(d0, x0, d1);
            x0 = BBE_PACKVNX40(a, d1);
            x0 = BBE_MULNX16PACKQ(k_norm, x0);
            BBE_SVNX16_XP(x0, pv, 2 * BBE_SIMD_WIDTH*(size_V - len_x + 1));
        }
    }
}


/*
    update columns of matrix R by Housholder vectors
    */
static void qrUpdateR(int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int len_x, int L)
{
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r, v, vr, _0x4000 = 0x4000;
    int i, l, j, sizeV = (SIZE_OF_V(M, N));
    NASSERT(len_x > 1);
    if (len_x == 2)
    {
        for (j = 0; j < N - m; j++)
        {
            pv = (xb_vecNx16*)(V);
            pr = (xb_vecNx16*)(R);
            pw = (xb_vecNx16*)(R);
            for (l = 0; l < L; l += BBE_SIMD_WIDTH)
            {
                xb_vecNx16 v0, v1, r0, r1;
                r1 = BBE_LVNX16_X(pr, 2 * N*L);
                BBE_LVNX16_IP(r0, pr, 2 * BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(v1, pv, 2 * BBE_SIMD_WIDTH*(sizeV - 1));
                acc = BBE_MULRNX16(r0, v0, q14);
                BBE_MULANX16(acc, r1, v1);
                vr = BBE_PACKVNX40(acc, q14);
                acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
                acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
                BBE_SVNX16_X(r1, pw, 2 * N*L);
                BBE_SVNX16_IP(r0, pw, 2 * BBE_SIMD_WIDTH);
            }
            R = (int16_t*)XT_ADDX2(L, (uintptr_t)R);
        }
    }
    else
    {
        for (j = 0; j < N - m; j++)
        {
            pv = (xb_vecNx16*)(V);
            pr = (xb_vecNx16*)(R);
            for (l = 0; l < L; l += BBE_SIMD_WIDTH)
            {
                pw = (xb_vecNx16*)pr;
                BBE_LVNX16_XP(r, pr, 2 * N*L);
                BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                acc = BBE_MULRNX16(r, v, q14);
                __Pragma("loop_count min=2")
                    for (i = 1; i < len_x; i++)
                    {
                        BBE_LVNX16_XP(r, pr, 2 * N*L);
                        BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                        BBE_MULANX16(acc, r, v);
                    }
                vr = BBE_PACKVNX40(acc, q14);
                pv = (const xb_vecNx16*)XT_ADD(-2 * BBE_SIMD_WIDTH*(len_x), (uintptr_t)pv);
                pr = (const xb_vecNx16*)pw;
                //__Pragma("loop_count min=3") : can not be used here due to the compiler bug
                for (i = 0; i < len_x; i++)
                {
                    BBE_LVNX16_XP(r, pr, 2 * N*L);
                    BBE_LVNX16_IP(v, pv, 2 * BBE_SIMD_WIDTH);
                    acc = BBE_MULRNX16(r, _0x4000, q14);
                    BBE_MULSNX16(acc, v, vr);
                    r = BBE_PACKVNX40(acc, q14);
                    BBE_SVNX16_XP(r, pw, 2 * N*L);
                }
                pv = (const xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH*(sizeV - len_x), (uintptr_t)pv);
                pr = (const xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH - len_x * 2 * N*L, (uintptr_t)pr);
            }
            R = (int16_t*)XT_ADDX2(L, (uintptr_t)R);
        }
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
void  qr_build_rmxns ( void* pScr, int16_t * restrict V, int16_t * restrict R, int M, int N,int L)
{
    int m;
    int Ncolumns = M > N ? N : N - 1;  // update N columns for overdetermined linear equation system
    int16_t * pV = V;
    int16_t* pR;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH == 0);
    NASSERT(M >= N);
    (void)pScr;

    pR = R;
    for (m = 0; m < Ncolumns; m++)
    {
        qrHouseholder(pR, pV, M, N, (M - m), L);
        qrUpdateR(pR, pV, m, M, N, (M - m), L);
        pR += (N + 1)*L;
        pV += (M - m)*BBE_SIMD_WIDTH;
    }
} /* qr_build_rmxns() */
#endif

size_t qr_build_rmxns_getScratchSize (int M, int N,int L)
{
    (void)M;
    (void)N;
    (void)L;
    return 0;
} /* qr_build_rmxns_getScratchSize() */

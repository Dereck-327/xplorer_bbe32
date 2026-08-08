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

// final rotation of matrices R(M,N)xL
static void cqrRotateR7x7(int16_t* R, const int16_t* pFi, int L)
#if 1
{
    vsaN q14 = BBE_MOVVSA32(14);
    int k, i;
    xb_vecNx16  fi, r;
    xb_vecNx40  acc0;
    const xb_vecNx16* restrict _pFi;
    const xb_vecNx16* restrict _pr;
    xb_vecNx16* restrict _pw;

    _pr = (const xb_vecNx16*)(R);
    _pw = (xb_vecNx16*)(R);
    _pFi = (const xb_vecNx16*)pFi;
    for (k = i = 0; i < 7 * (L >> (LOG2_BBE_SIMD_WIDTH - 1)); i++)
    {
        int strideFi, strideR;
        k = BBE_ADDMOD16U(k, (7 << 16) | 1);
        strideFi = 2 * BBE_SIMD_WIDTH;
        strideR = 4 * L;
        XT_MOVEQZ(strideFi, 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(7, 7) + 1), k);
        XT_MOVEQZ(strideR, -(7 * 7 - 1) * 4 * L + 2 * BBE_SIMD_WIDTH, k);
        BBE_LVNX16_XP(fi, _pFi, strideFi);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, 4 * L);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, 4 * L);
        BBE_LVNX16_XP(r, _pr, strideR);
        acc0 = BBE_MULRNX16J(r, fi, q14);
        r = BBE_PACKVNX40(acc0, q14);
        BBE_SVNX16_XP(r, _pw, strideR);
    }
}
#else
{
    int ALIGN(8) temp[4];
    xb_vecNx16 idx,idx0,idx1,zero=0;
    vsaN q14=BBE_MOVVSA32(14);
    int k,i; 
    vboolN b;
    xb_vecNx16  fi, r;
    xb_vecNx40  acc0;
    const xb_vecNx16* restrict _pFi; 
    const xb_vecNx16* restrict _pr; 
    xb_vecNx16* restrict _pw; 

    temp[0]=2*BBE_SIMD_WIDTH;
    temp[1]=4*L;
    temp[2]=2*BBE_SIMD_WIDTH*(SIZE_OF_V(7,7)+1);
    temp[3]=-(7*7-1)*4*L+2*BBE_SIMD_WIDTH;
    idx0=BBE_LV4X16_I(temp,0);
    idx1=BBE_LV4X16_I(temp,8);

    _pr =  (const xb_vecNx16*)(R);
    _pw =  (      xb_vecNx16*)(R);
    _pFi = (const xb_vecNx16*)pFi; 
    for(k=i=0; i<7*(L>>(LOG2_BBE_SIMD_WIDTH-1)); i++)
    {
        int strideFi,strideR;
        k=BBE_ADDMOD16U(k,(7<<16)|1);
        idx=BBE_MOVVA16(k);
        b=BBE_EQNX16(idx,zero);
        idx=BBE_MOVNX16T(idx1,idx0,b);
        strideFi=BBE_EXTRNX16C(idx,0);
        strideR =BBE_EXTRNX16C(idx,1);

        BBE_LVNX16_XP(fi,_pFi, strideFi);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, 4*L);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,4*L);
        BBE_LVNX16_XP(r,_pr, strideR);
        acc0 = BBE_MULRNX16J(r,fi,q14);
        r=BBE_PACKVNX40(acc0,q14);
        BBE_SVNX16_XP(r, _pw,strideR);
    }
}
#endif

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
void cqr_build_r7x7s  ( void* pScr, complex_fract16 * restrict _V, complex_fract16 * restrict _R, int L)
{
    int16_t * restrict V=(int16_t *)_V;
    int16_t * restrict R=(int16_t *)_R;
    int16_t * restrict pV;
    int16_t * restrict pFi;
    int16_t * restrict pR;

    NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(V, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(R, BBE_SIMD_WIDTH * 2);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    pFi = V + BBE_SIMD_WIDTH* SIZE_OF_V(7, 7);  pR = R; pV = V;

    cqrHouseholder7(pScr, pR, pV, pFi, 7, 7, L);
    cqrUpdateR7(pR, pV, 0, 7, 7, L);
    pFi += BBE_SIMD_WIDTH; pV += BBE_SIMD_WIDTH * 7; pR += (7 + 1)*L * 2;
    cqrHouseholder6(pScr, pR, pV, pFi, 7, 7, L);
    cqrUpdateR6(pR, pV, 1, 7, 7, L);
    pFi += BBE_SIMD_WIDTH; pV += BBE_SIMD_WIDTH * 6; pR += (7 + 1)*L * 2;
    cqrHouseholder5(pScr, pR, pV, pFi, 7, 7, L);
    cqrUpdateR5(pR, pV, 2, 7, 7, L);
    pFi += BBE_SIMD_WIDTH; pV += BBE_SIMD_WIDTH * 5; pR += (7 + 1)*L * 2;
    cqrHouseholder4(pR, pV, pFi, 7, 7, L);
    cqrUpdateR4(pR, pV, 3, 7, 7, L);
    pFi += BBE_SIMD_WIDTH; pV += BBE_SIMD_WIDTH * 4; pR += (7 + 1)*L * 2;
    cqrHouseholder3(pR, pV, pFi, 7, 7, L);
    cqrUpdateR3(pR, pV, 4, 7, 7, L);
    pFi += BBE_SIMD_WIDTH; pV += BBE_SIMD_WIDTH * 3; pR += (7 + 1)*L * 2;
    cqrHouseholder2(pR, pV, pFi, 7, 7, L);
    cqrUpdateR2(pR, pV, 5, 7, 7, L);
    // for square matrices we have to calculate last rotation explicitely
    cqrComputeLastRot(V + BBE_SIMD_WIDTH* SIZE_OF_V(7, 7), R, 7, L);
    // final rotation
    pFi = V + BBE_SIMD_WIDTH* SIZE_OF_V(7, 7);
    cqrRotateR7x7(R, pFi, L);
}

size_t cqr_build_r7x7s_getScratchSize (int M, int N,int L)
{
    (void)M; (void)N;
    return 4 * L;
} /* cqr_build_r7x7s_getScratchSize() */
#else
DISCARD_FUN(void, cqr_build_r7x7s, ( void* pScr,complex_fract16 * restrict V, complex_fract16 * restrict R, int L))
size_t cqr_build_r7x7s_getScratchSize (int M, int N,int L) { (void)M;(void)N;(void)L; return 0; }
#endif

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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, real data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "qrsf_common.h"

#if HAVE_VFPU
/*
Reference code:
% compute Q'B matrix
% input:
% V  - sequence of Householder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% R  - upper triangle decomposition
function [B] = cqr_calcQB(B,V,Fi)
[M, P] = size(B); 
[N, t] = size(Fi);
Z=zeros(M,P);
im=1;
for m=1:N
    v=V(im:im+M-m);
    im=im+(M-m+1);
    Bm=B(m:end,:);
    Bm=(Bm-2*v*v'*Bm);
    B(m:end,:)=Bm;
end
B=diag([Fi;ones(M-N,1)])'*B;
*/

/*---------------------------------------------------------
    update matrix B[L][MxP] by housholder vectors V[L][SV]

    Input:
    M,N,P,L     dimensions
    V[L][SV]    Householder vectors
    Input/output:
    B[L][MxP]   B matrices MxP
  
---------------------------------------------------------*/
void qrsfUpdateB(float32_t* B,const float32_t* V,int M,int P,int L)
{
#if 0
    float32_t v_re, z_re, b_re;
    int m, l, p, _PL = P*L;
    float32_t A_re;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 0 && M >= 1 && L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));
    for (p = 0; p < P; p++)
    {
        for (l = 0; l < L; l++)
        {
            // compute v'B first
            A_re = 0.f;
            for (m = 0; m < M; m++)
            {
                v_re = V[l + m*L];
                b_re = B[m*_PL + l];
                A_re += (v_re*b_re);
            }
            z_re = 2.f*A_re;
            for (m = 0; m < M; m++)
            {
                A_re = B[m*_PL + l];
                v_re = V[l + m*L];
                A_re -= (z_re*v_re);
                B[m*_PL + l] = A_re;
            }
        }
        B += L; // next column
    }
#endif // 0

    int m, l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pB0;
    const xb_vecN_2xf32 * restrict pV;
    const xb_vecN_2xf32 * restrict pV0;

    xb_vecN_2xf32 Acc, V0, B0, Z0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 0 && M >= 1 && L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));

    WUR_CBEGIN((uintptr_t)V);
    WUR_CEND((uintptr_t)V + L * sizeof(float32_t));

    pB0 = (xb_vecN_2xf32 *)B;
    pV0 = (const xb_vecN_2xf32 *)V;
    for (l = 0; l < P*L; l += (BBE_SIMD_WIDTH / 2))
    {
        pV = pV0;
        pBr = pB0;
        Acc = BBE_ZERON_2XF32();
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, P*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc, V0, B0);
        }
        pV = pV0;
        pB = pB0;
        pBr = pB0;
        Z0 = BBE_MULN_2XF32(Acc, 2.f);
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, P*L * sizeof(float32_t));
            BBE_MULSN_2XF32(B0, Z0, V0);
            BBE_SVN_2XF32_XP(B0, pB, P*L * sizeof(float32_t));
        }
        BBE_LVN_2XF32_IC(V0, pV0);
        pB0 = (xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pB0);
    }
}


/*
    rotate B[L][NxP] by diagonal matrix Fi'[L][N]
*/
void qrsfRotateBconj(float32_t* B,const float32_t* Fi,int N,int P,int L)
{
#if 0
    float32_t f_re, b_re;
    int n, p, l;
    float32_t A_re;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 0 && N > 1 && L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));
    for (n = 0; n < N; n++)
    {
        for (l = 0; l < L; l++)
        {
            f_re = Fi[l + n*L];
            for (p = 0; p < P; p++)
            {
                b_re = B[l + p*L];
                A_re = (f_re*b_re);
                B[l + p*L] = A_re;
            }
        }
        B += P*L;   // next row
    }
#endif // 0

    int n, l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pF;

    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(P > 0 && N > 1 && L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);

    pBr = (const xb_vecN_2xf32 *)B;
    pB = (xb_vecN_2xf32 *)B;
    pF = (const xb_vecN_2xf32 *)Fi;
    for (n = 0; n < N; n++)
    {
        WUR_CBEGIN((uintptr_t)pF);
        WUR_CEND((uintptr_t)pF + L * sizeof(float32_t));

        for (l = 0; l < P*L; l += (BBE_SIMD_WIDTH / 2))
        {
            BBE_LVN_2XF32_IC(F0, pF);
            BBE_LVN_2XF32_XP(B0, pBr, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULN_2XF32(F0, B0);
            BBE_SVN_2XF32_XP(Acc, pB, 2 * BBE_SIMD_WIDTH);
        }
        pF = (const xb_vecN_2xf32 *)XT_ADDX4(L, (uintptr_t)pF);
    }
}

/*-------------------------------------------------------------------------
[c]qr_calc_qbMxNxPsf

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
[c]qr_build_rMxNsf.

Data format: IEEE-754 Std single precision floating-point

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                Matrices B (L matrices of size MxP)
V[((2*M-N+1)*N/2+N)*L]  L sets of Householder vectors
Output:
B[M*P][L]               Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
void qr_calc_qbmxnxpsf  (void *pScr, float32_t *B, const float32_t* V , int M, int N, int P, int L)
{
    int m;
    const float32_t *pV;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(V,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    for (pV=V,m=0; m<N; m++)
    {
        qrsfUpdateB(B+m*P*L,pV,(M-m),P,L); 
        pV+=(M-m)*L;
    }
    qrsfRotateBconj(B,pV,N,P,L);
}

// scratch memory needed for cqbxxxs functions
size_t qr_calc_qbmxnxpsf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    return 0; 
}
#else
DISCARD_FUN(void, qr_calc_qbmxnxpsf, (void *pScr, float32_t *B, const float32_t* V , int M, int N, int P, int L))

size_t qr_calc_qbmxnxpsf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    return 0; 
}
#endif

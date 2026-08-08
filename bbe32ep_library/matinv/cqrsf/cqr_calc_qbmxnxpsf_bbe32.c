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
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cqrsf_common.h"

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
void cqrsfUpdateB(float32_t* B,const float32_t* V,int M,int P,int L)
{
#if 0
    float32_t v_re, v_im, z_re, z_im, b_re, b_im;
    int m, l, p, _2L = 2 * L, _PL = 2 * P*L;
    float32_t A_re, A_im;
    for (p = 0; p < P; p++)
    {
        for (l = 0; l < L; l++)
        {
            // compute v'B first
            A_re = A_im = 0.f;
            for (m = 0; m < M; m++)
            {
                v_re = V[l * 2 + m * 2 * L + 0];
                v_im = V[l * 2 + m * 2 * L + 1];
                b_re = B[m*_PL + l * 2 + 0];
                b_im = B[m*_PL + l * 2 + 1];
                A_re += (v_re*b_re) + (v_im*b_im);
                A_im += (v_re*b_im) - (v_im*b_re);
            }
            z_re = 2.f*A_re;
            z_im = 2.f*A_im;
            for (m = 0; m < M; m++)
            {
                A_re = B[m*_PL + l * 2 + 0];
                A_im = B[m*_PL + l * 2 + 1];
                v_re = V[l * 2 + m * 2 * L + 0];
                v_im = V[l * 2 + m * 2 * L + 1];
                A_re -= (z_re*v_re) - (z_im*v_im);
                A_im -= (z_re*v_im) + (z_im*v_re);
                B[m*_PL + l * 2 + 0] = A_re;
                B[m*_PL + l * 2 + 1] = A_im;
            }
        }
        B += _2L; // next column
    }
#endif // 0

    int m, l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
          xb_vecN_2xf32 * restrict pB0;
    const xb_vecN_2xf32 * restrict pV;
    const xb_vecN_2xf32 * restrict pV0;

    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, V0, B0, Z0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 4) == 0));

    WUR_CBEGIN((uintptr_t)V);
    WUR_CEND((uintptr_t)V + 2 * L * sizeof(float32_t));

    pB0 = (xb_vecN_2xf32 *)B;
    pV0 = (const xb_vecN_2xf32 *)V;
    for (l = 0; l < P*L; l += (BBE_SIMD_WIDTH / 4))
    {
        pV = pV0;
        pBr = pB0;
        Acc = BBE_ZERON_2XF32();
        Acc1 = BBE_ZERON_2XF32();
        Acc2 = BBE_ZERON_2XF32();
        Acc3 = BBE_ZERON_2XF32();
        for (m = 0; m < M >> 1; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, 2 * L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, 2 * P*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, V0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc1, V0, B0, 2, 11);

            BBE_LVN_2XF32_XP(V0, pV, 2 * L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, 2 * P*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc2, V0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc3, V0, B0, 2, 11);
        }
        if (M & 1)
        {
            BBE_LVN_2XF32_XP(V0, pV, 2 * L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, 2 * P*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, V0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc1, V0, B0, 2, 11);
        }
        pV = pV0;
        pB = pB0;
        pBr = pB0;
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Z0 = BBE_MULN_2XF32(Acc, 2.f);
        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(V0, pV, 2 * L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B0, pBr, 2 * P*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(B0, Z0, V0, 3, 4);
            BBE_MULMASN_2XF32(B0, Z0, V0, 2, 11);
            BBE_SVN_2XF32_XP(B0, pB, 2 * P*L * sizeof(float32_t));
        }
        BBE_LVN_2XF32_IC(V0, pV0);
        pB0 = (xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pB0);
    }
}

/*
    rotate B[L][NxP] by diagonal matrix Fi'[L][N]
*/
static void cqrsfRotateBconj(float32_t* B,const float32_t* Fi,int N,int P,int L)
{
#if 0
    float32_t f_re, f_im, b_re, b_im;
    int n, p, l, _2L = 2 * L;
    float32_t A_re, A_im;

    for (n = 0; n < N; n++)
    {
        for (l = 0; l < L; l++)
        {
            f_re = Fi[l * 2 + n*_2L + 0];
            f_im = Fi[l * 2 + n*_2L + 1];
            for (p = 0; p < P; p++)
            {
                b_re = B[l * 2 + p*_2L + 0];
                b_im = B[l * 2 + p*_2L + 1];
                A_re = (f_re*b_re) + (f_im*b_im);
                A_im = (f_re*b_im) - (f_im*b_re);
                B[l * 2 + p*_2L + 0] = A_re;
                B[l * 2 + p*_2L + 1] = A_im;
            }
        }
        B += 2 * P*L;   // next row
    }
#endif // 0

    int n, l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pF;

    xb_vecN_2xf32 Acc, F0, B0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 4) == 0);

    pBr = (const xb_vecN_2xf32 *)B;
    pB = (xb_vecN_2xf32 *)B;
    pF = (const xb_vecN_2xf32 *)Fi;
    for (n = 0; n < N; n++)
    {
        WUR_CBEGIN((uintptr_t)pF);
        WUR_CEND((uintptr_t)pF + 2 * L * sizeof(float32_t));

        for (l = 0; l < P*L; l += (BBE_SIMD_WIDTH / 4))
        {
            BBE_LVN_2XF32_IC(F0, pF);
            BBE_LVN_2XF32_IP(B0, pBr, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULMN_2XF32(F0, B0, 0, 4);
            BBE_MULMASN_2XF32(Acc, F0, B0, 2, 11);
            BBE_SVN_2XF32_IP(Acc, pB, 2 * BBE_SIMD_WIDTH);
        }
        pF = (const xb_vecN_2xf32 *)XT_ADDX4(2 * L, (uintptr_t)pF);
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
void cqr_calc_qbmxnxpsf  (void *pScr, complex_float *_B, const complex_float *_V , int M, int N, int P, int L)
{
    float32_t *B=(float32_t *)_B;
    float32_t *V=(float32_t *)_V;
    int m;
    const float32_t *pV;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(V,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    for (pV=V,m=0; m<N; m++)
    {
        cqrsfUpdateB(B+2*m*P*L,pV,(M-m),P,L); 
        pV+=2*(M-m)*L;
    }
    cqrsfRotateBconj(B,pV,N,P,L);
}

// scratch memory needed for cqbxxxs functions
size_t cqr_calc_qbmxnxpsf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return 0; 
}
#else
DISCARD_FUN(void,cqr_calc_qbmxnxpsf, (void *pScr, complex_float *_B, const complex_float *_V , int M, int N, int P, int L))
size_t cqr_calc_qbmxnxpsf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return 0; 
}
#endif

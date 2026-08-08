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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*===========================================================================
  complex Cholesky decomposition (streaming format, floating point) 
  C code optimized for BBE32

  Integrit 2006-2017
===========================================================================*/
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex- or real-valued least squares problem: A*X=B, 
where A is an MxN coefficient matrix with M >= N; X is an NxP matrix of
unknowns; and B is an MxP right hand matrix.

The decomposition results in an upper triangular NxN matrix R with real and
positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the (conjugate) transpose of a matrix, and sigma2*I is
an NxN identity matrix multiplied by the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the stream order.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see [r]cholfwdmxnxpsf() and [r]cholbkwnxpsf(), 
respectively.

Data format: IEEE-754 Std. single precision floating-point

Input:
  M, N           Dimensional parameters
  L              Number of matrices
  sigma2[L]      Regularization term
  A[M*N][L]      Sequence of L real/complex matrices A
Output:
  R[N*N][L]      Sequence of L upper triangular real/complex matrices R
  D[L*N]         Reciprocal of main diagonal written in a special format
Restrictions:
  1. A, R, D, sigma2 must not overlap
  2. A, R, D, sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 4 for complex-valued 
     functions, or a multiple of 8 for real-valued functions.
  4. Matrix sizes must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
#if !(HAVE_VFPU)
DISCARD_FUN(void, cholmxnsf,(
            complex_float * restrict _R,
            complex_float       * restrict _D,
            const complex_float * restrict _A, 
            const float32_t     * restrict sigma2,
            int M, int N, int L))
#else
/*
Reference Matlab code:
function R=chol2(A,sigma2)
sz=size(A); M=sz(1); N=sz(2);
R=zeros(N,N);
for m=1:N
    Rm=R(:,m);  % take m-th column of original and decomposing matrix
    Am=A(:,m);
    Amm=Am'*Am+sigma2;
    Rmm=Rm'*Rm;
    Rmm=sqrt(real(Amm-Rmm));
    x(1,1:m)=[zeros(1,m-1) Rmm];
    for k=m+1:N
        Akm=A(:,k)'*Am;
        Rkm=R(:,k)'*Rm;
        x(1,k)=(Akm-Rkm)/Rmm;
    end
    R(m,:)=conj(x);
end
*/

#define VECLEN (BBE_SIMD_WIDTH/4)

static void cleanR(float32_t *R,int M,int N,int L)
{
    int k;
    xb_vecN_2xf32 * restrict pR = (xb_vecN_2xf32 *)R;

    NASSERT_ALIGN(R, (BBE_SIMD_WIDTH * 2));

    for (k = 0; k < (M*N*L / (BBE_SIMD_WIDTH / 4)); k++)
    {
        BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(), pR, 2 * BBE_SIMD_WIDTH);
    }
}

/*
    one iteration of Cholesky (VECLEN matrices in parallel)
    Input:
    A[2*M*N*L]
    Input/Output:
    R[2*N*N*L]  (updated m-th row)

*/
static void iter2(
                const float32_t* restrict A,
                      float32_t* restrict R,
                      float32_t* restrict D,
                const float32_t* restrict sigma2,
                int m,int M,int N,int L)
{
#if 0
    int n,k,l;
    NASSERT_ALIGN(R, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(A, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(D, (BBE_SIMD_WIDTH * 2));
    NASSERT_ALIGN(sigma2, (BBE_SIMD_WIDTH));
    NASSERT(L % (BBE_SIMD_WIDTH / 4) == 0);
    NASSERT(L>0);

    // take colunms of A and R and calculate diagonal elements
    for (l = 0; l<VECLEN; l++)
    {
        float32_t Acc;
        Acc = sigma2[l];
        for (k = 0; k<M; k++)
        {
            float32_t re, im;
            re = A[(l + (m + k*N)*L) * 2 + 0];
            im = A[(l + (m + k*N)*L) * 2 + 1];
            Acc += re*re + im*im;
        }
        for (k = 0; k<m; k++) // NOTE: this loop may be from 0 to N: for (k=0; k<N; k++)
        {
            float32_t re, im;
            re = R[(l + (m + k*N)*L) * 2 + 0];
            im = R[(l + (m + k*N)*L) * 2 + 1];
            Acc -= re*re + im*im;
        }
        if (Acc <= 0.f)
        {
            Acc = 1.f;
        }
        // first calculate 1/sqrt(Acc)
        D[2 * l + 0] = D[2 * l + 1] = 1.f / sqrtf(Acc);
    }
    // find elements in row from m+1 to N
    for (l = 0; l<VECLEN; l++)
    {
        for (k = m; k<N; k++)
        {
            float32_t A_re, A_im;
            A_re = A_im = 0.f;
            A_re = sigma2[l];
            for (n = 0; n<M; n++)
            {
                float32_t ak_re, ak_im;
                float32_t am_re, am_im;
                ak_re = A[(l + (k + n*N)*L) * 2 + 0];
                ak_im = A[(l + (k + n*N)*L) * 2 + 1];
                am_re = A[(l + (m + n*N)*L) * 2 + 0];
                am_im = A[(l + (m + n*N)*L) * 2 + 1];
                A_re += ak_re*am_re + ak_im*am_im;
                A_im -= ak_re*am_im - ak_im*am_re;
            }

            for (n = 0; n<m; n++) // this loop also may be to the N: for (n=0; n<N; n++)
            {
                float32_t rk_re, rk_im;
                float32_t rm_re, rm_im;
                rk_re = R[(l + (k + n*N)*L) * 2 + 0];
                rk_im = R[(l + (k + n*N)*L) * 2 + 1];
                rm_re = R[(l + (m + n*N)*L) * 2 + 0];
                rm_im = R[(l + (m + n*N)*L) * 2 + 1];
                A_re -= rk_re*rm_re + rk_im*rm_im;
                A_im += rk_re*rm_im - rk_im*rm_re;
            }
            R[(l + (k + m*N)*L) * 2 + 0] = A_re*D[2 * l + 0];
            R[(l + (k + m*N)*L) * 2 + 1] = A_im*D[2 * l + 1];
        }
    }
#endif

    int n, k;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pAk;
    const xb_vecN_2xf32 * restrict pAm;
          xb_vecN_2xf32 * restrict pRk;
    const xb_vecN_2xf32 * restrict pRm;
    const xb_vecN_2xf32 * restrict pS;

    valign vs;
    vboolN_2 vpred;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, rsqrt, sigma;
    xb_vecN_2xf32 A0, A1;
    xb_vecN_2xf32 R0, R1;

    NASSERT_ALIGN(R,(BBE_SIMD_WIDTH*2));
    NASSERT_ALIGN(A,(BBE_SIMD_WIDTH*2));
    NASSERT_ALIGN(D,(BBE_SIMD_WIDTH*2));
    NASSERT_ALIGN(sigma2,(BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(L > 0);

    // take colunms of A and R and calculate diagonal elements
    pD = (xb_vecN_2xf32 *)D;
    pAm = (const xb_vecN_2xf32 *)(A + 2 * m*L);
    pRm = (const xb_vecN_2xf32 *)(R + 2 * m*L);
    pS = (const xb_vecN_2xf32 *)sigma2;
    vs = BBE_LAN_2XF32_PP(pS);

    BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
    Acc = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);
    for (k = 0; k < M; k++)
    {
        BBE_LVN_2XF32_XP(A0, pAm, 2 * N*L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A0);
    }
    for (k = 0; k < m; k++) // NOTE: this loop may be from 0 to N: for (k=0; k<N; k++)
    {
        BBE_LVN_2XF32_XP(R0, pRm, 2 * N*L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc, R0, R0);
    }
    A0 = BBE_SHFLN_2XF32I(Acc, BBE_SHFLI_SWAP_2);
    Acc = BBE_ADDN_2XF32(Acc, A0);

    vpred = BBE_OGTN_2XF32(Acc, BBE_ZERON_2XF32());
    rsqrt = 1.f;
    BBE_RSQRTN_2XF32T(rsqrt, Acc, vpred);
    BBE_SVN_2XF32_IP(rsqrt, pD, 2 * BBE_SIMD_WIDTH);

    // find elements in row from m+1 to N
    for (k = m; k < N; k++)
    {
        Acc = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);
        Acc1 = BBE_ZERON_2XF32();
        Acc2 = BBE_ZERON_2XF32();
        Acc3 = BBE_ZERON_2XF32();
        pAm = (const xb_vecN_2xf32 *)(A + 2 * m*L);
        pRm = (const xb_vecN_2xf32 *)(R + 2 * m*L);
        pAk = (const xb_vecN_2xf32 *)(A + 2 * k*L);
        pRk = (      xb_vecN_2xf32 *)(R + 2 * k*L);
        for (n = 0; n < M >> 1; n++)
        {
            BBE_LVN_2XF32_XP(A0, pAm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(A1, pAk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc1, A0, A1, 2, 11);
            BBE_LVN_2XF32_XP(A0, pAm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(A1, pAk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc2, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        }
        if (M & 1)
        {
            BBE_LVN_2XF32_XP(A0, pAm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(A1, pAk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc2, A0, A1, 0, 4);
            BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        }
        for (n = 0; n < m >> 1; n++) // this loop also may be to the N: for (n=0; n<N; n++)
        {
            BBE_LVN_2XF32_XP(R0, pRm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(R1, pRk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, R0, R1, 3, 4);
            BBE_MULMASN_2XF32(Acc1, R0, R1, 1, 11);
            BBE_LVN_2XF32_XP(R0, pRm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(R1, pRk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc2, R0, R1, 3, 4);
            BBE_MULMASN_2XF32(Acc3, R0, R1, 1, 11);
        }
        if (m & 1)
        {
            BBE_LVN_2XF32_XP(R0, pRm, 2 * N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(R1, pRk, 2 * N*L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc, R0, R1, 3, 4);
            BBE_MULMASN_2XF32(Acc1, R0, R1, 1, 11);
        }
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }
}

void cholmxnsf(
            complex_float * restrict _R,
            complex_float       * restrict _D,
            const complex_float * restrict _A, 
            const float32_t     * restrict sigma2,
            int M, int N, int L)
{
    float32_t *       restrict R=(float32_t *      )_R;
    float32_t *       restrict D=(float32_t *      )_D;
    const float32_t * restrict A=(const float32_t *)_A;
    int l,m;
    NASSERT_ALIGN(R     ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A     ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D     ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2,(2*BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(L>0);
    NASSERT(M>=N);
    
    cleanR(R,N,N,L);
    for (l=0; l<L; l+=VECLEN)
    {
        for (m=0; m<N; m++)
        {
            iter2(A,R,D,sigma2,m,M,N,L);
            D+=2*VECLEN;
        }
        R+=2*VECLEN;
        A+=2*VECLEN;
        sigma2+=VECLEN;
    }
}
#endif

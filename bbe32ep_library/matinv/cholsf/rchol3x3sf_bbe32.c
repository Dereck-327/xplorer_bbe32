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
  real Cholesky decomposition (streaming format, floating point) 
  and related routines:
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
DISCARD_FUN(void, rchol3x3sf,(
                  float32_t * restrict R,
                  float32_t * restrict D,
            const float32_t * restrict A, 
            const float32_t * restrict sigma2,
            int L))
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

#define VECLEN (BBE_SIMD_WIDTH/2)

void rchol3x3sf(
                  float32_t * restrict R,
                  float32_t * restrict D,
            const float32_t * restrict A, 
            const float32_t * restrict sigma2,
            int L)
{
    int l;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pAk;
    const xb_vecN_2xf32 * restrict pA0;
    const xb_vecN_2xf32 * restrict pAm;
          xb_vecN_2xf32 * restrict pRk;
          xb_vecN_2xf32 * restrict pR = (xb_vecN_2xf32 *)R;
    const xb_vecN_2xf32 * restrict pRm;
    const xb_vecN_2xf32 * restrict pS;

 vboolN_2 vpred;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, rsqrt, sigma;
    xb_vecN_2xf32 A0, A1, A2;
    xb_vecN_2xf32 R0, R1;

    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2, (2 * BBE_SIMD_WIDTH));
    NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);
    if (L <= 0) return;

    for (l = 0; l < (3 * 3 * L / (BBE_SIMD_WIDTH / 2)); l++)
    {
        BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(), pR, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)D;
    pS = (const xb_vecN_2xf32 *)sigma2;
    pAm = (const xb_vecN_2xf32 *)(A);
    pRk = (xb_vecN_2xf32 *)(R);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(sigma, pS, 2 * BBE_SIMD_WIDTH);
        Acc = sigma;

        A0 = BBE_LVN_2XF32_X(pAm, 6 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 3 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A0);
        BBE_LVN_2XF32_XP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, A0);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);

        vpred = BBE_OGTN_2XF32(Acc, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 6 * BBE_SIMD_WIDTH);

        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)D;
    pS = (const xb_vecN_2xf32 *)sigma2;
    pAm = (const xb_vecN_2xf32 *)(A);
    pAk = (const xb_vecN_2xf32 *)(A + L);
    pA0 = (const xb_vecN_2xf32 *)(A + 2 * L);
    pRk = (xb_vecN_2xf32 *)(R + L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(sigma, pS, 2 * BBE_SIMD_WIDTH);
        Acc = Acc2 = sigma;

        A0 = BBE_LVN_2XF32_X(pAm, 6 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 6 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA0, 6 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A1);
        Acc3 = BBE_MULN_2XF32(A0, A2);
        A0 = BBE_LVN_2XF32_X(pAm, 3 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 3 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA0, 3 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A1);
        BBE_MULAN_2XF32(Acc2, A0, A2);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A1, pAk, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A2, pA0, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, A1);
        BBE_MULAN_2XF32(Acc3, A0, A2);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);

        BBE_LVN_2XF32_IP(rsqrt, pD, 6 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        Acc2 = BBE_MULN_2XF32(Acc2, rsqrt);

        BBE_SVN_2XF32_X(Acc2, pRk, L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)D + 1;
    pS = (const xb_vecN_2xf32 *)sigma2;
    pRm = (const xb_vecN_2xf32 *)(R + L);
    pRk = (xb_vecN_2xf32 *)(R + 4 * L);
    pAm = (const xb_vecN_2xf32 *)(A + L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(Acc, pS, 2 * BBE_SIMD_WIDTH);

        A0 = BBE_LVN_2XF32_X(pAm, 2 * 3 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 3 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A0);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, A0);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);

        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc, R0, R0);

        vpred = BBE_OGTN_2XF32(Acc, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 6 * BBE_SIMD_WIDTH);

        // find elements in row from m+1 to N
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)D + 1;
    pS = (const xb_vecN_2xf32 *)sigma2;
    pRm = (const xb_vecN_2xf32 *)(R + L);
    pRk = (xb_vecN_2xf32 *)(R + 5 * L);
    pAm = (const xb_vecN_2xf32 *)(A + L);
    pAk = (const xb_vecN_2xf32 *)(A + 2 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(Acc, pS, 2 * BBE_SIMD_WIDTH);

        A0 = BBE_LVN_2XF32_X(pAm, 2 * 3 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 2 * 3 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A1);
        A0 = BBE_LVN_2XF32_X(pAm, 3 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 3 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A1);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A1, pAk, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, A1);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);

        R1 = BBE_LVN_2XF32_X(pRm, L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc, R0, R1);

        BBE_LVN_2XF32_IP(rsqrt, pD, 6 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)D + 2;
    pS = (const xb_vecN_2xf32 *)sigma2;
    pRk = (xb_vecN_2xf32 *)(R + 8 * L);
    pRm = (const xb_vecN_2xf32 *)(R + 2 * L);
    pAm = (const xb_vecN_2xf32 *)(A + 2 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(Acc, pS, 2 * BBE_SIMD_WIDTH);

        A0 = BBE_LVN_2XF32_X(pAm, 2 * 3 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 3 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, A0, A0);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, A0);

        R0 = BBE_LVN_2XF32_X(pRm, 3 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc, R0, R0);
        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc1, R0, R0);

        Acc = BBE_ADDN_2XF32(Acc, Acc1);

        vpred = BBE_OGTN_2XF32(Acc, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 6 * BBE_SIMD_WIDTH);

        Acc = BBE_MULN_2XF32(Acc, rsqrt);
        BBE_SVN_2XF32_IP(Acc, pRk, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

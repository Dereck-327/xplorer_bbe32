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
DISCARD_FUN(void, chol4x4sf,(
            complex_float * restrict _R,
            complex_float       * restrict _D,
            const complex_float * restrict _A, 
            const float32_t     * restrict sigma2,
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

#define VECLEN (BBE_SIMD_WIDTH/4)

void chol4x4sf(
            complex_float * restrict _R,
            complex_float       * restrict _D,
            const complex_float * restrict _A, 
            const float32_t     * restrict sigma2,
            int L)
{
    int l;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pAk;
    const xb_vecN_2xf32 * restrict pA1;
    const xb_vecN_2xf32 * restrict pA2;
    const xb_vecN_2xf32 * restrict pA3;
    const xb_vecN_2xf32 * restrict pAm;
          xb_vecN_2xf32 * restrict pRk;
          xb_vecN_2xf32 * restrict pR = (xb_vecN_2xf32 *)_R;
    const xb_vecN_2xf32 * restrict pRm;
    const xb_vecN_2xf32 * restrict pS;

    valign vs;
  vboolN_2 vpred;
    xb_vecN_2xf32 Acc0, Acc1, Acc2, Acc3, Acc4, Acc5, rsqrt, sigma;
    xb_vecN_2xf32 A0, A1, A2, A3;
    xb_vecN_2xf32 R0, R1, R2;

    NASSERT_ALIGN(_R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_D, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2, (2 * BBE_SIMD_WIDTH));
    NASSERT(L % (BBE_SIMD_WIDTH / 4) == 0);
    if (L<=0) return;

    for (l = 0; l < (4 * 4 * L / (BBE_SIMD_WIDTH / 4)); l++)
    {
        BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(), pR, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)_D;
    pAm = (const xb_vecN_2xf32 *)(_A);
    pRk = (xb_vecN_2xf32 *)(_R);
    for (l = 0; l < L; l += VECLEN)
    {
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);
        
        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc0, A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc1, A0, A0);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc0, A0, A0);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);

        A0 = BBE_SHFLN_2XF32I(Acc0, BBE_SHFLI_SWAP_2);
        Acc0 = BBE_ADDN_2XF32(Acc0, A0);

        vpred = BBE_OGTN_2XF32(Acc0, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc0, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);

        // find elements in row from m+1 to N
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), Acc0, BBE_SELI_INTERLEAVE_2_EVENODD);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        sigma2 += VECLEN;
    }
    sigma2 -= L;

    pD = (xb_vecN_2xf32 *)_D;
    pS = (const xb_vecN_2xf32 *)sigma2;
    vs = BBE_LAN_2XF32_PP(pS);
    pAm = (const xb_vecN_2xf32 *)(_A);
    pA1 = (const xb_vecN_2xf32 *)(_A + L);
    pA2 = (const xb_vecN_2xf32 *)(_A + 2 * L);
    pA3 = (const xb_vecN_2xf32 *)(_A + 3 * L);
    pRk = (xb_vecN_2xf32 *)(_R + L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc3 = Acc4 = Acc5 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 6 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 6 * 4 * L * sizeof(float32_t));
        A3 = BBE_LVN_2XF32_X(pA3, 6 * 4 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        Acc1 = BBE_MULMN_2XF32(A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc4, A0, A2, 2, 11);
        Acc2 = BBE_MULMN_2XF32(A0, A3, 0, 4);
        BBE_MULMASN_2XF32(Acc5, A0, A3, 2, 11);

        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 4 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 4 * 4 * L * sizeof(float32_t));
        A3 = BBE_LVN_2XF32_X(pA3, 4 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc4, A0, A2, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A0, A3, 0, 4);
        BBE_MULMASN_2XF32(Acc5, A0, A3, 2, 11);

        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 2 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 2 * 4 * L * sizeof(float32_t));
        A3 = BBE_LVN_2XF32_X(pA3, 2 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc4, A0, A2, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A0, A3, 0, 4);
        BBE_MULMASN_2XF32(Acc5, A0, A3, 2, 11);

        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A1, pA1, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A2, pA2, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A3, pA3, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc4, A0, A2, 2, 11);
        BBE_MULMASN_2XF32(Acc2, A0, A3, 0, 4);
        BBE_MULMASN_2XF32(Acc5, A0, A3, 2, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc3);
        Acc1 = BBE_ADDN_2XF32(Acc1, Acc4);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc5);
        BBE_LVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        Acc1 = BBE_MULN_2XF32(Acc1, rsqrt);
        Acc2 = BBE_MULN_2XF32(Acc2, rsqrt);

        BBE_SVN_2XF32_X(Acc2, pRk, 4 * L * sizeof(float32_t));
        BBE_SVN_2XF32_X(Acc1, pRk, 2 * L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);
    }

    pD = (xb_vecN_2xf32 *)_D + 1;
    pAm = (const xb_vecN_2xf32 *)(_A + L);
    pRm = (const xb_vecN_2xf32 *)(_R + L);
    pRk = (xb_vecN_2xf32 *)(_R + 5 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc0, A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc1, A0, A0);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc0, A0, A0);

        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc1, R0, R0);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);

        A0 = BBE_SHFLN_2XF32I(Acc0, BBE_SHFLI_SWAP_2);
        Acc0 = BBE_ADDN_2XF32(Acc0, A0);

        vpred = BBE_OGTN_2XF32(Acc0, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc0, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);

        // find elements in row from m+1 to N
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), Acc0, BBE_SELI_INTERLEAVE_2_EVENODD);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        sigma2 += VECLEN;
    }
    sigma2 -= L;

    pD = (xb_vecN_2xf32 *)_D + 1;
    pAm = (const xb_vecN_2xf32 *)(_A + L);
    pA1 = (const xb_vecN_2xf32 *)(_A + 2 * L);
    pA2 = (const xb_vecN_2xf32 *)(_A + 3 * L);
    pRm = (const xb_vecN_2xf32 *)(_R + L);
    pRk = (xb_vecN_2xf32 *)(_R + 6 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc2 = Acc3 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 6 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 6 * 4 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A0, A1, 2, 11);
        Acc1 = BBE_MULMN_2XF32(A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A2, 2, 11);

        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 4 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 4 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A2, 2, 11);

        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pA1, 2 * 4 * L * sizeof(float32_t));
        A2 = BBE_LVN_2XF32_X(pA2, 2 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A2, 2, 11);

        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A1, pA1, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A2, pA2, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc2, A0, A1, 2, 11);
        BBE_MULMASN_2XF32(Acc1, A0, A2, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A2, 2, 11);

        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        R1 = BBE_LVN_2XF32_X(pRk, -2 * 4 * L * sizeof(float32_t));
        R2 = BBE_LVN_2XF32_X(pRk, -2 * 3 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, R0, R1, 3, 4);
        BBE_MULMASN_2XF32(Acc2, R0, R1, 1, 11);
        BBE_MULMASN_2XF32(Acc1, R0, R2, 3, 4);
        BBE_MULMASN_2XF32(Acc3, R0, R2, 1, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
        Acc1 = BBE_ADDN_2XF32(Acc1, Acc3);
        BBE_LVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        Acc1 = BBE_MULN_2XF32(Acc1, rsqrt);

        BBE_SVN_2XF32_X(Acc1, pRk, 2 * L * sizeof(float32_t));
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        sigma2 += VECLEN;
    }
    sigma2 -= L;


    pD = (xb_vecN_2xf32 *)_D + 2;
    pAm = (const xb_vecN_2xf32 *)(_A + 2 * L);
    pRm = (const xb_vecN_2xf32 *)(_R + 2 * L);
    pRk = (xb_vecN_2xf32 *)(_R + 10 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc0, A0, A0);
        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc1, A0, A0);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc0, A0, A0);

        R0 = BBE_LVN_2XF32_X(pRm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc1, R0, R0);
        BBE_LVN_2XF32_IP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc0, R0, R0);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);

        A0 = BBE_SHFLN_2XF32I(Acc0, BBE_SHFLI_SWAP_2);
        Acc0 = BBE_ADDN_2XF32(Acc0, A0);

        vpred = BBE_OGTN_2XF32(Acc0, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc0, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);

        // find elements in row from m+1 to N
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), Acc0, BBE_SELI_INTERLEAVE_2_EVENODD);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        sigma2 += VECLEN;
    }
    sigma2 -= L;

    pD = (xb_vecN_2xf32 *)_D + 2;
    pAm = (const xb_vecN_2xf32 *)(_A + 2 * L);
    pAk = (const xb_vecN_2xf32 *)(_A + 3 * L);
    pRm = (const xb_vecN_2xf32 *)(_R + 6 * L);
    pRk = (xb_vecN_2xf32 *)(_R + 11 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);
        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc3 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        A0 = BBE_LVN_2XF32_X(pAm, 6 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 6 * 4 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, A1, 0, 4);
        Acc1 = BBE_MULMN_2XF32(A0, A1, 2, 11);
        A0 = BBE_LVN_2XF32_X(pAm, 4 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 4 * 4 * L * sizeof(float32_t));
        Acc2 = BBE_MULMN_2XF32(A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
        A0 = BBE_LVN_2XF32_X(pAm, 2 * 4 * L * sizeof(float32_t));
        A1 = BBE_LVN_2XF32_X(pAk, 2 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A0, A1, 2, 11);
        BBE_LVN_2XF32_IP(A0, pAm, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(A1, pAk, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, A0, A1, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);

        R0 = BBE_LVN_2XF32_X(pRm, -2 * 4 * L * sizeof(float32_t));
        R1 = BBE_LVN_2XF32_X(pRk, -4 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, R0, R1, 3, 4);
        BBE_MULMASN_2XF32(Acc1, R0, R1, 1, 11);
        BBE_LVN_2XF32_XP(R0, pRm, 2 * BBE_SIMD_WIDTH);
        R1 = BBE_LVN_2XF32_X(pRk, -2 * 4 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc2, R0, R1, 3, 4);
        BBE_MULMASN_2XF32(Acc3, R0, R1, 1, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);

        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        sigma2 += VECLEN;
    }
    sigma2 -= L;

    pD = (xb_vecN_2xf32 *)_D + 3;
    pRk = (xb_vecN_2xf32 *)(_R + 15 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        pAm = (const xb_vecN_2xf32 *)(_A + 3 * L);
        pRm = (const xb_vecN_2xf32 *)(_R + 3 * L);
        pS = (const xb_vecN_2xf32 *)sigma2;
        vs = BBE_LAN_2XF32_PP(pS);

        BBE_LAVN_2XF32_XP(sigma, vs, pS, BBE_SIMD_WIDTH);
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), sigma, BBE_SELI_INTERLEAVE_2_LO);

        BBE_LVN_2XF32_XP(A0, pAm, 2 * 4 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(A0, A0);
        BBE_LVN_2XF32_XP(A0, pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc0, A0, A0);
        BBE_LVN_2XF32_XP(A0, pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc1, A0, A0);
        BBE_LVN_2XF32_XP(A0, pAm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc0, A0, A0);

        BBE_LVN_2XF32_XP(R0, pRm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc1, R0, R0);
        BBE_LVN_2XF32_XP(R0, pRm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc0, R0, R0);
        BBE_LVN_2XF32_XP(R0, pRm, 2 * 4 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc1, R0, R0);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);

        A0 = BBE_SHFLN_2XF32I(Acc0, BBE_SHFLI_SWAP_2);
        Acc0 = BBE_ADDN_2XF32(Acc0, A0);

        vpred = BBE_OGTN_2XF32(Acc0, BBE_ZERON_2XF32());
        rsqrt = 1.f;
        BBE_RSQRTN_2XF32T(rsqrt, Acc0, vpred);
        BBE_SVN_2XF32_IP(rsqrt, pD, 8 * BBE_SIMD_WIDTH);

        // find elements in row from m+1 to N
        Acc0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), Acc0, BBE_SELI_INTERLEAVE_2_EVENODD);
        Acc0 = BBE_MULN_2XF32(Acc0, rsqrt);
        BBE_SVN_2XF32_IP(Acc0, pRk, 2 * BBE_SIMD_WIDTH);

        _R += VECLEN;
        _A += VECLEN;
        sigma2 += VECLEN;
    }
}
#endif

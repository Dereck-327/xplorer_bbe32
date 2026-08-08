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
  NatureDSP_Baseband library. Cholesky decomposition for a complex-valued pseudo-inversion:
    Apply the Cholesky decomposition to the matrix of normal equations system
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, cholmxns,(
                  complex_fract16 * restrict R,
                  complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t* restrict sigma2,
            int M, int N, int L))
#else

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the streaming order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxps() and cholbkwnxps(), respectively.

Input:
  M, N           Dimensional parameters
  L              Number of matrices
  sigma2[L]      Regularization term; fixed point position is twice the
                 number of fractional bits for matrices A, R
  A[M*N][L]      sequence of L complex matrices A
Output:
  R[N*N][L]      Sequence of L upper triangular complex matrices R
  D[L/4][N][8]   Reciprocal of main diagonal (mantissa, exponent) in the 
                 special format
Restrictions:
  1. A, R, D, sigma2 must not overlap
  2. A, R, D, sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 8
  4. Matrix sizes must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/

void cholmxns (
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
            const int32_t * restrict sigma2,
            int M,int N,
            int L)
{
    int16_t *       restrict R=(int16_t *      )_R;
    int16_t *       restrict D=(int16_t *      )_D;
    const int16_t * restrict A=(const int16_t *)_A;

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
          xb_vecNx16* restrict pD;
          xb_vecNx16* restrict pRm;
          xb_vecNx16* restrict pRk;
    const xb_vecNx16* restrict pAm;
    const xb_vecNx16* restrict pAk;
    const xb_vecNx16* restrict pSigma;
    xb_vecNx16 t,lo,hi,AK,AM,RK,RM,D0;
    xb_vecNx40 wA,wB;
    vsaN sh,D1,_16=BBE_MOVVSA32(16);
    int n,k,l,m;
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2,2*BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    NASSERT(L>0);
    NASSERT(M>=N);
    
    // clean R matrix 
    RK=0;
    pRk=(xb_vecNx16*)R;
    for (k=0; k<(2*N*N*L)/BBE_SIMD_WIDTH; k++) BBE_SVNX16_IP(RK,pRk,2*BBE_SIMD_WIDTH);
    pD=(       xb_vecNx16*)(D);
    pSigma=(const xb_vecNx16*)sigma2;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        for (m=0; m<N; m++)
        {
            // take colunms of A and R and calculate diagonal elements
            pAk=(const xb_vecNx16*)(A);
            pRk=(      xb_vecNx16*)(R);
            BBE_LVNX16_IP(t,pSigma,0);
            BBE_DSELNX16I(hi,lo,t,t,BBE_DSELI_INTERLEAVE_2);
            wA=BBE_MOVSWV(hi,lo);
            wB=BBE_MOVQINT40(0);
            #ifdef COMPILER_XTENSA
            #pragma loop_count min=2
            #endif
            for (k=0; k<M; k++)
            {
                BBE_LVNX16_XP(AM,pAk,4*N*L);
                BBE_MULANX16J(wB,AM,AM);
            }
            for (k=0; k<m; k++)
            {
                BBE_LVNX16_XP(RM,pRk,4*N*L);
                BBE_MULSNX16J(wB,RM,RM);
            }
            wA=BBE_ADDNX40(wA,wB);
            wA=BBE_ADDNX40(wA,wA);
            D1=BBE_NSAENX40(wA);
            wA=BBE_SLLNX40(wA,D1);
            BBE_RSQRTLUNX40_0(wA,t, AK, wA);
            BBE_MULUUSNX16(wA, AK,  t);
            wA=BBE_SRAINX40(wA,23);
            D0=BBE_PACKLNX40(wA);
            D0=BBE_SELNX16I(D0,D0,BBE_SELI_INTERLEAVE_1_EVEN);
            D1=BBE_ADDSAVSN(-16,D1);
            t=BBE_MOVVVS(D1);
            t=BBE_SELNX16I(t,t,BBE_SELI_INTERLEAVE_1_EVEN);
            D1=BBE_MOVVSV(t,1);
            BBE_SVNX16_IP(D0,pD,2*BBE_SIMD_WIDTH);
            t=BBE_MOVVVS(D1);
            BBE_SVNX16_IP(t ,pD,2*BBE_SIMD_WIDTH);
            // compute diagonal element
            sh=BBE_SUBSAVSN(11,D1);
            t=BBE_PACKVNX40(wB,sh);
            wA=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(wA);
            BBE_SVNX16_XP(t,pRk,4*(L-m*N*L));
            // find elements in row from m+1 to N
            pAk=(const xb_vecNx16*)XT_ADDX4(L-M*N*L,(uintptr_t)pAk);
            #ifdef COMPILER_XTENSA
            #pragma concurrent
            #endif
            for (k=1; k<N-m; k++)
            {
                pAm=(const xb_vecNx16*)(A);
                pRm=(      xb_vecNx16*)(R);
                wA=BBE_MOVQINT40(0);
                #ifdef COMPILER_XTENSA
                #pragma loop_count min=2
                #endif
                for (n=0; n<M; n++)
                {
                    BBE_LVNX16_XP(AK,pAk,4*N*L);
                    BBE_LVNX16_XP(AM,pAm,4*N*L);
                    BBE_MULANX16J(wA,AK,AM);
                }
                pAk=(const xb_vecNx16*)XT_ADDX4(L-M*N*L,(uintptr_t)pAk);
                for (n=0; n<m; n++)
                {
                    BBE_LVNX16_XP(RK,pRk,4*N*L);
                    BBE_LVNX16_XP(RM,pRm,4*N*L);
                    BBE_MULSNX16J(wA,RK,RM);
                }
                wA=BBE_SLSNX40(wA,D1);
                lo=BBE_PACKLNX40(wA);
                hi=BBE_PACKVNX40(wA,_16);
                wA=BBE_MULUUNX16(D0,lo);
                wA=BBE_SRAINX40(wA,16);
                BBE_MULUSANX16(wA,D0,hi);
                RK=BBE_PACKPNX40(wA);
                BBE_SVNX16_XP(RK,pRk,4*(L-m*N*L));
            }
            A=(const int16_t*)XT_ADDX4(L,(uintptr_t)A);
            R=(      int16_t*)XT_ADDX4(L,(uintptr_t)R);
        }
        A=(const int16_t*)XT_ADDX4(BBE_SIMD_WIDTH/2-N*L,(uintptr_t)A);
        R=(      int16_t*)XT_ADDX4(BBE_SIMD_WIDTH/2-N*L,(uintptr_t)R);
        BBE_LVNX16_IP(t,pSigma,2*BBE_SIMD_WIDTH);
    }
} /* cholmxns() */
#endif

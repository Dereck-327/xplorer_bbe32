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
  Cholesky forward recursion for pseudo-inversion API (real floating point 
  data)
  C code optimized for BBE32

  Integrit 2006-2017
===========================================================================*/
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Data format: IEEE-754 Std. single precision floating-point

Input:
M               Matrix dimension (number of rows in matrices A)
N               Matrix dimension (number of columns and rows in matrices 
                R)
P               Number of columns in right-side matrices B
L               Number of matrices
R[N*N][L]       Cholesky upper triangular matrices R
A[M*N][L]       Original left-side matrices A
B[M*P][L]       Original right-side matrices B
D[L*N]          Reciprocal of main diagonal written in a special format
Output:
y[N*P][L]       Decision matrix y

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex-valued 
   functions, or a multiple of 8 for real-valued functions.
3. Matrix sizes (M,N) must be greater than 1, P must be >=1
4. M >= N
---------------------------------------------------------------------------*/
#if !(HAVE_VFPU)
DISCARD_FUN(void, rcholfwdmxnxpsf,(
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
            int M, int N, int P, int L))
#else
/*
  Reference Matlab code
  function [Y]=cholfwd(R,A,B)
  sz=size(A); M=sz(1); N=sz(2); 
  sz=size(B); P=sz(2); 
  D=real(1./diag(R));
  AB=A'*B;
  Y=zeros(N,P);
  for n=1:N
    Rn=R(:,n); 
    Bn=AB(n,:);
    y=(Bn-Rn'*Y)*D(n);
    Y(n,:)=y;
  end
*/

#define VECLEN (BBE_SIMD_WIDTH/2)

static void iter2( 
                float32_t* restrict _Y, 
          const float32_t* restrict _R, 
          const float32_t* restrict _A, 
          const float32_t* restrict _B, 
          const float32_t* restrict _D, 
          int n,int M,int N,int P,int L)
{
#if 0
    int k, m, p;
    NASSERT_ALIGN(Y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
    NASSERT(M>0 && P>0 && L>0);
    // calculate A(:,n)'*B-Rn'*Y, 1xP
    for (k = 0; k<VECLEN; k++)
    {
        for (p = 0; p<P; p++)
        {
            float32_t A_re;
            float32_t B_re = 0.f;
            float32_t C_re = 0.f;
            for (m = 0; m<M; m++)
            {
                float32_t a_re, b_re;
                a_re = A[(k + (n + m*N)*L)];
                b_re = B[(k + (p + m*P)*L)];
                B_re += (a_re*b_re);
            }
            for (m = 0; m<n; m++)   // this loop may be to M: for (m=0; m<M; m++) 
            {
                float32_t r_re, y_re;
                r_re = R[(k + (n + m*N)*L)];
                y_re = Y[(k + (p + m*P)*L)];
                C_re += (y_re*r_re);
            }
            A_re = B_re - C_re;
            Y[(k + (p + n*P)*L)] = A_re*D[k];
        }
    }
#endif

    int m, p;
          xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pA;
    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pD = (const xb_vecN_2xf32 *)_D;

    xb_vecN_2xf32 A, B, R, Y, D;
    xb_vecN_2xf32 Acc0, Acc1;

    NASSERT_ALIGN(_Y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_A,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(_D,(2*BBE_SIMD_WIDTH));     
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    NASSERT(M > 0 && P > 0 && L > 0);

    // calculate A(:,n)'*B-Rn'*Y, 1xP
    for (p = 0; p < P; p++)
    {
        Acc0 = Acc1 = BBE_ZERON_2XF32();
        pA = (const xb_vecN_2xf32 *)(_A + n*L);
        pB = (const xb_vecN_2xf32 *)(_B + p*L);
        pR = (const xb_vecN_2xf32 *)(_R + n*L);
        pY = (      xb_vecN_2xf32 *)(_Y + p*L);

        for (m = 0; m < M; m++)
        {
            BBE_LVN_2XF32_XP(A, pA, N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(B, pB, P*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc0, A, B);
        }
        for (m = 0; m < n; m++) // this loop may be to M: for (m=0; m<M; m++) 
        {
            BBE_LVN_2XF32_XP(R, pR, N*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(Y, pY, P*L * sizeof(float32_t));
            BBE_MULAN_2XF32(Acc1, R, Y);
        }
        Acc0 = BBE_SUBN_2XF32(Acc0, Acc1);
        BBE_LVN_2XF32_IP(D, pD, 0);
        Acc0 = BBE_MULN_2XF32(Acc0, D);
        BBE_SVN_2XF32_IP(Acc0, pY, 2 * BBE_SIMD_WIDTH);
    }
}

void rcholfwdmxnxpsf (
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
            int M, int N, int P, int L)
{
    int l,n;

    // check alignment
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT(M>1);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);

    for (l=0; l<L; l+=VECLEN)
    {
        for (n=0; n<N; n++)
        {
            iter2(y, R, A, B, D, n,M,N,P,L);
            D+=VECLEN;
        }
        A+=VECLEN;  
        R+=VECLEN;
        B+=VECLEN;
        y+=VECLEN;
    }
}
#endif

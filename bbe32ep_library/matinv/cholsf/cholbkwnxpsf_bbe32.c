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
  Cholesky backward recursion for pseudo-inversion API (complex floating point 
  data)
  C code optimized for BBE32

  Integrit 2006-2017
===========================================================================*/
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward 
recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Data format: IEEE-754 Std. single precision floating-point

Input:
N                Matrix dimension (number of columns and rows in 
                 matrices R)
P                Number of columns in right-side matrices B
L                Number of matrices
R[N*N][L]        Cholesky upper triangular matrices R
D[L*N]           Reciprocal of main diagonal written in a special format
y[N*P][L]        Results of forward recursion stage
Output:
x[N*P][L]        Decision matrix x

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex-valued 
   functions, or a multiple of 8 for real-valued functions.
3. N must be greater than 1, P must be >=1
---------------------------------------------------------------------------*/
#if !(HAVE_VFPU)
DISCARD_FUN(void, cholbkwnxpsf,(
                  complex_float* restrict x, 
            const complex_float* restrict   R, 
            const complex_float* restrict  D,
            const complex_float* restrict   y, 
            int N, int P, int L))
#else
/*
    Reference code:
    function [X]=cholbkw(R,y, extraBits, isFixedPoint, isShow)
    sz=size(y); N=sz(1);P=sz(2);
    X=zeros(N,P);
    D=real(1./diag(R));
    q=ceil(log2(min(abs(D))));
    q=15-q;
    q=q-extraBits;
    for m=N:-1:1
        Rm=R(m,:); 
        ym=y(m,:);
        x=(ym-Rm*X)*D(m);
        X(m,:)=x;
    end
*/

#define VECLEN (BBE_SIMD_WIDTH/4)

/*
    clean output
*/
//static void cleanX(complex_float *_X,int N,int P,int L)
//{
//    int k;
//    xb_vecN_2xf32 * restrict pX = (xb_vecN_2xf32 *)_X;
//
//    NASSERT_ALIGN(_X, (2 * BBE_SIMD_WIDTH));
//    NASSERT((L&(BBE_SIMD_WIDTH / 4 - 1)) == 0);
//    NASSERT(N > 0 && P > 0 && L > 0);
//
//    for (k = 0; k < (N*P*L / (BBE_SIMD_WIDTH / 4)); k++)
//    {
//        BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(), pX, 2 * BBE_SIMD_WIDTH);
//    }
//}

/*
    one recursion iteration:
    Input/output:
    X[8*2][N*P]  input/output decision
    Input:
    R[L][2*N*N]  upper triangle matrix
    y[L][2*N*P]  right side of equation
    D[L][N]      reciprocal of main diagonal
    m            iteration index (from M-1 downwards to 0)
    scratch:
    Bm[P*2*8]

*/
static void iter( 
                complex_float* restrict X, 
          const complex_float* restrict R, 
          const complex_float* restrict Y, 
          const complex_float* restrict D,
          int m,int N,int P,int L)
{
#if 0
    float32_t* restrict X = (float32_t*)_X;
    const float32_t* restrict D = (const float32_t*)_D;
    const float32_t* restrict R = (const float32_t*)_R;
    const float32_t* restrict y = (const float32_t*)_y;
    int k, n, p;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT((L&(BBE_SIMD_WIDTH / 4 - 1)) == 0);
    NASSERT(N>0 && P>0 && L>0);
    y += m * 2 * P*L;
    R += m * 2 * N*L;
    // calculate y(m,:)-R(m,:)*X, 1xP
    for (k = 0; k<VECLEN; k++)
    {
        for (p = 0; p<P; p++)
        {
            float32_t A_re = 0.f, A_im = 0.f;
            float32_t B_re = 0.f, B_im = 0.f;
            A_re = (y[(k + p*L) * 2 + 0]);
            A_im = (y[(k + p*L) * 2 + 1]);
            for (n = m + 1; n<N; n++)   // NOTE: loop may begin from 0: for (n=0; n<N; n++) 
            {
                float32_t r_re, r_im;
                float32_t x_re, x_im;
                r_re = R[(k + n*L) * 2 + 0];
                r_im = R[(k + n*L) * 2 + 1];
                x_re = X[(k + (p + n*P)*L) * 2 + 0];
                x_im = X[(k + (p + n*P)*L) * 2 + 1];
                B_re += (x_re * r_re) - (x_im * r_im);    // ->qX+qR
                B_im += (x_re * r_im) + (x_im * r_re);
            }
            A_re = A_re - B_re;
            A_im = A_im - B_im;
            X[(k + (p + m*P)*L) * 2 + 0] = A_re*D[2 * k + 0];
            X[(k + (p + m*P)*L) * 2 + 1] = A_im*D[2 * k + 1];
        }
    }
#endif

    int n, p;
          xb_vecN_2xf32 * restrict pXw = (xb_vecN_2xf32 *)(X + m*P*L);
    const xb_vecN_2xf32 * restrict pXr;
    const xb_vecN_2xf32 * restrict pY = (const xb_vecN_2xf32 *)(Y + m*P*L);
    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pD = (const xb_vecN_2xf32 *)D;

    xb_vecN_2xf32 Acc0, Acc1, Acc2, Acc3;
    xb_vecN_2xf32 x0, r0, d0;

    NASSERT_ALIGN(X,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT((L&(BBE_SIMD_WIDTH/4-1))==0);
    NASSERT(N > 1 && P > 0 && L > 0);

    // calculate y(m,:)-R(m,:)*X, 1xP
    for (p = 0; p < P; p++)
    {
        pR = (const xb_vecN_2xf32 *)(R + m*N*L + (m + 1)*L);
        pXr = (const xb_vecN_2xf32 *)(X + (p + (m + 1)*P)*L);
        BBE_LVN_2XF32_XP(Acc0, pY, 2 * L * sizeof(float32_t));
        Acc1 = Acc2 = Acc3 = BBE_ZERON_2XF32();

        for (n = 0; n < (N - m - 1) >> 1; n++)
        {
            BBE_LVN_2XF32_XP(x0, pXr, 2 * P*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(r0, pR, 2 * L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc0, x0, r0, 3, 4);
            BBE_MULMASN_2XF32(Acc1, x0, r0, 2, 11);
            BBE_LVN_2XF32_XP(x0, pXr, 2 * P*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(r0, pR, 2 * L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc2, x0, r0, 3, 4);
            BBE_MULMASN_2XF32(Acc3, x0, r0, 2, 11);
        }
        if ((N - m - 1) & 1)
        {
            BBE_LVN_2XF32_XP(x0, pXr, 2 * P*L * sizeof(float32_t));
            BBE_LVN_2XF32_XP(r0, pR, 2 * L * sizeof(float32_t));
            BBE_MULMASN_2XF32(Acc0, x0, r0, 3, 4);
            BBE_MULMASN_2XF32(Acc1, x0, r0, 2, 11);
        }
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);

        BBE_LVN_2XF32_IP(d0, pD, 0);
        Acc0 = BBE_MULN_2XF32(Acc0, d0);
        BBE_SVN_2XF32_XP(Acc0, pXw, 2 * L * sizeof(float32_t));
    }
}

void cholbkwnxpsf (
                  complex_float* restrict x, 
            const complex_float* restrict R, 
            const complex_float* restrict D,
            const complex_float* restrict y, 
            int N, int P, int L)
{
    int m,l;

    // check alignment
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(N>1);
    NASSERT(P>0);
    NASSERT((L&(BBE_SIMD_WIDTH/4-1))==0);
    
    D+=VECLEN*(N-1);
    //cleanX(x,N,P,L);
    for (l=0; l<L; l+=VECLEN)
    {
        for (m=N-1; m>=0; m--)
        {
            iter(x,R,y,D-(N-1-m)*VECLEN,m,N,P,L);
        }
        D+=VECLEN*N;
        // go to the next matrices
        R+=VECLEN;
        y+=VECLEN;
        x+=VECLEN;
    }
}
#endif

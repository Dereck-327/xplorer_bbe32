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
  Cholesky backward recursion for pseudo-inversion API (real floating point 
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
DISCARD_FUN(void, rcholbkw4x1sf,(
                  float32_t* restrict x, 
            const float32_t* restrict R, 
            const float32_t* restrict D,
            const float32_t* restrict y, 
            int L))
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

#define VECLEN (BBE_SIMD_WIDTH/2)

void rcholbkw4x1sf(
                  float32_t* restrict x, 
            const float32_t* restrict R, 
            const float32_t* restrict D,
            const float32_t* restrict y, 
            int L)
{
    int l;
          xb_vecN_2xf32 * restrict pXw;
    const xb_vecN_2xf32 * restrict pXr;
    const xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pD;

    xb_vecN_2xf32 Acc, Acc1;
    xb_vecN_2xf32 x0, r0, d0;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
    if (L<=0) return;

    pR = (const xb_vecN_2xf32 *)(R + 11 * L);
    pD = (const xb_vecN_2xf32 *)(D + VECLEN * 3);
    pY = (const xb_vecN_2xf32 *)(y + 2 * L);
    pXw = (xb_vecN_2xf32 *)(x + 2 * L);
    for (l = 0; l < L; l += VECLEN)
    {
        Acc = BBE_LVN_2XF32_X(pY, L * sizeof(float32_t));

        BBE_LVN_2XF32_IP(d0, pD, -2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, d0);
        BBE_SVN_2XF32_X(Acc, pXw, L * sizeof(float32_t));

        x0 = Acc;
        BBE_LVN_2XF32_IP(Acc, pY, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(r0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc, x0, r0);

        BBE_LVN_2XF32_IP(d0, pD, 10 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, d0);
        BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
    }

    pR = (const xb_vecN_2xf32 *)(R + 6 * L);
    pD = (const xb_vecN_2xf32 *)(D + VECLEN);
    pY = (const xb_vecN_2xf32 *)(y + L);
    pXr = (const xb_vecN_2xf32 *)(x + 2 * L);
    pXw = (xb_vecN_2xf32 *)(x + L);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(Acc, pY, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_XP(x0, pXr, L * sizeof(float32_t));
        BBE_LVN_2XF32_XP(r0, pR, L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc, x0, r0);
        BBE_LVN_2XF32_XP(x0, pXr, -L * sizeof(float32_t) + 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(r0, pR, -L * sizeof(float32_t) + 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc, x0, r0);

        BBE_LVN_2XF32_IP(d0, pD, 8 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, d0);
        BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
    }

    pR = (const xb_vecN_2xf32 *)(R + L);
    pD = (const xb_vecN_2xf32 *)(D);
    pY = (const xb_vecN_2xf32 *)(y);
    pXr = (const xb_vecN_2xf32 *)(x + L);
    pXw = (xb_vecN_2xf32 *)(x);
    for (l = 0; l < L; l += VECLEN)
    {
        BBE_LVN_2XF32_IP(Acc, pY, 2 * BBE_SIMD_WIDTH);

        x0 = BBE_LVN_2XF32_X(pXr, 2 * L * sizeof(float32_t));
        r0 = BBE_LVN_2XF32_X(pR, 2 * L * sizeof(float32_t));
        Acc1 = BBE_MULN_2XF32(x0, r0);

        x0 = BBE_LVN_2XF32_X(pXr, L * sizeof(float32_t));
        r0 = BBE_LVN_2XF32_X(pR, L * sizeof(float32_t));
        BBE_MULSN_2XF32(Acc, x0, r0);

        BBE_LVN_2XF32_IP(x0, pXr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(r0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, x0, r0);

        Acc = BBE_SUBN_2XF32(Acc, Acc1);


        BBE_LVN_2XF32_IP(d0, pD, 8 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, d0);
        BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

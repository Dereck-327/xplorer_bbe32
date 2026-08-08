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
  Cholesky forward recursion for pseudo-inversion API (complex floating point 
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
DISCARD_FUN(void, cholfwd4x4x1sf,(
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
            int L))
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

#define VECLEN (BBE_SIMD_WIDTH/4)

void cholfwd4x4x1sf(
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
            int L)
{
    int l;
          xb_vecN_2xf32 * restrict pY;
    const xb_vecN_2xf32 * restrict pYr;
    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pA;
    const xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pD;

    xb_vecN_2xf32 A0, B0, R0, Y0, D0;
    xb_vecN_2xf32 Acc0, Acc1, Acc2, Acc3;

    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B, (2 * BBE_SIMD_WIDTH));
    NASSERT((L&(BBE_SIMD_WIDTH / 4 - 1)) == 0);
    if (L<=0) return;

    pA = (const xb_vecN_2xf32 *)(A);
    pB = (const xb_vecN_2xf32 *)(B);
    pY = (xb_vecN_2xf32 *)(y);
    pD = (const xb_vecN_2xf32 *)D;
    for (l = 0; l < L; l += VECLEN)
    {
        A0 = BBE_LVN_2XF32_X(pA, 6 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 3 * 2 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 4 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * 2 * L * sizeof(float32_t));
        Acc2 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc3 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 2 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, B0, 2, 11);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);

        BBE_LVN_2XF32_IP(D0, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, D0);
        BBE_SVN_2XF32_IP(Acc0, pY, 2 * BBE_SIMD_WIDTH);
    }

    pA = (const xb_vecN_2xf32 *)(A + L);
    pB = (const xb_vecN_2xf32 *)(B);
    pR = (const xb_vecN_2xf32 *)(R + L);
    pYr = (const xb_vecN_2xf32 *)(y);
    pY = (xb_vecN_2xf32 *)(y + L);
    pD = (const xb_vecN_2xf32 *)D + 1;
    for (l = 0; l < L; l += VECLEN)
    {
        A0 = BBE_LVN_2XF32_X(pA, 6 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 3 * 2 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 4 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * 2 * L * sizeof(float32_t));
        Acc2 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc3 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 2 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, B0, 2, 11);

        BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Y0, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc0, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc1, R0, Y0, 1, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
        BBE_LVN_2XF32_IP(D0, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, D0);
        BBE_SVN_2XF32_IP(Acc0, pY, 2 * BBE_SIMD_WIDTH);
    }

    pA = (const xb_vecN_2xf32 *)(A + 2 * L);
    pB = (const xb_vecN_2xf32 *)(B);
    pR = (const xb_vecN_2xf32 *)(R + 2 * L);
    pYr = (const xb_vecN_2xf32 *)(y);
    pY = (xb_vecN_2xf32 *)(y + 2 * L);
    pD = (const xb_vecN_2xf32 *)D + 2;
    for (l = 0; l < L; l += VECLEN)
    {
        A0 = BBE_LVN_2XF32_X(pA, 6 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 3 * 2 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 4 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * 2 * L * sizeof(float32_t));
        Acc2 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc3 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 2 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, B0, 2, 11);

        R0 = BBE_LVN_2XF32_X(pR, 2 * 4 * L * sizeof(float32_t));
        Y0 = BBE_LVN_2XF32_X(pYr, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc1, R0, Y0, 1, 11);
        BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Y0, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc3, R0, Y0, 1, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
        BBE_LVN_2XF32_IP(D0, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, D0);
        BBE_SVN_2XF32_IP(Acc0, pY, 2 * BBE_SIMD_WIDTH);
    }

    pA = (const xb_vecN_2xf32 *)(A + 3 * L);
    pB = (const xb_vecN_2xf32 *)(B);
    pR = (const xb_vecN_2xf32 *)(R + 3 * L);
    pYr = (const xb_vecN_2xf32 *)(y);
    pY = (xb_vecN_2xf32 *)(y + 3 * L);
    pD = (const xb_vecN_2xf32 *)D + 3;
    for (l = 0; l < L; l += VECLEN)
    {
        A0 = BBE_LVN_2XF32_X(pA, 6 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 3 * 2 * L * sizeof(float32_t));
        Acc0 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc1 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 4 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * 2 * L * sizeof(float32_t));
        Acc2 = BBE_MULMN_2XF32(A0, B0, 0, 4);
        Acc3 = BBE_MULMN_2XF32(A0, B0, 2, 11);
        A0 = BBE_LVN_2XF32_X(pA, 2 * 4 * L * sizeof(float32_t));
        B0 = BBE_LVN_2XF32_X(pB, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc1, A0, B0, 2, 11);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc2, A0, B0, 0, 4);
        BBE_MULMASN_2XF32(Acc3, A0, B0, 2, 11);

        R0 = BBE_LVN_2XF32_X(pR, 4 * 4 * L * sizeof(float32_t));
        Y0 = BBE_LVN_2XF32_X(pYr, 2 * 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc0, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc1, R0, Y0, 1, 11);
        R0 = BBE_LVN_2XF32_X(pR, 2 * 4 * L * sizeof(float32_t));
        Y0 = BBE_LVN_2XF32_X(pYr, 2 * L * sizeof(float32_t));
        BBE_MULMASN_2XF32(Acc2, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc3, R0, Y0, 1, 11);
        BBE_LVN_2XF32_XP(R0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Y0, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_MULMASN_2XF32(Acc0, R0, Y0, 3, 4);
        BBE_MULMASN_2XF32(Acc1, R0, Y0, 1, 11);

        Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
        BBE_LVN_2XF32_IP(D0, pD, 8 * BBE_SIMD_WIDTH);
        Acc0 = BBE_MULN_2XF32(Acc0, D0);
        BBE_SVN_2XF32_IP(Acc0, pY, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

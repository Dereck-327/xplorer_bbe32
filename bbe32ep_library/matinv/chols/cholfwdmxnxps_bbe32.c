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
  NatureDSP_Baseband library. Cholesky forward recursion for pseudo-inversion API (complex data)
    These functions make forward recursion stage of pseudo-inversion. They use
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
M               Matrix dimension (number of rows in matrices A)
N               Matrix dimension (number of columns and rows in matrices 
                R)
P               Number of columns in right-side matrices B
L               Number of matrices
R[N*N][L]       Cholesky upper triangular matrices R
A[M*N][L]       Original left-side matrices A
B[M*P][L]       Original right-side matrices B
D[L/4][N][8]    Reciprocal of main diagonal (mantissa, exponent) in the 
                special format
qA,qB,qY        Fixed point representation of matrices A (or R which is 
                the same), B and y
Output:
y[N*P][L]       Decision matrix y

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8
3. Matrix sizes (M,N) must be greater than 1, P must be >=1
4. M >= N
---------------------------------------------------------------------------*/

void cholfwdmxnxps (
                  complex_fract16* restrict _y,
            const complex_fract16* restrict _R, 
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _A, 
            const complex_fract16* restrict _B, 
            int qA,int qB,int qY,
            int M,int N, int P,
            int L)
{
          int16_t* restrict y=(      int16_t*)_y;
    const int16_t* restrict R=(const int16_t*)_R; 
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict A=(const int16_t*)_A; 
    const int16_t* restrict B=(const int16_t*)_B; 
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
    vsaN qYB=BBE_MOVVSA32(qY-qB),D1,_16=BBE_MOVVSA32(16);
    xb_vecNx16 A0,B0,Y0,R0,D0,t,lo,hi;
    const xb_vecNx16* restrict pD;
    const xb_vecNx16* restrict pA;
    const xb_vecNx16* restrict pB;
    const xb_vecNx16* restrict pR;
          xb_vecNx16* restrict pY;
    xb_vecNx40 wA; 

    int l,n;
    int m,p;

    // check alignment
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT(M>1);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);

    pD=(const xb_vecNx16*)D;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        for (n=0; n<N; n++)
        {
            BBE_LVNX16_IP(D0,pD,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t ,pD,2*BBE_SIMD_WIDTH);
            D1=BBE_MOVVSV(t,0);
            #ifdef COMPILER_XTENSA
            #pragma concurrent
            #endif
            for(p=0; p<P; p++)
            {
                wA=0;
                pA=(const xb_vecNx16*)A;
                pB=(const xb_vecNx16*)B;
                #ifdef COMPILER_XTENSA
                #pragma loop_count min=2
                #endif
                for (m=0; m<M; m++) 
                {
                    BBE_LVNX16_XP(A0,pA,4*N*L);
                    BBE_LVNX16_XP(B0,pB,4*P*L);
                    BBE_MULANX16J(wA,B0,A0);
                }
                wA=BBE_SLSNX40(wA,qYB);
                pR=(const xb_vecNx16*)R;
                pY=(      xb_vecNx16*)y;
                for (m=0; m<n; m++) 
                {
                    BBE_LVNX16_XP(R0,pR,4*N*L);
                    BBE_LVNX16_XP(Y0,pY,4*P*L);
                    BBE_MULSNX16J(wA,Y0,R0);
                }
                lo=BBE_PACKLNX40(wA);
                hi=BBE_PACKVNX40(wA,_16);
                wA=BBE_MULUUNX16(D0,lo);
                wA=BBE_SRAINX40(wA,16);
                BBE_MULUSANX16(wA,D0,hi);
                wA=BBE_SLSNX40(wA,D1);
                Y0=BBE_PACKPNX40(wA);
                BBE_SVNX16_XP(Y0,pY,4*P*L);
                B=(const int16_t*)XT_ADDX4(L,(uintptr_t)B);
                y=(      int16_t*)XT_ADDX4(L,(uintptr_t)y);
            }
            B=(const int16_t*)XT_ADDX4(-P*L,(uintptr_t)B);
            y=(      int16_t*)XT_ADDX4(-P*L,(uintptr_t)y);
            A=(const int16_t*)XT_ADDX4(L   ,(uintptr_t)A);
            R=(const int16_t*)XT_ADDX4(L   ,(uintptr_t)R);
        }
        A=(const int16_t*)XT_ADDX4(BBE_SIMD_WIDTH/2-N*L,(uintptr_t)A);
        R=(const int16_t*)XT_ADDX4(BBE_SIMD_WIDTH/2-N*L,(uintptr_t)R);
        // go to the next matrices
        B+=BBE_SIMD_WIDTH;
        y+=BBE_SIMD_WIDTH;
    }
} /* cholfwdmxnxps() */

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
  NatureDSP_Baseband library. Cholesky backward recursion for pseudo-inversion API (complex data)
    These functions make backward recursion stage of pseudo-inversion. They
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if !(HAVE_VSAMATH && 1)
DISCARD_FUN( void, cholbkwnxps,(
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int N, int P, int L))
#else

/*-------------------------------------------------------------------------
   These functions make backward recursion stage of pseudo-inversion. They
   use Cholesky decomposition of original matrices and results of forward 
   recursion. 
   NOTE:
   Data layout for matrices is selected as for other matrices written in a 
   streaming order. 

   Input:
   N                Matrix dimension (number of columns and rows in 
                    matrices R)
   P                Number of columns in right-side matrices B
   L                Number of matrices
   R[N*N][L]        Cholesky upper triangular matrices R
   D[L/4][N][8]     reciprocal of main diagonal (mantissa, exponent) in the 
                    special format
   y[N*P][L]        Results of forward recursion stage
   qA,qX,qY         fixed point representation of matrices A(or R which is 
                    the same), x and y
   Output:
   x[N*P][L]        Decision matrix x


   Restrictions:
   1. All matrices must not overlap and be aligned on 32-byte boundary 
   2. Number of matrices L must be a multiple of 8
   3. Matrix sizes (M,N) must be greater than 1, P must be >=1
   4. qX+qA-qY must be <=16 
---------------------------------------------------------------------------*/

void cholbkwnxps (
                  complex_fract16* restrict _x, 
            const complex_fract16* restrict _R,
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _y, 
            int qA, int qY, int qX,
            int N, int P, int L)
{
          int16_t* restrict x=(      int16_t*)_x;
    const int16_t* restrict R=(const int16_t*)_R;
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict y=(const int16_t*)_y;
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
    vsaN qXAY=BBE_MOVVSA32(qX+qA-qY),D1;
    const xb_vecNx16* restrict pD;
    const xb_vecNx16* restrict pR;
    const xb_vecNx16* restrict pY;
          xb_vecNx16* restrict pX;
    xb_vecNx40 wA;
    xb_vecNx16 X0,Y0,R0,D0,t;
    int m,l;

    // check alignment
    NASSERT_ALIGN(x,BBE_SIMD_WIDTH*2);
    NASSERT_ALIGN(R,BBE_SIMD_WIDTH*2);
    NASSERT_ALIGN(y,BBE_SIMD_WIDTH*2);
    NASSERT_ALIGN(D,BBE_SIMD_WIDTH*2);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    
    X0=0;
    pX=(xb_vecNx16*)x;
    for (l=0; l<(N*P*L)>>(LOG2_BBE_SIMD_WIDTH-1); l++)
    {
        BBE_SVNX16_IP(X0,pX,2*BBE_SIMD_WIDTH);
    }

    D+=(2*BBE_SIMD_WIDTH)*(N-1)+BBE_SIMD_WIDTH;
    R+=(N-1)*L*2;
    x+=(N*P-1)*L*2;
    y+=(N*P-1)*L*2;
    pD=(const xb_vecNx16*)(D);
    for (l=0; l<L>>(LOG2_BBE_SIMD_WIDTH-1); l++)
    {
        pY=(const xb_vecNx16*)(y);
        R+=N*2*N*L;
        for (m=0; m<N; m++)
        {
          int n,p;
          R=(const int16_t*)XT_ADDX4(-N*L,(uintptr_t)R);
          BBE_LVNX16_IP(t ,pD,-2*BBE_SIMD_WIDTH);
          BBE_LVNX16_IP(D0,pD,-2*BBE_SIMD_WIDTH);
          D1=BBE_MOVVSV(t,0);
          D1=BBE_SUBSAVSN(11,D1);
          pX=(      xb_vecNx16*)x;
          #ifdef COMPILER_XTENSA
          #pragma concurrent
          #endif
          for(p=P-1; p>=0; p--)
          {
              BBE_LVNX16_XP(Y0,pY,-4*L);
              wA=BBE_UNPKSNX16(Y0);
              wA=BBE_SLSNX40(wA,qXAY);
              pR=(const xb_vecNx16*)R;
              for (n=0;n<m;n++)
              {
                  BBE_LVNX16_XP(R0,pR,-4*L);
                  BBE_LVNX16_XP(X0,pX,-4*P*L);
                  BBE_MULSNX16C(wA,R0,X0);
              }
              wA=BBE_RNDSADJNX40(wA,D1);
              X0=BBE_PACKVNX40(wA,D1);
              wA=BBE_MULUSNX16(D0,X0);
              X0=BBE_PACKQNX40(wA);
              BBE_SVNX16_XP(X0,pX,4*(m*P-1)*L);
          }
        }
        // go to the next matrices
        pD+=4*N;
        R+=BBE_SIMD_WIDTH;
        y+=BBE_SIMD_WIDTH;
        x+=BBE_SIMD_WIDTH;
    }
} /* cholbkwnxps() */
#endif

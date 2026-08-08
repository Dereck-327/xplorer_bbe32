/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*
  NatureDSP_Baseband library. Banded Cholesky forward recursion for pseudo-inversion API (complex data)
    These functions make forward recursion stage of pseudo-inversion. They use
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "bcholn_common.h"

#if HAVE_BCHOLN
/*
  Reference Matlab code
%--------------------------------------------------------------------------
% forward recursion (substitution) for banded Cholesky decomposition where
% matrices R and A written in the compact form:
% y=R'\(A'*B), 
% R[WxN], A[WxN], B[MxP], y[NxP], M=N+W-1
%--------------------------------------------------------------------------
    function y=bcholfwd(R,A,B)
    sz=size(A); W=sz(1); N=sz(2); M=W+N-1; 
    sz=size(B); P=sz(2);
    y=zeros(N,P);
    D=1./R(W,:);    % diagonal elements

    x=A(:,1)'*B(1:W,:);
    y(1,:)=x*D(1);
    for n=2:W-1
        x=A(:,n)'*B(n:n+W-1,:);
        y(n,:)=(x-R(W-n+1:W-1,n)'*y(1:n-1,:))*D(n); % lengths from 1 to W-1
    end
    for n=W:N
        x=A(:,n)'*B(n:n+W-1,:);
        y(n,:)=(x-R(1:W-1,n)'*y(n-W+1:n-1,:))*D(n);     % length W
    end
*/

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the  next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion for 
specific band width W. They use Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD,SB,SY are selected as 
usual for complex block ordered matrix sequencies, i.e. total size is 
rounded up to the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)
SB=size((W+N-1)*P)
SY=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR][2]  Cholesky upper-triangle matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagonal elements are in the last raw)
D[L][SD][2]   Sequence of L reciprocals of main diagonal A represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
At[L][SA][2]  Original left-side matrices A represented in the compact 
              form (only band)
Bt[L][SB][2]  Original right-side matrices B. SB is computed as for complex 
              block ordered matrices of size (W+N-1)*P
qA,qB,qY      Fixed point representation of matrices A (or R which is the 
              same),B and y

Output:
Yt[L][SY][2]  Decision matrix y. SY is computed as for complex 
              block ordered matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
// W==4
int bcholfwd4xnx1n (
            int16_t* restrict Yt,
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict At, 
            const int16_t* restrict Bt, 
            int qYB,
            int N, int L)
{
    const int W=4;
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAt;
    const xb_vecNx16* restrict pB;
    const xb_vecNx16* restrict pY;

    const int16_t* restrict _D;
          int16_t* restrict _Yt;
    const int16_t* restrict _Bt;


    xb_vecNx16 t,AN,RN,BP,D0,YP,lo,hi;
    xb_vecNx40 W0;
    xb_c40 I0;
    valign va,vr,vb,vy;
    vsaN D1,vsa_qYB=BBE_MOVVSA32(qYB);
    int l,n,M;
    int SB=2*getSpace((W+N-1));
    int SY=2*getSpace(N);
    int SA=2*N*W;

    // check alignment
    NASSERT_ALIGN(Yt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Bt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);

    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(L>0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);

    {
        /* Zero the output array */
        xb_vecNx16 * py_st;
        py_st = (xb_vecNx16*)(Yt);
        YP = BBE_ZERONX16();
        for (n = 0; n < ((L*SY)>>LOG2_BBE_SIMD_WIDTH); n++)
        {
            BBE_SVNX16_IP(YP, py_st, sizeof(int16_t)*BBE_SIMD_WIDTH);
        }
        if ((L*SY)&(BBE_SIMD_WIDTH/2))
        {
            vy = BBE_ZALIGN();
            BBE_SAVNX16_XP(YP, vy, py_st, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
            BBE_SANX16POS_FP(vy, py_st);
        }
    }

    pRt=(const xb_vecNx16*)Rt;
    pAt=(const xb_vecNx16*)At;
    _D=D;
    _Yt=Yt;
    _Bt=Bt;
    for(n=0; n<N; n++)
    {
        M = XT_MIN((W-1),n);
        Yt=_Yt;
        D =_D;
        Bt=_Bt;
        pB =(const xb_vecNx16*)_Bt;
        pAt=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        pY =(const xb_vecNx16*)(Yt);
        pY =(      xb_vecNx16*)XT_ADDX4(-M,(uintptr_t)pY);
        __Pragma("concurrent");
        __Pragma("loop_count min=1");
        for (l=0; l<L; l++)
        {
            va=BBE_LA_PP(pAt);
            pRt=(const xb_vecNx16*)XT_ADDX4((W-1-M),(uintptr_t)pRt);
            vr=BBE_LA_PP(pRt);
            BBE_LAVNX16_XP(RN,vr,pRt,M*4);
            pRt=(const xb_vecNx16*)XT_ADDI_N((uintptr_t)pRt,4);
            BBE_LAVNX16_XP(AN,va,pAt,W*4);
            pAt=(const xb_vecNx16*)XT_ADDX4(SA/2-W,(uintptr_t)pAt);
            pRt=(const xb_vecNx16*)XT_ADDX4(SA/2-W,(uintptr_t)pRt);
            // load diagonals
            BBE_LPNX16_XP(D0,D,2*SY);
            t=BBE_REPNX16(D0,1);
            D1=BBE_MOVVSV(t,0);
            D0=BBE_REPNX16(D0,0);

            // calculate A(:,n)'*B-Rn'*Y, 1xP
            vb=BBE_LA_PP(pB);
            BBE_LAVNX16_XP(BP,vb,pB,W*4);
            pB=(const xb_vecNx16*)XT_ADDX4(SB/2-W,(uintptr_t)pB);
            W0=BBE_MULNX16J(BP,AN);
            //R(W-n+1:W-1,n)  y(1    :n-1,:)
            //R(    1:W-1,n)  y(n-W+1:n-1,:)
            W0=BBE_SLSNX40(W0,vsa_qYB); // representation qA+qB->qA+qY-16
            vy=BBE_LA_PP(pY);
            BBE_LAVNX16_XP(YP,vy,pY,M*4);
            pY=(const xb_vecNx16*)XT_ADDX4(SY/2-M,(uintptr_t)pY);
            BBE_MULSNX16J(W0,YP,RN);
            I0=BBE_RADDNX40C(W0);
            W0=BBE_MOVNX40_FROMC40(I0);
            lo=BBE_PACKLNX40(W0);
            W0=BBE_SRAINX40(W0,16);
            hi=BBE_PACKSNX40(W0);
            W0=BBE_MULUUNX16(lo,D0);
            W0=BBE_SRAINX40(W0,16);
            BBE_MULUSANX16(W0,D0,hi);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            BBE_SPNX16_XP(t,Yt,2*SY);
        }
        _Yt+=2;
        _D +=2;
        _Bt+=2;
        At=(const int16_t*)XT_ADDX4(W,(uintptr_t)At);
        Rt=(const int16_t*)XT_ADDX4(W,(uintptr_t)Rt);
    }

    return 0;
}
#endif

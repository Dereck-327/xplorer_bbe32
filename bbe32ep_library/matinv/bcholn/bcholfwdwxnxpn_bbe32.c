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

#if !HAVE_BCHOLN
DISCARD_FUN(void,bcholfwdwxnxpn ,(
                  complex_fract16* restrict Yt,
            const complex_fract16* restrict Rt, 
            const complex_fract16* restrict D, 
            const complex_fract16* restrict At, 
            const complex_fract16* restrict Bt, 
            int qB,int qY,
            int W, int N, int P, int L))
#else
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
inline_ int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the  next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}


/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
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
Rt[L][SR]     Cholesky upper triangular matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagoanal elements are in the last raw)
D[L][SD]      Sequence of L reciprocals of main diagonal A represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
At[L][SA]     Original left-side matrices A represented in the compact 
              form (only band)
Bt[L][SB]     Original right-side matrices B. SB is computed as for complex 
              block ordered matrices of size (W+N-1)*P
qA,qB,qY      Fixed point representation of matrices A (or R which is the 
              same),B and y

Output:
Yt[L][SY]     Decision matrix y. SY is computed as for complex 
              block ordered matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
void bcholfwdwxnxpn (
                  complex_fract16* restrict _Yt,
            const complex_fract16* restrict _Rt, 
            const complex_fract16* restrict _D, 
            const complex_fract16* restrict _At, 
            const complex_fract16* restrict _Bt, 
            int qB,int qY,
            int W, int N, int P, int L)
{
          int16_t* restrict Yt=(      int16_t*)_Yt;
    const int16_t* restrict Rt=(const int16_t*)_Rt; 
    const int16_t* restrict D =(const int16_t*)_D ; 
    const int16_t* restrict At=(const int16_t*)_At; 
    const int16_t* restrict Bt=(const int16_t*)_Bt; 

    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAt;
    const xb_vecNx16* restrict pB;
    const xb_vecNx16* restrict pY;

    xb_vecNx16 t,AN0,AN1,RN0,RN1,BP0,BP1,D0,YP0,YP1,lo,hi;
    xb_vecNx40 W0;
    xb_c40 I0;
    valign va,vr,vb,vy;
    vsaN D1,qYB=BBE_MOVVSA32(qY-qB);
    int l,n,p,M;
    int SB,SD,SY;

    // check alignment
    NASSERT_ALIGN(Yt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Bt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);

    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);
    NASSERT(P>=1);
    
    if (L<=0) return;
    if (P==1)
    {
        typedef int (*fnfwd)(int16_t* ,const int16_t* , const int16_t* , const int16_t* , const int16_t* , int ,int, int );
        static const fnfwd fwd[]={bcholfwd4xnx1n,bcholfwd8xnx1n,bcholfwd12xnx1n,bcholfwd16xnx1n};
        fwd[(W>>2)-1](Yt,Rt, D, At, Bt, qY-qB,N,L);
        return;
    }

    SB = 2*getSpace((W+N-1)*P);
    SD = 2*getSpace(N);
    SY = 2*getSpace(N*P);
    

    {
        /* Zero the output array */
        xb_vecNx16 * py_st;
        py_st = (xb_vecNx16*)(Yt);
        YP0 = BBE_ZERONX16();
        for (n = 0; n < ((L*SY)>>LOG2_BBE_SIMD_WIDTH); n++)
        {
            BBE_SVNX16_IP(YP0, py_st, sizeof(int16_t)*BBE_SIMD_WIDTH);
        }
        if ((L*SY)&(BBE_SIMD_WIDTH/2))
        {
            vy = BBE_ZALIGN();
            BBE_SAVNX16_XP(YP0, vy, py_st, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
            BBE_SANX16POS_FP(vy, py_st);
        }
    }

    pRt=(const xb_vecNx16*)Rt;
    pAt=(const xb_vecNx16*)At;
    for (l=0; l<L; l++)
    {
        //iter(Yt, Rt,D, At, Bt, qB,qY, W,N,P);
        NASSERT_ALIGN(pAt,2*BBE_SIMD_WIDTH);
        NASSERT_ALIGN(pRt,2*BBE_SIMD_WIDTH);
        // these will be possible aligned on 16 bytes !
        NASSERT_ALIGN(Yt,16);
        NASSERT_ALIGN(Bt,16);
        NASSERT_ALIGN(D,16); 
        NASSERT(W>1 && P>0 );
        for(n=0; n<N; n++)
        {
            M = XT_MIN(W-1,n);
            pRt=(const xb_vecNx16*)XT_ADDX4((W-1-M),(uintptr_t)pRt);
            vr=BBE_LA_PP(pRt);
            BBE_LAVNX16_XP(RN0,vr,pRt,M*4);
            BBE_LAVNX16_XP(RN1,vr,pRt,M*4-2*BBE_SIMD_WIDTH);
            pRt=(const xb_vecNx16*)XT_ADDI_N((uintptr_t)pRt,4);
            va=BBE_LA_PP(pAt);
            BBE_LAVNX16_XP(AN0,va,pAt,W*4);
            BBE_LAVNX16_XP(AN1,va,pAt,W*4-2*BBE_SIMD_WIDTH);
            // load diagonals
            BBE_LPNX16_IP(D0,D,4);
            t=BBE_REPNX16(D0,1);
            D1=BBE_MOVVSV(t,0);
            D0=BBE_REPNX16(D0,0);

            // calculate A(:,n)'*B-Rn'*Y, 1xP
            pY=(xb_vecNx16*)Yt;
            pB=(xb_vecNx16*)Bt;
            pY=(      xb_vecNx16*)XT_ADDX4(-N,(uintptr_t)pY);
            __Pragma("concurrent");
            __Pragma("loop_count min=1");
            for(p=0; p<P; p++)
            {
                vb=BBE_LA_PP(pB);
                BBE_LAVNX16_XP(BP0,vb,pB,W*4);
                BBE_LAVNX16_XP(BP1,vb,pB,W*4-2*BBE_SIMD_WIDTH);
                pB=(const xb_vecNx16*)XT_ADDX4((N-1),(uintptr_t)pB);
                pY=(      xb_vecNx16*)XT_ADDX4((N-M),(uintptr_t)pY);
                vy=BBE_LA_PP(pY);
                BBE_LAVNX16_XP(YP0,vy,pY,M*4);
                BBE_LAVNX16_XP(YP1,vy,pY,M*4-2*BBE_SIMD_WIDTH);
                W0=BBE_MULNX16J(BP0,AN0);
                BBE_MULANX16J(W0,BP1,AN1);
                //R(W-n+1:W-1,n)  y(1    :n-1,:)
                //R(    1:W-1,n)  y(n-W+1:n-1,:)
                W0=BBE_SLSNX40(W0,qYB); // representation qA+qB->qA+qY-16
                BBE_MULSNX16J(W0,YP0,RN0);
                BBE_MULSNX16J(W0,YP1,RN1);
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
                BBE_SPNX16_XP(t,Yt,N*4);
            }
            Yt=(      int16_t*)XT_ADDX4((-N*P+1),(uintptr_t)Yt);
            Bt=(const int16_t*)XT_ADDI_N((uintptr_t)Bt,4);
        }
        // go to the next matrices
        Bt+=SB-2*N;
        Yt+=SY-2*N;
        D +=SD-2*N;
    }
} /* bcholfwdwxnxpn() */
#endif

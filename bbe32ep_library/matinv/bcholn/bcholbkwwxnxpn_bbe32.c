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
  NatureDSP_Baseband library. Banded Cholesky backward recursion for pseudo-inversion API (complex data)
    These functions make backward recursion stage of pseudo-inversion. They use
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
DISCARD_FUN(void,bcholbkwwxnxpn,(
                  complex_fract16* restrict Xt, 
            const complex_fract16* restrict Rt, 
            const complex_fract16* restrict D, 
            const complex_fract16* restrict Yt, 
            int qA, int qY, int qX,
            int W,int N, int P, int L))
#else

/*
    Reference code:
    %--------------------------------------------------------------------------
    % backward recursion (substitution) for banded Cholesky decomposition where
    % matrix R written in the compact form:
    % x=R\y, 
    % R[WxN], y[NxP], x[NxP], M=N+W-1
    %--------------------------------------------------------------------------
    function x=bcholbkw(R,y)
    sz=size(R); W=sz(1); N=sz(2); M=W+N-1; 
    sz=size(y); P=sz(2);
    x=zeros(N,P);
    D=1./R(W,:);    % diagonal elements

    n=N;
    x(n,:)=y(n,:)*D(n);     
    % lengths from 1 to W-1
    for m=1:W-1
    n=N-m;
    S=zeros(1,m);
    for k=1:m
        S(1,k)=R(W-k,N-m+k);
    end
    x(n,:)=(y(n,:)-S*x(n+1:N,:))*D(n);     
    end
    % fixed length W-1
    for m=W:N-1
    n=N-m;
    S=zeros(1,W-1);
    for k=1:W-1
        S(1,k)=R(W-k,N-m+k);
    end
    x(n,:)=(y(n,:)-S*x(n+1:n+W-1,:))*D(n);     
    end

*/

/*
    one recursion iteration:
    Input/output:
    X[L][2*N*P]  input/output decision
    Input:
    R[L][2*W*N]  upper triangle matrix
    y[L][2*N*P]  right side of equation
    D[L][N]      reciprocal of main diagonal
    m            iteration index (from M-1 downwards to 0)
    scratch:
    Bm[P*2*8]

    returns nonzero if overflow is detected
*/
static void iter1( 
          int16_t* restrict X, 
          const int16_t* restrict R, 
          const int16_t* restrict D, 
          const int16_t* restrict y, 
          int qXRY,
          int W,int N,int P)
{
    valign vx;
    xb_vecNx40 W0;
    xb_c40 I0;
    xb_vecNx16 D0,t,YP,RN0,RN1,XN0,XN1,lo,hi,cqXRY;
    vsaN D1,vqXRY=BBE_MOVVSA32(qXRY);
    const xb_vecNx16 * restrict pX;
    const int16_t* restrict pYp;
    int n,m,p,M;
    NASSERT_ALIGN(R,32);
    // these will be possible aligned on 16 bytes !
    NASSERT_ALIGN(X,16);
    NASSERT_ALIGN(y,16);
    NASSERT_ALIGN(D,16);
    NASSERT(N>0 && P>0);
    t=BBE_MOVVA16(1);
    cqXRY=BBE_SLLNX16(t,vqXRY);

    y+=N*2;
    X+=N*2;
    R+=(N*W-1)*2;
    D+=N*2-2;
    for (n=0; n<N; n++)
    {
        int idx;
        pX=(xb_vecNx16*)X;
        X-=2;
        y-=2;
        M=XT_MIN(n,(W-1));
        idx=0;
        RN0=RN1=0;
        for (m=0; m<M; m++)   
        {
            t=BBE_LPNX16_X(R,idx);
            RN1=BBE_SELNX16I(RN1,RN0,BBE_SELI_ROTATE_LEFT_2);
            RN0=BBE_SELNX16I(RN0,t,BBE_SELI_PACK_2);
            idx=XT_ADDX4(1-W,idx);
        }
        idx=-1;
        XT_MOVGEZ(idx,-W,n+(1-W));
        R=(int16_t*)XT_ADDX4(idx,(uintptr_t)R);
        // load diagonals
        BBE_LPNX16_IP(D0,D,-4);
        t=BBE_REPNX16(D0,1);
        D1=BBE_MOVVSV(t,0);
        D0=BBE_REPNX16(D0,0);

        // calculate y(m,:)-R(m,:)*X, 1xP
        pYp=y;
        #ifdef COMPILER_XTENSA
        #pragma concurrent
        #pragma loop_count min=1
        #endif
        for(p=0; p<P; p++)
        {
            YP=BBE_LPNX16_I(pYp,0);
            //W0=BBE_SLSNX40(YP,vqXRY);
            W0=BBE_MULUSNX16(cqXRY,YP);
            vx=BBE_LA_PP(pX);
            BBE_LAVNX16_XP(XN0,vx,pX,M*4);
            BBE_LAVNX16_XP(XN1,vx,pX,M*4-2*BBE_SIMD_WIDTH);
            BBE_MULSNX16C(W0,RN0,XN0);
            BBE_MULSNX16C(W0,RN1,XN1);
            I0=BBE_RADDNX40C(W0);
            W0=BBE_MOVNX40_FROMC40(I0);
            // NOTE: having this 32x16 multiple is critical !
            lo=BBE_PACKLNX40(W0);
            W0=BBE_SRAINX40(W0,16);
            hi=BBE_PACKSNX40(W0);
            W0=BBE_MULUUNX16(lo,D0);
            W0=BBE_SRAINX40(W0,16);
            BBE_MULUSANX16(W0,D0,hi);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            BBE_SPNX16_XP(t,X,N*4);
            //BBE_SPNX16_X(t,pX,-(M+1)*4);
            pYp=(const int16_t*)XT_ADDX4(N,(uintptr_t)pYp);
            pX =(   xb_vecNx16*)XT_ADDX4((N-M),(uintptr_t)pX);
        }
        X =(int16_t*)XT_ADDX4(-N*P,(uintptr_t)X);
    }
}

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
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SR,SD,SY,SX are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SR=size(W*N)
SD=size(N)
SY=size(N*P)
SX=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR]     Cholesky upper triangular matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagoanal elements are in the last raw)
D[L][SD]      Sequence of L reciprocals of main diagonal R represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
Yt[L][SY]     Results of forward recursion stage. SY is computed as for complex 
              block ordered matrices of size N*P
qA,qX,qY      Fixed point representation of matrices A(or R which is the same), 
              x and y
Output:
Xt[L][SX]     Decision matrix x. SX is computed as for complex block ordered 
              matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
void bcholbkwwxnxpn (
                    complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, 
                    const complex_fract16* restrict Yt, int qA, int qY, int qX,
                    int W,int N, int P, int L)
{
    int l;
    int qXYA;
    int SY,SX,SD;

    qXYA=qX-qY+qA;
    // check alignment
    NASSERT_ALIGN(Xt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Yt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D ,(2*BBE_SIMD_WIDTH));

    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);
    NASSERT(P>=1);
    if (L<=0) return;
    if (P == 1)
    {
        if (W == 4) { bcholbkw4xnx1n((int16_t*)Xt, (const int16_t*)Rt, (const int16_t*)D, (const int16_t*)Yt, qXYA, N, L);   return; }
        if (W == 8) { bcholbkw8xnx1n((int16_t*)Xt, (const int16_t*)Rt, (const int16_t*)D, (const int16_t*)Yt, qXYA, N, L);   return; }
        //if (W == 12)return bcholbkw12xnx1n(Xt, Rt, D, Yt, qXYA, N, L);
        if (W == 16){ bcholbkw16xnx1n((int16_t*)Xt, (const int16_t*)Rt, (const int16_t*)D, (const int16_t*)Yt, qXYA, N, L); return; }
    }
    
    SY=getSpace(N*P);
    SX=getSpace(N*P);
    SD=getSpace(N);
    
    {
        /* Zero the output array */
        xb_vecNx16 XN;
        xb_vecNx16 * px_st;
        int n;

        px_st = (xb_vecNx16*)(Xt);
        XN = BBE_ZERONX16();
        for (n = 0; n < ((L*SX*2)>>LOG2_BBE_SIMD_WIDTH); n++)
        {
            BBE_SVNX16_IP(XN, px_st, sizeof(int16_t)*BBE_SIMD_WIDTH);
        }
        if ((L*SX*2)&(BBE_SIMD_WIDTH/2))
        {
            valign vx;
            vx = BBE_ZALIGN();
            BBE_SAVNX16_XP(XN, vx, px_st, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
            BBE_SANX16POS_FP(vx, px_st);
        }
    }

    for (l=0; l<L; l++)
    {
        iter1((int16_t*)Xt,(const int16_t*)Rt,(const int16_t*)D,(const int16_t*)Yt,qXYA,W,N,P);
        // go to the next matrices
        Rt+=W*N;
        Yt+=SY;
        Xt+=SX;
        D +=SD;
    }
} /* bcholbkwwxnxpn() */
#endif

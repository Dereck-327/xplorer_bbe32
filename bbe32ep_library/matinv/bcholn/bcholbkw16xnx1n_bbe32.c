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

#if HAVE_BCHOLN

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
These functions make backward recursion stage of pseudo-inversion for 
specific band width W. They use Cholesky decomposition
of original matrices and results of forward recursion. 
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
Rt[L][SR][2]  Cholesky upper-triangle matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagonal elements are in the last raw)
D[L][SD][2]   Sequence of L reciprocals of main diagonal R represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
Yt[L][SY][2]  Results of forward recursion stage. SY is computed as for complex 
              block ordered matrices of size N*P
qA,qX,qY      Fixed point representation of matrices A(or R which is the same), 
              x and y
Output:
Xt[L][SX][2]  Decision matrix x. SX is computed as for complex block ordered 
              matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
// W==16
int bcholbkw16xnx1n  (
            int16_t* restrict Xt, 
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict Yt, 
            int qXYA, int N, int L)
{
    static const int16_t ALIGN(32) combine[32]={ 14,15,28,29,12,13,26,27,10,11,24,25,8,9,22,23,6,7,20,21,4,5,18,19,2,3,16,17,0,1,   0,0};
          int16_t*     restrict x;
    const int16_t*     restrict R;
    const int16_t*     restrict pD; 
    const int16_t*     restrict y; 
    const xb_vecNx16 * restrict pX;
    const xb_vecNx16 * restrict _pX;
    const int W=16;

    valign vx;
    xb_vecNx40 W0;
    xb_c40 I0;
    xb_vecNx16 D0,t,YP,RN0,RN1,XN0,XN1,lo,hi,cqXRY;
    vsaN D1,vqXRY=BBE_MOVVSA32(qXYA);
    int n,M,l,SN=getSpace(N);
    uint32_t idx,mod;
    vselN vcombine0,vcombine1;    
    // check alignment
    NASSERT_ALIGN(Xt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Yt,2*BBE_SIMD_WIDTH);
    
    NASSERT(L>0);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);
    
    {
        /* Zero the output array */
        xb_vecNx16 * px_st;
        px_st = (xb_vecNx16*)(Xt);
        XN0 = BBE_ZERONX16();
        for (n = 0; n < ((L*2*SN)>>LOG2_BBE_SIMD_WIDTH); n++)
        {
            BBE_SVNX16_IP(XN0, px_st, sizeof(int16_t)*BBE_SIMD_WIDTH);
        }
        if ((L*2*SN)&(BBE_SIMD_WIDTH/2))
        {
            vx = BBE_ZALIGN();
            BBE_SAVNX16_XP(XN0, vx, px_st, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
            BBE_SANX16POS_FP(vx, px_st);
        }
    }

    // prepare scale coefficient
    t=BBE_MOVVA16(1);
    cqXRY=BBE_SLLNX16(t,vqXRY);
    // form shuffle register for reversing W-1 elements
    t=BBE_LVNX16_I((const xb_vecNx16*)combine,0);
    vcombine0=BBE_MOVVSELNX16(t,0);
    t=BBE_LVNX16_I((const xb_vecNx16*)combine,2*BBE_SIMD_WIDTH);
    vcombine1=BBE_MOVVSELNX16(t,0);

    mod=((N*W)<<18)|((W-1)<<2);
    RN0=RN1=0;
    _pX=(const xb_vecNx16*)(Xt);
    idx=(N*W-1)<<2;
    Yt+=(N*2-2);
    Xt+=(N*2-2);
    D +=(N*2-2);
    _pX =(const xb_vecNx16*)XT_ADDI_N((uintptr_t)Xt,4);
    for (n=0; n<N; n++)
    {
        R=Rt;
        y=Yt;
        x=Xt;
        pD=D;
        pX =_pX;
        M=XT_MIN(n,(W-1));
        #ifdef COMPILER_XTENSA
        #pragma concurrent
        #pragma loop_count min=1
        #endif
        for (l=0; l<L; l++)
        {
            int k;
            xb_vecNx16 AN,RN;
            // load column
            k=idx;
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); AN=BBE_SELNX16I(AN,t,BBE_SELI_PACK_2);
            k=BBE_ADDMOD16U(k,mod); t=BBE_LPNX16_X(R,k); RN=BBE_SELNX16I(RN,t,BBE_SELI_PACK_2);
            RN0=BBE_SELNX16(AN,RN,vcombine0);
            RN1=BBE_SELNX16(AN,RN,vcombine1);
            R =(const int16_t*)XT_ADDX4(W*N,(uintptr_t)R);
            // load diagonals
            BBE_LPNX16_XP(D0,pD,4*SN);
            t=BBE_REPNX16(D0,1);
            D1=BBE_MOVVSV(t,0);
            D0=BBE_REPNX16(D0,0);

            // calculate y(m,:)-R(m,:)*X, 1xP
            BBE_LPNX16_XP(YP,y,4*SN);
            W0=BBE_MULUSNX16(cqXRY,YP);
            vx=BBE_LA_PP(pX);
            BBE_LAVNX16_XP(XN0,vx,pX,M*4);
            BBE_LAVNX16_XP(XN1,vx,pX,M*4-2*BBE_SIMD_WIDTH);
            pX =(const xb_vecNx16*)XT_ADDX4(-M+SN,(uintptr_t)pX);
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
            BBE_SPNX16_XP(t,x,4*SN);
        }
        // goto to the next elements
        idx=XT_ADDX4(-W,idx);   
        Yt-=2;
        Xt-=2;
        D -=2;
        _pX =(const xb_vecNx16*)XT_ADD((uintptr_t)_pX,-4);
    }
    return 0;
}
#endif

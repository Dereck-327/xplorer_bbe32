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
  NatureDSP_Baseband library. Banded Cholesky decomposition for a complex-valued pseudo-inversion:
    Apply the Cholesky decomposition to the matrix of normal equations system
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

#define RADD_BEGIN(W0,AN1,AN0,RN1,RN0,AN,RN,pAk,pRk,nbytes,inc) \
{                                                     \
    xb_vecNx40 W1;                                    \
    xb_c40 I0;                                        \
    BBE_LAVNX16_XP(AN,va,pAk,nbytes);                 \
    BBE_LAVNX16_XP(RN,vr,pRt,nbytes);                 \
    BBE_SELPCNX16I(AN1,AN0,ZERO,AN,0);                \
    BBE_SELPCNX16I(RN1,RN0,ZERO,RN,0);                \
    nbytes=XT_ADDX4(inc,nbytes);                      \
    W1=BBE_MULNX16J(AN,AN); BBE_MULSNX16J(W1,RN,RN);  \
    I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0); \
    W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2); \
}

#define RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN,RN)   \
{                                                     \
    xb_vecNx16 AK,RK;                                 \
    xb_vecNx40 W1;                                    \
    xb_c40 I0;                                        \
    BBE_LAVNX16_XP(AK,va,pAk,nbytes);                 \
    BBE_LAVNX16_XP(RK,vr,pRt,nbytes);                 \
    nbytes=XT_ADDX4(inc,nbytes);                      \
    W1=BBE_MULNX16J(AK,AN); BBE_MULSNX16J(W1,RK,RN);  \
    I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0); \
    W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2); \
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
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right-hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxpn() and cholbkwnxpn(), respectively.

The code for banded matrices is intended for cases where matrix A contains 
W non-zero elements on the main diagonal and below. So, size M is N+W-1. 
Matrix A may be stored in the compact form of size WxN. Cholesky matrix R 
also has WxN non-zero elements

NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)


Input:
  W             Band width
  N             Dimensional parameters
  L             Number of matrices
  sigma2[L]     Regularization term; fixed point position is twice the
                number of fractional bits for matrices A, R
  At[L][SA][2]  Sequence of L complex matrices A represented in the 
                compact form (only band)
Output:
  Rt[L][SR][2]  Sequence of L upper triangular complex matrices R 
                represented in the compact form (saved only elements on 
                the main diagonal and above in such a way that diagonal
                elements are in the last raw)
  D[L][SD][2]   Sequence of L reciprocals of main diagonal A represented 
                in the  block floating point (mantissa and exponent).

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
Note:
  Function may speculatively read up to (W-2)*W complex elements
  beyond the upper bound of At and Rt.
---------------------------------------------------------------------------*/
int bchol8xnn(
            int16_t * restrict Rt, 
            int16_t * restrict D, 
            const int16_t * restrict At, 
            const int32_t * restrict sigma2,
            int N,int L)
#if 0 // variant with convolution using transposed matrices
{
    int cache[8];
    const int W=8;
    valign va,vr;
    xb_vecNx16 t,AN,RN,ZERO=0,D0;
    xb_vecNx40 W0;
    vsaN D1;
    xb_c40 I0;
    int l,Sd=2*getSpace(N);
    int16_t * restrict _D; 
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAk;
          xb_vecNx16* restrict pRk;
          xb_vecNx16* restrict pRk0;
    int n,modinc;
    int idx,M,nbytes,inc,inc1;

    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);

    // clean Rt
    {
        int k;
        pRk=(xb_vecNx16*)Rt;
        for (k=0; k<((W*N*2)/BBE_SIMD_WIDTH)*L; k++) 
        {
            BBE_SVNX16_IP(ZERO,pRk,2*BBE_SIMD_WIDTH);
        }
    }

    W0=0;
    _D=D;
    modinc=((N*W)<<18)|((W-1)<<2);
    pRk0=(xb_vecNx16*)(Rt);
    cache[0]=W*N;
    cache[1]=2*Sd;
    cache[3]=modinc;
    cache[6]=W*N;
    cache[7]=-W;
    for (n=0; n<N; n++)
    {
        // take colunms of A and R and calculate diagonal elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            BBE_LVNX16_XP(AN,pAk,4*W*N);
            BBE_LVNX16_XP(RN,pRt,4*W*N);
            BBE_LPNX16_IP(t,sigma2,4);
            W0=BBE_MOVWVL(t);
            BBE_MULANX16J(W0,AN,AN);
            BBE_MULSNX16J(W0,RN,RN);
            I0=BBE_RADDNX40C(W0);
            W0=BBE_MOVNX40_FROMC40(I0);
            // calculate 1/sqrt(Acc)
            W0 =BBE_ADDNX40(W0,W0);
            D1=BBE_NSAENX40(W0);
            W0=BBE_SLLNX40(W0,D1);
            BBE_RSQRTLUNX40_0(W0,D0, t, W0);
            BBE_MULUUSNX16( W0, t, D0);
            W0=BBE_SRAINX40(W0,23);
            D0=BBE_PACKLNX40(W0);
            D1=BBE_SUBSR1SAVSN(18,D1);
            t=BBE_MOVVVS(D1);
            BBE_SSNX16_I (t,D,2);
            BBE_SSNX16_XP(D0,D,2*Sd);
        }
        sigma2=(const int32_t*)XT_ADDX4(-L,(uintptr_t)sigma2);
        // compute columnar elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        pRk=pRk0;
        M=XT_MIN(W,N-n);
        cache[2]=((W-1)+n*W)<<2;
        cache[4]=W*(N-M);
        cache[5]=(M*W)<<2;
        va=BBE_LA_PP(pAk);
        vr=BBE_LA_PP(pRt);
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            xb_vecNx16 AN0,AN1,RN0,RN1;
            xb_vecNx16 A0,A1,A2,A3,A4,A5,A6,A7;
            xb_vecNx16 R0,R1,R2,R3,R4,R5,R6,R7;
            idx=XT_L32I_N(cache,1*4); BBE_LPNX16_XP(t,D,idx); 
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);
            D1=BBE_ADDSAVSN(1,D1);

            nbytes=XT_L32I_N(cache,5*4);
            inc=XT_L32I_N(cache,7*4);
            BBE_LAVNX16_XP(A0,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A1,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A2,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A3,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A4,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A5,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A6,va,pAk,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(A7,va,pAk,nbytes);
            AN=A0;
            INTLV(A1,A0,A1,A0);
            INTLV(A3,A2,A3,A2);
            INTLV(A5,A4,A5,A4);
            INTLV(A7,A6,A7,A6);
            INTLV(A2,A0,A2,A0);
            INTLV(A3,A1,A3,A1);
            INTLV(A6,A4,A6,A4);
            INTLV(A7,A5,A7,A5);
            INTLV(A4,A0,A4,A0);
            INTLV(A5,A1,A5,A1);
            INTLV(A6,A2,A6,A2);
            INTLV(A7,A3,A7,A3);
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,0);
            W0=BBE_MULRNX16J(A0,AN0,D1);
            BBE_MULANX16J(W0,A1,AN1);
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,2);
            BBE_MULANX16J(W0,A2,AN0);
            BBE_MULANX16J(W0,A3,AN1);
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,4);
            BBE_MULANX16J(W0,A4,AN0);
            BBE_MULANX16J(W0,A5,AN1);
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,6);
            BBE_MULANX16J(W0,A6,AN0);
            BBE_MULANX16J(W0,A7,AN1);
            pAk=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pAk);

            nbytes=XT_L32I_N(cache,5*4);
            inc=XT_L32I_N(cache,7*4);
            BBE_LAVNX16_XP(R0,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R1,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R2,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R3,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R4,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R5,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R6,va,pRt,nbytes);nbytes-=32;
            BBE_LAVNX16_XP(R7,va,pRt,nbytes);
            RN=R0;
            INTLV(R1,R0,R1,R0);
            INTLV(R3,R2,R3,R2);
            INTLV(R5,R4,R5,R4);
            INTLV(R7,R6,R7,R6);
            INTLV(R2,R0,R2,R0);
            INTLV(R3,R1,R3,R1);
            INTLV(R6,R4,R6,R4);
            INTLV(R7,R5,R7,R5);
            INTLV(R4,R0,R4,R0);
            INTLV(R5,R1,R5,R1);
            INTLV(R6,R2,R6,R2);
            INTLV(R7,R3,R7,R3);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,0);
            BBE_MULSNX16J(W0,R0,RN0);
            BBE_MULSNX16J(W0,R1,RN1);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,2);
            BBE_MULSNX16J(W0,R2,RN0);
            BBE_MULSNX16J(W0,R3,RN1);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,4);
            BBE_MULSNX16J(W0,R4,RN0);
            BBE_MULSNX16J(W0,R5,RN1);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,6);
            BBE_MULSNX16J(W0,R6,RN0);
            BBE_MULSNX16J(W0,R7,RN1);
            pRt=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pRt);

            idx=XT_L32I_N(cache,2*4); //idx=((W-1)+n*W)<<2;
            modinc=XT_L32I_N(cache,3*4);

            // normalize wide result
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);
            // element-wise save using circular indexing
            BBE_SELPCNX16I(AN1,AN0,t,t,0);
            BBE_SPNX16_X(  t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SPNX16_X(AN1,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SELPCNX16I(AN1,AN0,t,t,2);
            BBE_SPNX16_X(AN0,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AN1,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            inc1=XT_L32I_N(cache,6*4);
            BBE_SELPCNX16I(AN1,AN0,t,t,4);
            BBE_SPNX16_X(AN0,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AN1,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SELPCNX16I(AN1,AN0,t,t,6);
            BBE_SPNX16_X(AN0,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AN1,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            pRk=(xb_vecNx16*)XT_ADDX4(inc1,(uintptr_t)pRk);
        }
        _D+=2;
        At+=2*W;
        Rt+=2*W;
    }

    return 0;
}
#elif 1
// variant with RADD based convolution and manual sw pipelining
{
    int cache[8];
    const int W=8;
    valign va,vr;
    xb_vecNx16 t,AN,RN,ZERO=0,D0;
    xb_vecNx16 AN0,AN1,RN0,RN1,BN,BK;
    xb_vecNx40 W0;
    vsaN D1;
    int l,Sd=2*getSpace(N);
    int inc,inc1;
    int16_t * restrict _D; 
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAk;
          xb_vecNx16* restrict pRk;
          xb_vecNx16* restrict pRk0;
    int n,modinc;
    int idx,M,nbytes;

    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);

    if (L <= 0) return 0;

    // clean Rt
    {
        int k;
        pRk=(xb_vecNx16*)Rt;
        for (k=0; k<((W*N*2)/BBE_SIMD_WIDTH)*L; k++) 
        {
            BBE_SVNX16_IP(ZERO,pRk,2*BBE_SIMD_WIDTH);
        }
    }

    W0=0;
    _D=D;
    modinc=((N*W)<<18)|((W-1)<<2);
    pRk0=(xb_vecNx16*)(Rt);
    cache[0]=W*N;
    cache[1]=2*Sd;
    cache[3]=modinc;
    cache[6]=W*N;
    cache[7]=-W;
    for (n=0; n<N; n++)
    {
        // take colunms of A and R and calculate diagonal elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            xb_c40 I0;
            BBE_LVNX16_XP(AN,pAk,4*W*N);
            BBE_LVNX16_XP(RN,pRt,4*W*N);

            BBE_LPNX16_IP(t,sigma2,4);
            W0=BBE_MOVWVL(t);
            BBE_MULANX16J(W0,AN,AN);
            BBE_MULSNX16J(W0,RN,RN);
            I0=BBE_RADDNX40C(W0);
            W0=BBE_MOVNX40_FROMC40(I0);
            // calculate 1/sqrt(Acc)
            W0 =BBE_ADDNX40(W0,W0);
            D1=BBE_NSAENX40(W0);
            W0=BBE_SLLNX40(W0,D1);
            BBE_RSQRTLUNX40_0(W0,D0, t, W0);
            BBE_MULUUSNX16( W0, t, D0);
            W0=BBE_SRAINX40(W0,23);
            D0=BBE_PACKLNX40(W0);
            D1=BBE_SUBSR1SAVSN(18,D1);
            t=BBE_MOVVVS(D1);
            BBE_SSNX16_I (t,D,2);
            BBE_SSNX16_XP(D0,D,2*Sd);
        }
        sigma2=(const int32_t*)XT_ADDX4(-L,(uintptr_t)sigma2);
        // compute columnar elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        pRk=pRk0;
        M=XT_MIN(W,N-n);
        cache[2]=((W-1)+n*W)<<2;
        cache[4]=W*(N-M);
        cache[5]=(M*W)<<2;
        va=BBE_LA_PP(pAk);
        vr=BBE_LA_PP(pRt);

        // sw pipelining prologue
        nbytes=XT_L32I_N(cache,5*4);
        inc=XT_L32I_N(cache,7*4);
        RADD_BEGIN(W0,AN1,AN0,RN1,RN0,AN,RN,pAk,pRk,nbytes,inc)
        RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
        BBE_SELPCNX16I(AN1,AN0,ZERO,AN,2);
        BBE_SELPCNX16I(RN1,RN0,ZERO,RN,2);
        RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
        RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)

        #ifdef COMPILER_XTENSA
        #pragma concurrent
        #endif
        for (l=0; l<L-1; l++)
        {
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,4);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,4);
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,6);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,6);
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
            inc=XT_L32I_N(cache,4*4);
            idx=XT_L32I_N(cache,1*4);
            BBE_LPNX16_XP(t,D,idx); 
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);
            pAk=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pRt);

            idx=XT_L32I_N(cache,2*4); //idx=((W-1)+n*W)<<2;
            modinc=XT_L32I_N(cache,3*4);
            
            // normalize wide result
            D1=BBE_ADDSAVSN(1,D1);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);

            // sw pipelining - early begin of next iteration
            nbytes=XT_L32I_N(cache,5*4);
            inc=XT_L32I_N(cache,7*4);
            RADD_BEGIN(W0,AN1,AN0,RN1,RN0,AN,RN,pAk,pRk,nbytes,inc)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,2);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,2);
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)

            // element-wise save using circular indexing
            BBE_SELPCNX16I(BK,BN,t,t,0);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SELPCNX16I(BK,BN,t,t,2);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            inc1=XT_L32I_N(cache,6*4);
            BBE_SELPCNX16I(BK,BN,t,t,4);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SELPCNX16I(BK,BN,t,t,6);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            pRk=(xb_vecNx16*)XT_ADDX4(inc1,(uintptr_t)pRk);
        }
        // sw pipelining - epilogue
        {
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,4);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,4);
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,6);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,6);
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN0,RN0)
            RADD_ITERATION(W0,pAk,pRk,nbytes,inc,AN1,RN1)
            inc=XT_L32I_N(cache,4*4);
            idx=XT_L32I_N(cache,1*4);
            BBE_LPNX16_XP(t,D,idx); 
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);
            pAk=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pRt);

            idx=XT_L32I_N(cache,2*4); //idx=((W-1)+n*W)<<2;
            modinc=XT_L32I_N(cache,3*4);
            
            // normalize wide result
            D1=BBE_ADDSAVSN(1,D1);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);

            // element-wise save using circular indexing
            BBE_SELPCNX16I(BK,BN,t,t,0);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SELPCNX16I(BK,BN,t,t,2);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            inc1=XT_L32I_N(cache,6*4);
            BBE_SELPCNX16I(BK,BN,t,t,4);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SELPCNX16I(BK,BN,t,t,6);
            BBE_SPNX16_X(BN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(BK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
         }
        _D+=2;
        At+=2*W;
        Rt+=2*W;
    }

    return 0;
}
#else
// variant with RADD based convolution
{
    int cache[8];
    const int W=8;
    valign va,vr;
    xb_vecNx16 t,AN,RN,AK,RK,ZERO=0,D0;
    xb_vecNx40 W0,W1;
    vsaN D1;
    xb_c40 I0;
    int l,Sd=2*getSpace(N);
    int16_t * restrict _D; 
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAk;
          xb_vecNx16* restrict pRk;
          xb_vecNx16* restrict pRk0;
    int n,modinc;
    int idx,M,nbytes;

    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);

    // clean Rt
    {
        int k;
        pRk=(xb_vecNx16*)Rt;
        for (k=0; k<((W*N*2)/BBE_SIMD_WIDTH)*L; k++) 
        {
            BBE_SVNX16_IP(ZERO,pRk,2*BBE_SIMD_WIDTH);
        }
    }

    W0=0;
    _D=D;
    modinc=((N*W)<<18)|((W-1)<<2);
    pRk0=(xb_vecNx16*)(Rt);
    cache[0]=W*N;
    cache[1]=2*Sd;
    cache[3]=modinc;
    cache[6]=W*N;
    cache[7]=-W;
    for (n=0; n<N; n++)
    {
        // take colunms of A and R and calculate diagonal elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {/*
            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            BBE_LAVNX16_XP(AN,va,pAk,W*4);
            BBE_LAVNX16_XP(RN,vr,pRt,W*4);
            pAk=(const xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pAk);
            pRt=(      xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pRt);
            */
            BBE_LVNX16_XP(AN,pAk,4*W*N);
            BBE_LVNX16_XP(RN,pRt,4*W*N);

            BBE_LPNX16_IP(t,sigma2,4);
            W0=BBE_MOVWVL(t);
            BBE_MULANX16J(W0,AN,AN);
            BBE_MULSNX16J(W0,RN,RN);
            I0=BBE_RADDNX40C(W0);
            W0=BBE_MOVNX40_FROMC40(I0);
            // calculate 1/sqrt(Acc)
            W0 =BBE_ADDNX40(W0,W0);
            D1=BBE_NSAENX40(W0);
            W0=BBE_SLLNX40(W0,D1);
            BBE_RSQRTLUNX40_0(W0,D0, t, W0);
            BBE_MULUUSNX16( W0, t, D0);
            W0=BBE_SRAINX40(W0,23);
            D0=BBE_PACKLNX40(W0);
            D1=BBE_SUBSR1SAVSN(18,D1);
            t=BBE_MOVVVS(D1);
            BBE_SSNX16_I (t,D,2);
            BBE_SSNX16_XP(D0,D,2*Sd);
        }
        sigma2=(const int32_t*)XT_ADDX4(-L,(uintptr_t)sigma2);
        // compute columnar elements
        D = _D;
        pAk=(const xb_vecNx16*)(At);
        pRt=(const xb_vecNx16*)(Rt);
        pRk=pRk0;
        M=XT_MIN(W,N-n);
        cache[2]=((W-1)+n*W)<<2;
        cache[4]=W*(N-M);
        cache[5]=(M*W)<<2;

        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            int inc;
            xb_vecNx16 AN0,AN1,RN0,RN1;
            idx=XT_L32I_N(cache,1*4);BBE_LPNX16_XP(t,D,idx); //BBE_LPNX16_XP(t,D,2*Sd);
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);

            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            //nbytes=(M*W)<<2;
            nbytes=XT_L32I_N(cache,5*4);
            inc=XT_L32I_N(cache,7*4);
            {
                BBE_LAVNX16_XP(AN,va,pAk,nbytes);
                BBE_LAVNX16_XP(RN,vr,pRt,nbytes);
                BBE_SELPCNX16I(AN1,AN0,ZERO,AN,0);
                BBE_SELPCNX16I(RN1,RN0,ZERO,RN,0);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AN,AN); BBE_MULSNX16J(W1,RN,RN);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN1); BBE_MULSNX16J(W1,RK,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,2);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,2);
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN0); BBE_MULSNX16J(W1,RK,RN0);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN1); BBE_MULSNX16J(W1,RK,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,4);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,4);
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN0); BBE_MULSNX16J(W1,RK,RN0);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN1); BBE_MULSNX16J(W1,RK,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            BBE_SELPCNX16I(AN1,AN0,ZERO,AN,6);
            BBE_SELPCNX16I(RN1,RN0,ZERO,RN,6);
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN0); BBE_MULSNX16J(W1,RK,RN0);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK,vr,pRt,nbytes);
                nbytes=XT_ADDX4(inc,nbytes);
                W1=BBE_MULNX16J(AK,AN1); BBE_MULSNX16J(W1,RK,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            inc=XT_L32I_N(cache,4*4);
            pAk=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pRt);

            idx=XT_L32I_N(cache,2*4); //idx=((W-1)+n*W)<<2;
            modinc=XT_L32I_N(cache,3*4);
            
            // normalize wide result
            D1=BBE_ADDSAVSN(1,D1);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);
            // element-wise save using circular indexing
            BBE_SELPCNX16I(AK,AN,t,t,0);
            BBE_SPNX16_X(AN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SPNX16_X(AK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            BBE_SELPCNX16I(AK,AN,t,t,2);
            BBE_SPNX16_X(AN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            inc=XT_L32I_N(cache,6*4);
            BBE_SELPCNX16I(AK,AN,t,t,4);
            BBE_SPNX16_X(AN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SELPCNX16I(AK,AN,t,t,6);
            BBE_SPNX16_X(AN,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            BBE_SPNX16_X(AK,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); 
            pRk=(xb_vecNx16*)XT_ADDX4(inc,(uintptr_t)pRk);
        }
        _D+=2;
        At+=2*W;
        Rt+=2*W;
    }

    return 0;
}
#endif

#endif

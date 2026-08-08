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

int bchol4xnn( int16_t * restrict Rt, 
               int16_t * restrict D, 
         const int16_t * restrict At, 
         const int32_t * restrict sigma2,
               int N,int L)
#if 0
{
    const int W=4;
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
    vselN rotateW;
    int n,modinc;
    int idx,M,nbytes;

    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);
    // prepare select that rotates right by W elements
    t=BBE_SEQNX16();
    RN=BBE_MOVVA16(W*2);
    t=BBE_SUBNX16(t,RN);
    rotateW=BBE_MOVVSELNX16(t,0);

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
            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            BBE_LAVNX16_XP(AN,va,pAk,W*4);
            BBE_LAVNX16_XP(RN,vr,pRt,W*4);
            pAk=(const xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pAk);
            pRt=(      xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pRt);
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
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            BBE_LPNX16_XP(t,D,2*Sd);
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);

            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            M=XT_MIN(W,N-n);
            nbytes=(M*W)<<2;
            {
                BBE_LAVNX16_XP(AN,va,pAk,XT_MIN(nbytes,W*4));
                BBE_LAVNX16_XP(RN,vr,pRt,XT_MIN(nbytes,W*4));
                nbytes=XT_ADDX4(-W,nbytes);
                W1=BBE_MULNX16J(AN,AN); BBE_MULSNX16J(W1,RN,RN);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN=BBE_SELNX16I(ZERO,AN,BBE_SELI_ROTATE_RIGHT_2);
                RN=BBE_SELNX16I(ZERO,RN,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,XT_MIN(nbytes,W*4));
                BBE_LAVNX16_XP(RK,vr,pRt,XT_MIN(nbytes,W*4));
                nbytes=XT_ADDX4(-W,nbytes);
                W1=BBE_MULNX16J(AK,AN); BBE_MULSNX16J(W1,RK,RN);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN=BBE_SELNX16I(ZERO,AN,BBE_SELI_ROTATE_RIGHT_2);
                RN=BBE_SELNX16I(ZERO,RN,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,XT_MIN(nbytes,W*4));
                BBE_LAVNX16_XP(RK,vr,pRt,XT_MIN(nbytes,W*4));
                nbytes=XT_ADDX4(-W,nbytes);
                W1=BBE_MULNX16J(AK,AN); BBE_MULSNX16J(W1,RK,RN);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN=BBE_SELNX16I(ZERO,AN,BBE_SELI_ROTATE_RIGHT_2);
                RN=BBE_SELNX16I(ZERO,RN,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            {
                BBE_LAVNX16_XP(AK,va,pAk,XT_MIN(nbytes,W*4));
                BBE_LAVNX16_XP(RK,vr,pRt,XT_MIN(nbytes,W*4));
                nbytes=XT_ADDX4(-W,nbytes);
                W1=BBE_MULNX16J(AK,AN); BBE_MULSNX16J(W1,RK,RN);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN=BBE_SELNX16I(ZERO,AN,BBE_SELI_ROTATE_RIGHT_2);
                RN=BBE_SELNX16I(ZERO,RN,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_SELNX40I(W1,W0,BBE_W_SELI_ROTATE_RIGHT_2);
            }
            pAk=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pRt);

            // normalize wide result
            D1=BBE_ADDSAVSN(1,D1);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);
            t=BBE_SHFLNX16(t,rotateW);
            // element-wise save using circular indexing
            idx=((W-1)+n*W)<<2;
            BBE_SPNX16_X(t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); t=BBE_SELNX16I(t,t,BBE_SELI_ROTATE_RIGHT_2);
            BBE_SPNX16_X(t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); t=BBE_SELNX16I(t,t,BBE_SELI_ROTATE_RIGHT_2);
            BBE_SPNX16_X(t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); t=BBE_SELNX16I(t,t,BBE_SELI_ROTATE_RIGHT_2);
            BBE_SPNX16_X(t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); t=BBE_SELNX16I(t,t,BBE_SELI_ROTATE_RIGHT_2);
            pRk=(xb_vecNx16*)XT_ADDX4(W*N,(uintptr_t)pRk);
        }
        _D+=2;
        At+=2*W;
        Rt+=2*W;
    }

    return 0;
}
#else
{
    static const int16_t ALIGN(32) sel[16]=
    {
        0, 1, 8, 9, 16,17,24,25,
        0, 1, 8, 9, 16,17,24,25
    };
    const int W=4;
    int cache[4];
    valign va,vr;
    xb_vecNx16 t,AN,RN,AK,RK,ZERO=0,D0;
    xb_vecNx40 W0;
    vsaN D1;
    xb_c40 I0;
    int l,Sd=2*getSpace(N);
          int16_t   * restrict _D; 
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAk;
          xb_vecNx16* restrict pRk;
          xb_vecNx16* restrict pRk0;
    const xb_vecNx16* restrict tempptr;
    vselN selu0,selu;
    int n,modinc;
    int idx,M,nbytes;

    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(At,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
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

    t=BBE_LVNX16_I((const xb_vecNx16*)sel,0);
    selu0=BBE_MOVVSELNX16(t,0);
    cache[0]=W*N;
    cache[1]=2*Sd;
    cache[3]=modinc;
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
            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            BBE_LAVNX16_XP(AN,va,pAk,W*4);
            BBE_LAVNX16_XP(RN,vr,pRt,W*4);
            pAk=(const xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pAk);
            pRt=(      xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pRt);
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
        nbytes=(M*W)<<2;
        cache[2]=((W-1)+n*W)<<2;
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=1
        #pragma concurrent
        #endif
        for (l=0; l<L; l++)
        {
            xb_vecNx16 x0,z0,y0,x3,x2,z3,z2;
            idx=XT_L32I(cache,1*4);BBE_LPNX16_XP(t,D,idx); //BBE_LPNX16_XP(t,D,2*Sd);
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);

            tempptr=pAk; va=BBE_LA_PP(tempptr); BBE_LAVNX16_XP(AN,va,tempptr,nbytes);BBE_LAVNX16_XP(AK,va,tempptr,nbytes-32);
            tempptr=pRt; vr=BBE_LA_PP(tempptr); BBE_LAVNX16_XP(RN,vr,tempptr,nbytes);BBE_LAVNX16_XP(RK,vr,tempptr,nbytes-32);
            idx=XT_L32I(cache,1*4);
            pAk=(const xb_vecNx16*)XT_ADDX4(idx,(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(idx,(uintptr_t)pRt);
            x0=BBE_SELNX16I(ZERO,AN,BBE_SELI_PACK_8);   // xxxx0000
            z0=BBE_SELNX16I(ZERO,RN,BBE_SELI_PACK_8);   // xxxx0000
            BBE_SELPCNX16I(x3, x2, ZERO, x0, 2);
            BBE_SELPCNX16I(z3, z2, ZERO, z0, 2);
            selu=selu0;
            y0=BBE_SELNX16(AK,AN,selu); W0=BBE_MULNX16J(y0,x0); 
            BBE_SELUNX16(y0,RK,RN,selu,2); BBE_MULSNX16J(W0,y0,z0); 
            x0=BBE_SELNX16I(ZERO,x0,BBE_SELI_ROTATE_RIGHT_2);// xxx00000
            z0=BBE_SELNX16I(ZERO,z0,BBE_SELI_ROTATE_RIGHT_2);// xxx00000
            y0=BBE_SELNX16(AK,AN,selu);    BBE_MULANX16J(W0,y0,x0); 
            BBE_SELUNX16(y0,RK,RN,selu,2); BBE_MULSNX16J(W0,y0,z0); 
            y0=BBE_SELNX16(AK,AN,selu);    BBE_MULANX16J(W0,y0,x2); 
            BBE_SELUNX16(y0,RK,RN,selu,2); BBE_MULSNX16J(W0,y0,z2); 
            y0=BBE_SELNX16(AK,AN,selu);    BBE_MULANX16J(W0,y0,x3); 
            BBE_SELUNX16(y0,RK,RN,selu,2); BBE_MULSNX16J(W0,y0,z3); 

            // normalize wide result
            D1=BBE_ADDSAVSN(1,D1);
            W0=BBE_RNDADJNX40(W0,D1);
            t=BBE_PACKVNX40(W0,D1);
            W0=BBE_MULUSNX16(D0,t);
            t=BBE_PACKQNX40(W0);

            // element-wise save using circular indexing
            idx=XT_L32I(cache,2*4); //idx=((W-1)+n*W)<<2;
            modinc=XT_L32I(cache,3*4);
            BBE_SPNX16_X(t ,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); x0=BBE_SHFLNX16I(t,BBE_SHFLI_DUPLICATE_2_ODD); // 2,3
            BBE_SPNX16_X(x0,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); x0=BBE_SHFLNX16I(t,BBE_SHFLI_MMC4X4X4X4_M1_STEP_3); // 4,5
            BBE_SPNX16_X(x0,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc); x0=BBE_SHFLNX16I(t,BBE_SHFLI_MMC4X4X4X4_M1_STEP_4); // 6,7 
            BBE_SPNX16_X(x0,pRk,idx);  
            //pRk=(xb_vecNx16*)XT_ADDX4(W*N,(uintptr_t)pRk);
            idx=XT_L32I(cache,1*4);
            pRk=(xb_vecNx16*)XT_ADDX4(idx,(uintptr_t)pRk);
        }
        _D+=2;
        At+=2*W;
        Rt+=2*W;
    }

    return 0;
}
#endif

#endif

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

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the streaming order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxps() and cholbkwnxps(), respectively.

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
  At[L][SA]     Sequence of L complex matrices A represented in the 
                compact form (only band)
Output:
  Rt[L][SR]     Sequence of L upper triangular complex matrices R 
                represented in the compact form (saved only elements on 
                the main diagonal and above in such a way that diagoanal
                elements are in the last raw)
  D[L][SD]      Sequence of L reciprocals of main diagonal A represented 
                in the  block floating point (mantissa and exponent).

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
Note:
  Function may speculatively read up to (W-2)*W complex elements
  beyond the upper bound of At and Rt.
---------------------------------------------------------------------------*/
#if !HAVE_BCHOLN
DISCARD_FUN(void,bchol16x32n,(
                  complex_fract16 * restrict Rt, 
                  complex_fract16 * restrict D, 
            const complex_fract16 * restrict At, 
            const int32_t * restrict sigma2,
            int L))
#else
void bchol16x32n (complex_fract16 * restrict Rt, complex_fract16 * restrict __D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L)
{
    int16_t * restrict D=(int16_t *)__D;
#define W 16
#define N 32
    valign va,vr;
    xb_vecNx16 t,AN0,RN0,AN1,RN1,AK0,RK0,AK1,RK1,ZERO=0,D0;
    xb_vecNx40 W0,W1;
    vsaN D1;
    xb_c40 I0;
    int l,Sd=2*N;
    int16_t * restrict _D; 
    const xb_vecNx16* restrict pRt;
    const xb_vecNx16* restrict pAk;
          xb_vecNx16* restrict pRk;
          xb_vecNx16* restrict pRk0;
    int n,modinc;
    int idx,M,nbytes;

    NASSERT_ALIGN(Rt,32);
    NASSERT_ALIGN(At,32);
    NASSERT((N&3)==0);
    NASSERT((W&3)==0);
    NASSERT(W>=4 && W<=BCHOLN_MAXW);
    
    if (L <= 0) return;
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
    _D=(int16_t*)D;
    modinc=((N*W)<<18)|((W-1)<<2);
    pRk0=(xb_vecNx16*)(Rt);
    #ifdef COMPILER_XTENSA
    #pragma loop_count min=4
    #endif
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
            BBE_LANX16_IP(AN0,va,pAk);
            BBE_LANX16_IP(AN1,va,pAk);
            BBE_LANX16_IP(RN0,vr,pRt);
            BBE_LANX16_IP(RN1,vr,pRt);
            pAk=(const xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pAk);
            pRt=(      xb_vecNx16*)XT_ADDX4(W*N-W,(uintptr_t)pRt);
            BBE_LPNX16_IP(t,sigma2,4);
            W0=BBE_MOVWVL(t);
            BBE_MULANX16J(W0,AN0,AN0);
            BBE_MULANX16J(W0,AN1,AN1);
            BBE_MULSNX16J(W0,RN0,RN0);
            BBE_MULSNX16J(W0,RN1,RN1);
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
            int k,mbytes;
            BBE_LPNX16_XP(t,D,2*Sd);
            D0=BBE_REPNX16(t,0);
            t =BBE_REPNX16(t,1);
            D1=BBE_MOVVSV(t,0);
            D1=BBE_ADDSAVSN(1,D1);

            va=BBE_LA_PP(pAk);
            vr=BBE_LA_PP(pRt);
            M=XT_MIN(W,N-n);
            nbytes=(M*W)<<2;
            mbytes=nbytes-2*BBE_SIMD_WIDTH;
            idx=((W-1)+n*W)<<2;
            {
                BBE_LAVNX16_XP(AN0,va,pAk,nbytes);
                BBE_LAVNX16_XP(RN0,vr,pRt,nbytes);
                BBE_LAVNX16_XP(AN1,va,pAk,mbytes);
                BBE_LAVNX16_XP(RN1,vr,pRt,mbytes);
                nbytes=XT_ADDX4(-W,nbytes);
                mbytes=XT_ADDX4(-W,mbytes);
                W1=BBE_MULNX16J(AN0,AN0); 
                BBE_MULANX16J(W1,AN1,AN1);
                BBE_MULSNX16J(W1,RN0,RN0);
                BBE_MULSNX16J(W1,RN1,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN0=BBE_SELNX16I( AN1,AN0,BBE_SELI_ROTATE_RIGHT_2);
                AN1=BBE_SELNX16I(ZERO,AN1,BBE_SELI_ROTATE_RIGHT_2);
                RN0=BBE_SELNX16I( RN1,RN0,BBE_SELI_ROTATE_RIGHT_2);
                RN1=BBE_SELNX16I(ZERO,RN1,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_RNDADJNX40(W1,D1);
                t=BBE_PACKVNX40(W0,D1);
                W0=BBE_MULUSNX16(D0,t);
                t=BBE_PACKQNX40(W0);
                // element-wise save using circular indexing
                BBE_SPNX16_X(  t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            }
            #ifdef COMPILER_XTENSA
            #pragma concurrent
            #endif
            for (k=1; k<M; k++)
            {
                BBE_LAVNX16_XP(AK0,va,pAk,nbytes);
                BBE_LAVNX16_XP(RK0,vr,pRt,nbytes);
                BBE_LAVNX16_XP(AK1,va,pAk,mbytes);
                BBE_LAVNX16_XP(RK1,vr,pRt,mbytes);
                nbytes=XT_ADDX4(-W,nbytes);
                mbytes=XT_ADDX4(-W,mbytes);
                W1=BBE_MULNX16J(AK0,AN0); 
                BBE_MULANX16J(W1,AK1,AN1); 
                BBE_MULSNX16J(W1,RK0,RN0);
                BBE_MULSNX16J(W1,RK1,RN1);
                I0=BBE_RADDNX40C(W1); W1=BBE_MOVNX40_FROMC40(I0);
                AN0=BBE_SELNX16I( AN1,AN0,BBE_SELI_ROTATE_RIGHT_2);
                AN1=BBE_SELNX16I(ZERO,AN1,BBE_SELI_ROTATE_RIGHT_2);
                RN0=BBE_SELNX16I( RN1,RN0,BBE_SELI_ROTATE_RIGHT_2);
                RN1=BBE_SELNX16I(ZERO,RN1,BBE_SELI_ROTATE_RIGHT_2);
                W0=BBE_RNDADJNX40(W1,D1);
                t=BBE_PACKVNX40(W0,D1);
                W0=BBE_MULUSNX16(D0,t);
                t=BBE_PACKQNX40(W0);
                // element-wise save using circular indexing
                BBE_SPNX16_X(  t,pRk,idx); idx=BBE_ADDMOD16U(idx,modinc);
            }
            pAk=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pAk);
            pRt=(const xb_vecNx16*)XT_ADDX4(W*(N-M),(uintptr_t)pRt);
            pRk=(xb_vecNx16*)XT_ADDMI((uintptr_t)pRk,4*W*N);
        }
        _D+=2;
        At+=W;
        Rt+=W;
    }

#undef N
#undef W
} /* bchol16x32n() */
#endif

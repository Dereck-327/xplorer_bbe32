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
  NatureDSP_Baseband library. Apply the QR decomposition to the matrix of normal equations system
    Apply backward recursion process for QR decomposition for block ordered
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqrn_common.h"

#if !HAVE_CQRN
DISCARD_FUN(void,cqr_bkwmxnxpn,(void *pScr,
                          complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int M, int N, int P,
                    int L))
#else
// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-1)) m=LOG2_BBE_SIMD_WIDTH-1;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*
    backward recursion: P!=1
*/
static void cqrnBkwnxp(int16_t* restrict X, 
                 const int16_t* restrict R,
                 const int16_t* restrict D,
                       int qABX,
                       int M,int N,int P, int L)
{
    NASSERT_ALIGN(X   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D   ,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0);

    const int SR=2*getSpace(M*N);
    const int SX=2*getSpace(M*P);
    const int SD=2*getSpace(N);

    int k;
    int l;
    int m;
    int p;

    int delta= 0;
    
    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    vsaN q;

    const xb_vecNx16* restrict pXrd_first;
    const xb_vecNx16* restrict pXrd_first_tmp;

    const xb_vecNx16* restrict pXrd;
    const xb_vecNx16* restrict pXrd_tmp;

    xb_vecNx16* restrict pXwr;
    xb_vecNx16* restrict pXwr_tmp;

    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    xb_vecNx16 xx;
    xb_vecNx16 rr;

    valign Xrd_first_align, Xrd_align, Xwr_align;

    const vsaN sh16= BBE_MOVVSA32(16);
    const vsaN vqABX= BBE_MOVVSA32(15-qABX);

    Xwr_align= BBE_ZALIGN();

    __Pragma("loop_count min=1")
    for (k=N-1; k>=0; k--)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for(p= P; p>0; p-=8)
        {
            pXrd_first= (const xb_vecNx16*)XT_ADDX2(delta,(uintptr_t)&X[2*k*P]);
            pXrd= (const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pXrd_first);
            pXwr= (xb_vecNx16*)pXrd_first;
            pRrd= (const xb_vecNx16*)(R+(k*N+(k+1))*2);
            //pDrd= (const xb_vecNx16*)&D[k+k];
            pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
            delta= XT_ADDI(delta, 2*8);
            __Pragma("loop_count min=1")
            for(l=0; l<L; l++)
            {
                // calculate y(m,:)-R(m,:)*X, 1xP
                // load 16 bit complex value (16 bit re and 16 bit im)
                pXrd_first_tmp= pXrd_first;
                Xrd_first_align= BBE_LA_PP(pXrd_first_tmp);
                BBE_LAVNX16_XP(b,Xrd_first_align,pXrd_first_tmp,p*4);
                pXrd_first= (const xb_vecNx16*)XT_ADDX2(SX, (uintptr_t)pXrd_first);

                B=BBE_UNPKNVNX16(b,vqABX);

                for (m=0; m<N-k-1; m++)
                {
                    pXrd_tmp= pXrd;
                    Xrd_align= BBE_LA_PP(pXrd_tmp);
                    BBE_LAVNX16_XP(xx,Xrd_align,pXrd_tmp,p*4);
                    pXrd= (const xb_vecNx16*)XT_ADDX2(2*P, (uintptr_t)pXrd);

                    BBE_LPNX16_IP(rr,pRrd,2*2);
                    rr= BBE_REPNX16C(rr,0);
                    BBE_MULSNX16C(B,xx,rr);
                }
                pXrd= (const xb_vecNx16*)XT_ADDX2(SX - (N-k-1)*2*P,(uintptr_t)pXrd);
                pRrd= (const xb_vecNx16*)XT_ADDX2(SR - (N-k-1)*2,(uintptr_t)pRrd);
                // get low 16 bits of B
                tl = BBE_PACKLNX40(B);
                // get high 16 bits of B
                th = BBE_PACKVNX40(B,sh16);
                // load D
                BBE_LPNX16_XP(dd,pDrd,2*SD);
                // replicate dd[0]
                d0 = BBE_REPNX16(dd,0);
                // mul low (unsigned*unsigned)
                B = BBE_MULUUNX16(d0,tl);
                // shift low left
                B = BBE_SRAINX40(B,16);
                // Bres+= mul high (unsigned*signed)
                BBE_MULUSANX16(B,d0,th);
                // make q
                d0= BBE_REPNX16(dd,1);
                q = BBE_MOVVSV(d0,0);
                // result rounding
                B = BBE_RNDADJNX40(B,q);
                // store res
                b_res = BBE_PACKVNX40(B,q);
                pXwr_tmp= pXwr;
                BBE_SAVNX16_XP(b_res,Xwr_align,pXwr_tmp,p*4);
                BBE_SAPOS_FP(Xwr_align,pXwr_tmp);
                pXwr= (xb_vecNx16*)XT_ADDX2(SX, (uintptr_t)pXwr);
            }
        }
    }
}

/*-------------------------------------------------------------------------
Apply backward recursion process for QR decomposition for block ordered 
matrices.
Matrix sizes SA,SB are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SB=size(M*P)
SD=size(N)
Scratch size in bytes is defined by cqr_bkwmxnxpn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      Dimensional parameters
 L            Number of matrices
 qABX         qA-qB+qX where qA,qB,qX - fixed point representations of 
              matrices A,B,X
Input/output:
 X[L][SB]     On input it is the sequence of L updated right parts Z=Q'B.
              They will be replaced with MMSE solution vectors X (only N*P 
              elements are used)
Input:
 R[L][SA]     Upper triangular matrices R (only N*N 
              elements of each matrix are used)
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format

Restrictions:
1. X, R, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/

void cqr_bkwmxnxpn (void *pScr,
                          complex_fract16* _X,
                    const complex_fract16* _R,
                    const complex_fract16* _D,
                    int qABX,
                    int M, int N, int P,
                    int L)
{
          int16_t* X=(      int16_t*)_X;
    const int16_t* R=(const int16_t*)_R;
    const int16_t* D=(const int16_t*)_D;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(X   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D   ,2*BBE_SIMD_WIDTH);
    NASSERT(M%4==0 && M>0);
    NASSERT(N%4==0 && N>0);
    NASSERT(L>0);
    if (P<=0) return;
    if (P==1)  cqrnBkwnx1(X,R,D,qABX,M,N,L);
    else       cqrnBkwnxp(X,R,D,qABX,M,N,P,L);
} /* cqr_bkwmxnxpn() */
#endif

size_t cqr_bkwmxnxpn_getScratchSize (int M, int N,int P,int L)
{
    (void)M;
    (void)N;
    (void)P;
    (void)L;
    return 0;
} /* cqr_bkwmxnxpn_getScratchSize() */

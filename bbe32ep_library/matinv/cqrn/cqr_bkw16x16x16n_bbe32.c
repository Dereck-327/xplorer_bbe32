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
DISCARD_FUN(void,cqr_bkw16x16x16n,(void *pScr,
                          complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L))
#else

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

void cqr_bkw16x16x16n (void *pScr,
                          complex_fract16* _X,
                    const complex_fract16* _R,
                    const complex_fract16* _D,
                    int qABX,
                    int L)
{
          int16_t* X=(      int16_t*)_X;
    const int16_t* R=(const int16_t*)_R;
    const int16_t* D=(const int16_t*)_D;
    #define N 16
    #define M 16
    #define P 16
    #define SR 512
    #define SX 512
    #define SD 32
    
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(X   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D   ,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    int k;
    int l;
    int m;
    
    xb_vecNx16 dd, d0, d1;
    xb_vecNx16 tl0, th0, tl1, th1, b_res0, b_res1, b0, b1;
    xb_vecNx40 B0, B1;
    vsaN q;

    const xb_vecNx16 * restrict pXrd_first;
    const xb_vecNx16 * restrict pXrd;
          xb_vecNx16 * restrict pXwr;
    const xb_vecNx16 * restrict pRrd;
    const xb_vecNx16 * restrict pDrd;

    xb_vecNx16 xx0, xx1;
    xb_vecNx16 rr;

    const vsaN sh16 = BBE_MOVVSA32(16);
    const vsaN vqABX= BBE_MOVVSA32(qABX);

    __Pragma("loop_count min=1")
    for (k=N-1; k>=0; k--)
    {
        pXrd_first= (const xb_vecNx16*)XT_ADDX4(k*P,(uintptr_t)X);
        pXrd= (const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pXrd_first);
        pXwr= (xb_vecNx16*)pXrd_first;
        pRrd= (const xb_vecNx16*)(R+(k*N+(k+1))*2);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);

        __Pragma("loop_count min=1")
        for(l=0; l<L; l++)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP
            // load P(=16) 16 bit complex value (16 bit re and 16 bit im)
            BBE_LVNX16_IP(b0, pXrd_first, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(b1, pXrd_first, 2*SX-2*BBE_SIMD_WIDTH);
            B0 = BBE_UNPKSNX16(b0);
            B1 = BBE_UNPKSNX16(b1);
            B0 = BBE_SLSNX40(B0, vqABX);
            B1 = BBE_SLSNX40(B1, vqABX);

            for (m=0; m<N-k-1; m++)
            {
                BBE_LVNX16_IP(xx0,pXrd,2*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(xx1,pXrd,2*BBE_SIMD_WIDTH);
                BBE_LPNX16_IP(rr,pRrd,2*2);
                rr= BBE_REPNX16C(rr,0);
                BBE_MULSNX16C(B0,xx0,rr);
                BBE_MULSNX16C(B1,xx1,rr);
            }
            pXrd= (const xb_vecNx16*)XT_ADDX2(SX - (N-k-1)*2*P,(uintptr_t)pXrd);
            pRrd= (const xb_vecNx16*)XT_ADDX2(SR - (N-k-1)*2  ,(uintptr_t)pRrd);
            // load D
            BBE_LPNX16_XP(dd,pDrd,2*SD);
            // replicate dd[0]
            d0 = BBE_REPNX16(dd,0);
            // make q
            d1= BBE_REPNX16(dd,1);
            q = BBE_MOVVSV(d1,0);
            //------------------------------------------------------------------
            // get low 16 bits of B
            tl0 = BBE_PACKLNX40(B0);
            tl1 = BBE_PACKLNX40(B1);
            // get high 16 bits of B
            th0 = BBE_PACKVNX40(B0,sh16);
            th1 = BBE_PACKVNX40(B1,sh16);
            // mul low (unsigned*unsigned)
            B0 = BBE_MULUUNX16(d0,tl0);
            B1 = BBE_MULUUNX16(d0,tl1);
            // shift low left
            B0 = BBE_SRAINX40(B0,16);
            B1 = BBE_SRAINX40(B1,16);
            // Bres+= mul high (unsigned*signed)
            BBE_MULUSANX16(B0,d0,th0);
            BBE_MULUSANX16(B1,d0,th1);
            // result rounding
            B0 = BBE_RNDADJNX40(B0,q);
            B1 = BBE_RNDADJNX40(B1,q);
            b_res0 = BBE_PACKVNX40(B0,q);
            b_res1 = BBE_PACKVNX40(B1,q);
            // store res
            BBE_SVNX16_IP(b_res0,pXwr,2*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(b_res1,pXwr,2*SX-2*BBE_SIMD_WIDTH);
            //-------------------------------------------------------------------
        }
    }
    #undef N
    #undef M
    #undef P
    #undef SR
    #undef SX
    #undef SD
} /* cqr_bkw16x16x16n() */
#endif

size_t cqr_bkw16x16x16n_getScratchSize (int M, int N,int P,int L)
{
    (void)M;
    (void)N;
    (void)P;
    (void)L;
    return 0;
} /* cqr_bkw16x16x16n_getScratchSize() */

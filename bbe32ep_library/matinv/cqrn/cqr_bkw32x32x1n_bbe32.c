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
DISCARD_FUN(void,cqr_bkw32x32x1n,(void *pScr,
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

void cqr_bkw32x32x1n (void *pScr,
                          complex_fract16* _X,
                    const complex_fract16* _R,
                    const complex_fract16* _D,
                    int qABX,
                    int L)
{
          int16_t* X=(      int16_t*)_X;
    const int16_t* R=(const int16_t*)_R;
    const int16_t* D=(const int16_t*)_D;
    #define SR 2048
    #define SX 64
    #define SD 64
    #define N 32

    NASSERT_ALIGN(X   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D   ,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0);

    int k;
    int l;

    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;

    valign R_align;
    valign X_align;

    const xb_vecNx16* restrict pXrd_first;
    const xb_vecNx16* restrict pXrd;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    xb_vecNx16 xx0, xx1, xx2, xx3;
    xb_vecNx16 rr0, rr1, rr2, rr3;

    const vsaN sh16= BBE_MOVVSA32(16);
    const vsaN vqABX= BBE_MOVVSA32(15-qABX);

    for (k=N-1; k>=0; k--)
    {
        pXrd_first= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)X);
        pXrd= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXrd_first,2*2);
        pRrd= (const xb_vecNx16*)(R + (2*k*N+(k+k+2)));
        pXwr= (xb_vecNx16*)pXrd_first;
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);

        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP
            // load 16 bit complex value (16 bit re and 16 bit im)
            BBE_LPNX16_XP(b,pXrd_first,SX*2);
            B= BBE_UNPKNVNX16(b,vqABX);

            // load xx (load only N-k-1 elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx0,X_align,pXrd,(N-k-1)*2*2);
            BBE_LAVNX16_XP(xx1,X_align,pXrd,(N-k-1)*2*2 - 2*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(xx2,X_align,pXrd,(N-k-1)*2*2 - 4*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(xx3,X_align,pXrd,(N-k-1)*2*2 - 6*BBE_SIMD_WIDTH);
            //pXrd++
            pXrd=(const xb_vecNx16*)XT_ADDX2(SX-((N-k-1)*2),(uintptr_t)pXrd);

            // load rr (load only N-k-1 elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr0,R_align,pRrd,(N-k-1)*2*2);
            BBE_LAVNX16_XP(rr1,R_align,pRrd,(N-k-1)*2*2 - 2*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(rr2,R_align,pRrd,(N-k-1)*2*2 - 4*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(rr3,R_align,pRrd,(N-k-1)*2*2 - 6*BBE_SIMD_WIDTH);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-((N-k-1)*2),(uintptr_t)pRrd);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx0,rr0);
            BBE_MULSNX16C(B,xx1,rr1);
            BBE_MULSNX16C(B,xx2,rr2);
            BBE_MULSNX16C(B,xx3,rr3);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
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
            BBE_SPNX16_XP(b_res,pXwr,2*SX);
        }
    }
    #undef SR
    #undef SX
    #undef SD
    #undef N
} /* cqr_bkw32x32x1n() */
#endif

size_t cqr_bkw32x32x1n_getScratchSize (int M, int N,int P,int L)
{
    (void)M;
    (void)N;
    (void)P;
    (void)L;
    return 0;
} /* cqr_bkw32x32x1n_getScratchSize() */

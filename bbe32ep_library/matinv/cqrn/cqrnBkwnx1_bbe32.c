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
    BBE32 code for QR decomposition (bkw part), block format
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

#if HAVE_CQRN

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

/*-------------------------------------------------------
    backward recursion: P==1

   Input:
    M, N, P      dimensional parameters
    L            number of matrices
    qABX         qA-qB+qX
   Input/output:
    x[L][SB][2]  at the input it is the sequence of L updated right parts Z=Q'B.
                 They will be replaced with MMSE solution vectors X (only N*P 
                 elements are used)
   Input:
    R[L][SA][2]  Upper triangle matrices R (only N*N 
                 elements of each matrix are used)
    D[L][SD][2]  reciprocal of main diagonal (mantissa, exponent) 
                 in the special format
-------------------------------------------------------*/
/*
    backward recursion: P==1
*/
void cqrnBkwnx1(int16_t* restrict X, 
          const int16_t* restrict R,
          const int16_t* restrict D,
                int qABX,
                int M,int N,int L)
{
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    const int SR=2*getSpace(M*N);
    const int SX=2*getSpace(M);
    const int SD=2*getSpace(N);

    int k;
    int l;
    int m;
    
    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;

    const xb_vecNx16* restrict pXrd_first;
    const xb_vecNx16* restrict pXrd;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    xb_vecNx16 xx;
    xb_vecNx16 rr;

    valign Xrd_align, Rrd_align;

    const vsaN sh16= BBE_MOVVSA32(16);
    const vsaN vqABX= BBE_MOVVSA32(qABX);

    __Pragma("loop_count min=1")
    for (k=N-1; k>=0; k--)
    {
        pXrd_first= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)X);
        pXrd= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXrd_first, 2*2);
        pXwr= (xb_vecNx16*)pXrd_first;
        pRrd= (const xb_vecNx16*)(R+(k*N+(k+1))*2);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);

        __Pragma("loop_count min=1")
        for(l=0; l<L; l++)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP
            // load 16 bit complex value (16 bit re and 16 bit im)
            BBE_LPNX16_XP(b,pXrd_first,2*SX);
            B= BBE_UNPKSNX16(b);
            B=BBE_SLSNX40(B,vqABX);
            Rrd_align= BBE_LA_PP(pRrd);
            Xrd_align= BBE_LA_PP(pXrd);
            for (m= N-k-1; m>0; m-= BBE_SIMD_WIDTH/2)
            {
                BBE_LAVNX16_XP(rr,Rrd_align,pRrd,m*4);
                BBE_LAVNX16_XP(xx,Xrd_align,pXrd,m*4);
                //B-= rr*xx
                BBE_MULSNX16C(B,xx,rr);
            }
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            pXrd= (const xb_vecNx16*)XT_ADDX2(SX - (N-k-1)*2,(uintptr_t)pXrd);
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
            BBE_SPNX16_XP(b_res,pXwr,2*SX);
        }
    }
}
#endif

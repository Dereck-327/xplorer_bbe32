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
  NatureDSP_Baseband library. Cholesky forward recursion for block ordered matrices:
    These functions make forward recursion stage of pseudo-inversion. They use
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
#include "choln_common.h"

#if HAVE_CHOLN

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


static void computeABmxnx1(int32_t * Z,
                     const int16_t * A,
                     const int16_t * B,
                     const int qYB,
                     const int M,
                     const int N,
                     const int L)
{
    int l, n, m;
    int delta= 0;
    xb_vecNx16 a0,a1, b; 
    xb_vecNx16 zl, zh;
    xb_vecNx40 acc0,acc1;
    const vsaN vqYB= BBE_MOVVSA32(qYB);

    valign A_align, Z_align;

    const xb_vecNx16* restrict pArd;
    const xb_vecNx16* restrict pBrd;
          xb_vecNx16* restrict pZwr= (xb_vecNx16*)Z;

    const xb_vecNx16* restrict pArd_this_row;
    const xb_vecNx16* restrict pArd_this_matrix= (const xb_vecNx16*)A;
    const xb_vecNx16* restrict pBrd_this_matrix= (const xb_vecNx16*)B;

    const int SA= 2*getSpace(M*N);
    const int SB= 2*getSpace(M);
    
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=1")
    for(l= 0; l<L; l++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for (n= N; n > 0; n-= 16)
        {
            // set pArd_this_row to (N-n)-th col
            pArd_this_row= (const xb_vecNx16*)XT_ADD(delta,(uintptr_t)pArd_this_matrix);
            // set pBrd to pBrd_this_matrix
            pBrd = pBrd_this_matrix;
            // load B[l][0]
            BBE_LPNX16_IP(b,pBrd,2*2);
            b= BBE_REPNX16C(b,0);

            // load n (max 16) elements of A[l][0] row
            pArd= pArd_this_row;
            A_align= BBE_LA_PP(pArd);
            BBE_LAVNX16_XP(a0,A_align,pArd,n*4);
            BBE_LAVNX16_XP(a1,A_align,pArd,n*4-2*BBE_SIMD_WIDTH);
            // multiplication
            acc0= BBE_MULNX16J(b, a0);
            acc1= BBE_MULNX16J(b, a1);
            // set pArd_this_row to next row
            //pArd_this_row=(const xb_vecNx16*)XT_ADD(N*2*2,(uintptr_t)pArd_this_row);
            pArd_this_row=(const xb_vecNx16*)XT_ADDX4(N,(uintptr_t)pArd_this_row);
            __Pragma("loop_count min=3")
            for (m=1; m < M; m++)
            {
                // load B[l][m]
                BBE_LPNX16_IP(b,pBrd,2*2);
                b= BBE_REPNX16C(b,0);

                // load n (max 16) elements of A[l][m] row
                pArd= pArd_this_row;
                A_align= BBE_LA_PP(pArd);
                BBE_LAVNX16_XP(a0,A_align,pArd,n*4);
                BBE_LAVNX16_XP(a1,A_align,pArd,n*4-2*BBE_SIMD_WIDTH);

                BBE_MULANX16J(acc0, b, a0);
                BBE_MULANX16J(acc1, b, a1);
                // set pArd_tis_row to next row
                //pArd_this_row= (const xb_vecNx16*)XT_ADD(N*2*2,(uintptr_t)pArd_this_row);
                pArd_this_row=(const xb_vecNx16*)XT_ADDX4(N,(uintptr_t)pArd_this_row);
            }
            // shift acc
            acc0= BBE_SLSNX40(acc0,vqYB);
            zl= BBE_MOVSVWL(acc0);
            zh= BBE_MOVSVWH(acc0);
            BBE_SANX16_IP(zl,Z_align,pZwr);
            BBE_SAVNX16_XP(zh,Z_align,pZwr,n*2*4 - 2*BBE_SIMD_WIDTH);
            acc1= BBE_SLSNX40(acc1,vqYB);
            zl= BBE_MOVSVWL(acc1);
            zh= BBE_MOVSVWH(acc1);
            BBE_SAVNX16_XP(zl,Z_align,pZwr,n*2*4 - 4*BBE_SIMD_WIDTH);
            BBE_SAVNX16_XP(zh,Z_align,pZwr,n*2*4 - 6*BBE_SIMD_WIDTH);
            delta= XT_ADD(16*2*2,delta);

        }
        // to next matrix
        pArd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SA,(uintptr_t)pArd_this_matrix);
        pBrd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SB,(uintptr_t)pBrd_this_matrix);
    }
    BBE_SAPOS_FP(Z_align,pZwr);
}

static void computeABmxnxp(int32_t* Z,const int16_t* A,const int16_t* B, int qYB, int M,int N,int P,int L)
#if 0
{
    #define VECLEN ((int)(sizeof(xb_vecNx16)/(2*sizeof(int16_t))))
    int l, n, m, p;
    int delta= 0;
    int delta2= 0;
    xb_vecNx16 a, b; 
    xb_vecNx16 zl, zh;
    xb_vecNx40 acc;
    const vsaN vqYB=BBE_MOVVSA32(qYB);

    valign B_align, Z_align;

    const xb_vecNx16* restrict pArd;
    const xb_vecNx16* restrict pBrd;
          xb_vecNx16* restrict pZwr;

    //const xb_vecNx16* restrict pArd_this_row;
    const xb_vecNx16* restrict pArd_this_matrix= (const xb_vecNx16*)A;
    const xb_vecNx16* restrict pBrd_this_col;
    const xb_vecNx16* restrict pBrd_this_row;
    const xb_vecNx16* restrict pBrd_this_matrix= (const xb_vecNx16*)B;
    xb_vecNx16* restrict pZwr_this_matrix= (      xb_vecNx16*)Z;
    xb_vecNx16* restrict pZwr_this_row;

    const int SA= 2*getSpace(M*N);
    const int SB= 2*getSpace(M*P);
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=1")
    for(l=0; l<L; l++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for(p=P; p > 0 ;p-= VECLEN)
        {
            // set pZwr to p-th col of Z
            pZwr= (xb_vecNx16*)XT_ADDX4(delta,(uintptr_t)pZwr_this_matrix);
            // set pBrd to p-th col of B
            pBrd_this_col= (xb_vecNx16*)XT_ADDX2(delta,(uintptr_t)pBrd_this_matrix);
            delta2= XT_XOR(delta2, delta2);
            __Pragma("loop_count min=4")
            for (n=0; n<N; n++)
            {
                 // set pArd_this_row to n-th col
                pArd= (const xb_vecNx16*)XT_ADDX2(delta2,(uintptr_t)pArd_this_matrix);
                // set pBrd to pBrd_this_col
                pBrd_this_row= pBrd_this_col;
                // load p elements from B[l][0][p]
                pBrd= pBrd_this_row; 
                B_align= BBE_LA_PP(pBrd);
                BBE_LAVNX16_XP(b,B_align,pBrd,p*4);
                // set pBrd_this_row to next row
                pBrd_this_row=(const xb_vecNx16*)XT_ADD(P*2*2,(uintptr_t)pBrd_this_row);
                //pBrd_this_row=(const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pBrd_this_row);
                // load A[l][0][n] and set pArd to next row
                BBE_LPNX16_XP(a,pArd,2*N*2);
                a= BBE_REPNX16C(a,0);
                // multiplication
                acc= BBE_MULNX16J(b, a);
                __Pragma("loop_count min=3")
                for (m= 1; m < M; m++)
                {
                    // load p elements from B[l][0][p]
                    pBrd= pBrd_this_row; 
                    B_align= BBE_LA_PP(pBrd);
                    BBE_LAVNX16_XP(b,B_align,pBrd,p*4);
                    // set pBrd_this_row to next row
                    pBrd_this_row=(const xb_vecNx16*)XT_ADD(P*2*2,(uintptr_t)pBrd_this_row);
                    //pBrd_this_row=(const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pBrd_this_row);
                    BBE_LPNX16_XP(a,pArd,2*N*2);
                    a= BBE_REPNX16C(a,0);
                    BBE_MULANX16J(acc, b, a);
                }
                // shift acc
                acc= BBE_SLSNX40(acc,vqYB);
                // moves 8 40-bit elements from the lower half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
                zl= BBE_MOVSVWL(acc);
                zh= BBE_MOVSVWH(acc);

                pZwr_this_row= pZwr;
                BBE_SAVNX16_XP(zl,Z_align,pZwr,p*2*4);
                BBE_SAVNX16_XP(zh,Z_align,pZwr,p*2*4 - 2*BBE_SIMD_WIDTH);
                BBE_SAPOS_FP(Z_align,pZwr);
                // set pZwr to next row
                pZwr= (xb_vecNx16*)XT_ADDX4(P*2,(uintptr_t)pZwr_this_row);
                delta2= XT_ADD(2, delta2);
            }
            delta= XT_ADDI(delta, 2*VECLEN);
        }
        // to next matrix
        pArd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SA,(uintptr_t)pArd_this_matrix);
        pBrd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SB,(uintptr_t)pBrd_this_matrix);
        pZwr_this_matrix= (xb_vecNx16*)XT_ADDX4(N*P*2,(uintptr_t)pZwr_this_matrix);   
    }
    #undef VECLEN
}
#else
{
    int l, n, m, p;
    int delta= 0;
    int delta2= 0;
    xb_vecNx16 a, b; 
    xb_vecNx16 zl, zh;
    xb_vecNx40 acc;
    const vsaN vqYB=BBE_MOVVSA32(qYB);

    valign B_align, Z_align;

    const xb_vecNx16* restrict pArd;
    const xb_vecNx16* restrict pBrd;
          xb_vecNx16* restrict pZwr;

    //const xb_vecNx16* restrict pArd_this_row;
    const xb_vecNx16* restrict pArd_this_matrix= (const xb_vecNx16*)A;
    const xb_vecNx16* restrict pBrd_this_col;
    const xb_vecNx16* restrict pBrd_this_row;
    const xb_vecNx16* restrict pBrd_this_matrix= (const xb_vecNx16*)B;
    xb_vecNx16* restrict pZwr_this_matrix= (      xb_vecNx16*)Z;
    xb_vecNx16* restrict pZwr_this_row;

    const int SA= 2*getSpace(M*N);
    const int SB= 2*getSpace(M*P);
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=1")
    for(l=0; l<L; l++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for(p=P; p > 0 ;p-= 8)
        {
            // set pZwr to p-th col of Z
            pZwr= (xb_vecNx16*)XT_ADDX4(delta,(uintptr_t)pZwr_this_matrix);
            // set pBrd to p-th col of B
            pBrd_this_col= (xb_vecNx16*)XT_ADDX2(delta,(uintptr_t)pBrd_this_matrix);
            delta2= XT_XOR(delta2, delta2);
            __Pragma("loop_count min=4")
            for (n=0; n<N; n++)
            {
                 // set pArd_this_row to n-th col
                pArd= (const xb_vecNx16*)XT_ADDX2(delta2,(uintptr_t)pArd_this_matrix);
                // set pBrd to pBrd_this_col
                pBrd_this_row= pBrd_this_col;
                // load p elements from B[l][0][p]
                pBrd= pBrd_this_row; 
                B_align= BBE_LA_PP(pBrd);
                BBE_LAVNX16_XP(b,B_align,pBrd,p*4);
                // set pBrd_this_row to next row
                pBrd_this_row=(const xb_vecNx16*)XT_ADD(P*2*2,(uintptr_t)pBrd_this_row);
                //pBrd_this_row=(const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pBrd_this_row);
                // load A[l][0][n] and set pArd to next row
                BBE_LPNX16_XP(a,pArd,2*N*2);
                a= BBE_REPNX16C(a,0);
                // multiplication
                acc= BBE_MULNX16J(b, a);
                __Pragma("loop_count min=3")
                for (m= 1; m < M; m++)
                {
                    // load p elements from B[l][0][p]
                    pBrd= pBrd_this_row; 
                    B_align= BBE_LA_PP(pBrd);
                    BBE_LAVNX16_XP(b,B_align,pBrd,p*4);
                    // set pBrd_this_row to next row
                    pBrd_this_row=(const xb_vecNx16*)XT_ADD(P*2*2,(uintptr_t)pBrd_this_row);
                    //pBrd_this_row=(const xb_vecNx16*)XT_ADDX4(P,(uintptr_t)pBrd_this_row);
                    BBE_LPNX16_XP(a,pArd,2*N*2);
                    a= BBE_REPNX16C(a,0);
                    BBE_MULANX16J(acc, b, a);
                }
                // shift acc
                acc= BBE_SLSNX40(acc,vqYB);
                // moves 8 40-bit elements from the lower half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
                zl= BBE_MOVVWL(acc);
                zh= BBE_MOVVWH(acc);

                pZwr_this_row= pZwr;
                BBE_SAVNX16_XP(zl,Z_align,pZwr,p*2*4);
                BBE_SAVNX16_XP(zh,Z_align,pZwr,p*2*4 - 2*BBE_SIMD_WIDTH);
                BBE_SAPOS_FP(Z_align,pZwr);
                // set pZwr to next row
                pZwr= (xb_vecNx16*)XT_ADDX4(P*2,(uintptr_t)pZwr_this_row);
                delta2= XT_ADD(2, delta2);
            }
            delta= XT_ADDI(delta, 2*8);
        }
        // to next matrix
        pArd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SA,(uintptr_t)pArd_this_matrix);
        pBrd_this_matrix= (const xb_vecNx16*)XT_ADDX2(SB,(uintptr_t)pBrd_this_matrix);
        pZwr_this_matrix= (xb_vecNx16*)XT_ADDX4(N*P*2,(uintptr_t)pZwr_this_matrix);   
    }
    BBE_SVNX16_I(BBE_ZERONX16(), pZwr_this_matrix, 0);
}
#endif


/*
    another algorithm for P!=1
*/
static void fwdnxp( int16_t* restrict y, 
              const int16_t* restrict R, 
              const int16_t* restrict D, 
              const int32_t* restrict Z, 
                    int N,int P,int L)
#if 0
{
    int SR=2*getSpace((N*(N+1))>>1);
    int SY=2*getSpace(N*P);
    int SD=2*getSpace(N);
    int l,n,m,p;
    const int16_t *pR;
    int64_t B_re,B_im;
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z,2*BBE_SIMD_WIDTH);
    NASSERT(P>1 && N>0 );

    for(n=0; n<N; n++)
    {
        for(p=0; p<P; p++)
        {
            for(l=0; l<L; l++)
            {
                pR=R+(n*(n+1))+l*SR;
                B_re=Z[2*l*N*P+2*n*P+2*p+0];
                B_im=Z[2*l*N*P+2*n*P+2*p+1];
                for (m=0; m<n; m++)   
                {
                    int16_t r_re,r_im;
                    int16_t y_re,y_im;
                    r_re=pR[2*m+0];
                    r_im=pR[2*m+1];
                    y_re=y[l*SY+m*P*2+p*2+0];
                    y_im=y[l*SY+m*P*2+p*2+1];
                    B_re-=L_mul_ss(y_re,r_re)+L_mul_ss(y_im,r_im); // representation qA+qY
                    B_im-=L_mul_ss(y_im,r_re)-L_mul_ss(y_re,r_im);
                }
                B_re=mul32x16su((int32_t)B_re,D[l*SD+2*n+0]); // representation (qA+qY),qD -> qA+qY+qD-31
                B_im=mul32x16su((int32_t)B_im,D[l*SD+2*n+0]); // 
                y[l*SY+m*P*2+p*2+0]=packr1x40(B_re,D[l*SD+2*n+1]);
                y[l*SY+m*P*2+p*2+1]=packr1x40(B_im,D[l*SD+2*n+1]);
            }
        }
    }
}
#elif 1
{
    int l, p;
    int n, m;
    int delta_pY= 0;
    int delta= 0;
    const int32_t* restrict pZrd;
    const xb_vecNx16* restrict pYrd;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;
    int16_t* restrict pYwr;

    xb_vecNx16* restrict pYwr_tmp;
    const xb_vecNx16* restrict pYrd_tmp;
    const xb_vecNx16* restrict pZrd_tmp;

    xb_vecNx16 dd,d0,rr,yy,tl,th,b_res, bl, bh;
    xb_vecNx40 B;
    vsaN q;
    valign B_align, Y_align, Z_align;

    const vsaN sh16= BBE_MOVVSA32(16);

    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));
    NASSERT(P>1 && N>0 );

    const int SR= 2*getSpace((N*(N+1))>>1);
    const int SY= 2*getSpace(N*P);
    const int SD= 2*getSpace(N);
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=4")
    for(n=0; n<N; n++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for(p= P; p>0; p-=8)
        {
            //pDrd = (const xb_vecNx16*)&D[n+n];
            pDrd= (const xb_vecNx16*)XT_ADDX4(n, (uintptr_t)D);
            pRrd= (const xb_vecNx16*)(R+(n*(n+1)));
            pYrd= (const xb_vecNx16*)XT_ADDX2(delta, (uintptr_t)y);
            pYwr= (int16_t*)XT_ADDX2(delta, (uintptr_t)&y[2*n*P]);
            pZrd= (int32_t*)XT_ADDX4(delta, (uintptr_t)&Z[2*n*P]);
            delta= XT_ADDI(delta,2*8);
            __Pragma("loop_count min=1")
            for(l=0; l<L; l++)
            {
                // load up to 4 32 complex elements
                pZrd_tmp =(const xb_vecNx16*)pZrd;
                B_align= BBE_LA_PP((const xb_vecNx16 *)pZrd_tmp);
                BBE_LAVNX16_XP(bl,B_align,pZrd_tmp,8*p);
                BBE_LAVNX16_XP(bh,B_align,pZrd_tmp,8*p - 2*BBE_SIMD_WIDTH);
                B= BBE_MOVSWV(bh,bl);
                // pZrd+=P*N*2
                pZrd= (const int32_t*)XT_ADDX4(P*N*2,(uintptr_t)pZrd);
                // load rr and yy (load only n elements. others = 0)
                delta_pY= XT_XOR(delta_pY, delta_pY);
                __Pragma("no_unroll");
                for (m= 0; m < n; m++)
                {
                    BBE_LPNX16_IP(rr, pRrd, 2*2);
                    rr= BBE_REPNX16C(rr,0);
                    pYrd_tmp= (const xb_vecNx16*)XT_ADD(delta_pY, (uintptr_t)pYrd);
                    Y_align= BBE_LA_PP(pYrd_tmp);
                    BBE_LAVNX16_XP(yy,Y_align,pYrd_tmp,4*p);
                    //B-= rr*yy
                    BBE_MULSNX16J(B,yy,rr);
                    delta_pY= XT_ADDX4(P, delta_pY);
                }
                pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
                pYrd=(const xb_vecNx16*)XT_ADDX2(SY,(uintptr_t)pYrd);
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
                pYwr_tmp= (xb_vecNx16*)pYwr;
                BBE_SAVNX16_XP(b_res,Z_align,pYwr_tmp,4*p);
                BBE_SAPOS_FP(Z_align,pYwr_tmp);
                pYwr= (int16_t*)XT_ADDX2(SY,(uintptr_t)pYwr);
            }
        }
    }
}
#else
{
    int l, p;
    int n, m;
    int delta_pY= 0;
    int delta= 0;
    const int32_t* restrict pZrd;
    const xb_vecNx16* restrict pYrd;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;
    int16_t* restrict pYwr;

    xb_vecNx16* restrict pYwr_tmp;
    const xb_vecNx16* restrict pYrd_tmp;
    const xb_vecNx16* restrict pZrd_tmp;

    xb_vecNx16 dd,d0,rr,yy0,yy1,tl,th,b_res, bl, bh;
    xb_vecNx40 B0,B1;
    vsaN q;
    valign B_align, Y_align, Z_align;

    const vsaN sh16= BBE_MOVVSA32(16);

    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));
    NASSERT(P>1 && N>0 );

    const int SR= 2*getSpace((N*(N+1))>>1);
    const int SY= 2*getSpace(N*P);
    const int SD= 2*getSpace(N);
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=4")
    for(n=0; n<N; n++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for(p= P; p>0; p-=8)
        {
            //pDrd = (const xb_vecNx16*)&D[n+n];
            pDrd= (const xb_vecNx16*)XT_ADDX4(n, (uintptr_t)D);
            pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
            pYrd= (const xb_vecNx16*)XT_ADDX2(delta, (uintptr_t)y);
            pYwr= (int16_t*)XT_ADDX2(delta, (uintptr_t)&y[2*n*P]);
            pZrd= (int32_t*)XT_ADDX4(delta, (uintptr_t)&Z[2*n*P]);
            delta= XT_ADDI(delta,2*8);
            __Pragma("loop_count min=1")
            for(l=0; l<L; l++)
            {
                // load up to 4 32 complex elements
                pZrd_tmp =(const xb_vecNx16*)pZrd;
                B_align= BBE_LA_PP((const xb_vecNx16 *)pZrd_tmp);
                BBE_LAVNX16_XP(bl,B_align,pZrd_tmp,8*p);
                BBE_LAVNX16_XP(bh,B_align,pZrd_tmp,8*p - 2*BBE_SIMD_WIDTH);
                B0= BBE_MOVSWV(bh,bl);
                BBE_LAVNX16_XP(bl,B_align,pZrd_tmp,8*p - 4*BBE_SIMD_WIDTH);
                BBE_LAVNX16_XP(bh,B_align,pZrd_tmp,8*p - 6*BBE_SIMD_WIDTH);
                B1= BBE_MOVSWV(bh,bl);
                // pZrd+=P*N*2
                pZrd= (const int32_t*)XT_ADDX4(P*N*2,(uintptr_t)pZrd);
                // load rr and yy (load only n elements. others = 0)
                delta_pY= XT_XOR(delta_pY, delta_pY);
                //__Pragma("loop_count min=1")
                for (m= 0; m < n; m++)
                {
                    BBE_LPNX16_IP(rr, pRrd, 2*2);
                    rr= BBE_REPNX16C(rr,0);
                    pYrd_tmp= (const xb_vecNx16*)XT_ADD(delta_pY, (uintptr_t)pYrd);
                    Y_align= BBE_LA_PP(pYrd_tmp);
                    BBE_LAVNX16_XP(yy0,Y_align,pYrd_tmp,4*p);
                    BBE_LAVNX16_XP(yy1,Y_align,pYrd_tmp,4*p-2*BBE_SIMD_WIDTH);
                    //B-= rr*yy
                    BBE_MULSNX16J(B0,yy0,rr);
                    BBE_MULSNX16J(B1,yy1,rr);
                    delta_pY= XT_ADD(2*P*2, delta_pY);
                }
                pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
                pYrd=(const xb_vecNx16*)XT_ADDX2(SY,(uintptr_t)pYrd);
                // load D
                BBE_LPNX16_XP(dd,pDrd,2*SD);
                // replicate dd[0]
                d0 = BBE_REPNX16(dd,0);
                // 32x16 multiply
                tl = BBE_PACKLNX40(B0);
                th = BBE_PACKVNX40(B0,sh16);
                B0 = BBE_MULUUNX16(d0,tl);
                B0 = BBE_SRAINX40(B0,16);
                BBE_MULUSANX16(B0,d0,th);
                tl = BBE_PACKLNX40(B1);
                th = BBE_PACKVNX40(B1,sh16);
                B1 = BBE_MULUUNX16(d0,tl);
                B1 = BBE_SRAINX40(B1,16);
                BBE_MULUSANX16(B1,d0,th);

                // make q
                d0= BBE_REPNX16(dd,1);
                q = BBE_MOVVSV(d0,0);
                // result rounding
                B0 = BBE_RNDADJNX40(B0,q);
                B1 = BBE_RNDADJNX40(B1,q);
                // store res
                pYwr_tmp= (xb_vecNx16*)pYwr;
                b_res = BBE_PACKVNX40(B0,q);
                BBE_SAVNX16_XP(b_res,Z_align,pYwr_tmp,4*p);
                b_res = BBE_PACKVNX40(B1,q);
                BBE_SAVNX16_XP(b_res,Z_align,pYwr_tmp,4*p-2*BBE_SIMD_WIDTH);
                BBE_SAPOS_FP(Z_align,pYwr_tmp);
                pYwr= (int16_t*)XT_ADDX2(SY,(uintptr_t)pYwr);
            }
        }
    }
}
#endif

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Matrix sizes SA,SR,SD,SB,SY are selected as usual for complex block ordered 
matrix sequencies, i.e. total size is rounded up to the closest bigger 
multiple of BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the 
closest bigger multiple of degree of 2. 
SA=size(M*N)
SR=size(((N+1)*N)/2)
SD=size(N)
SB=size(M*P)
SY=size(N*P)
Scratch size in bytes is defined by cholfwdxxx_getScratchSize()

Input:
 M            Matrix dimension (number of rows in matrices A)
 N            Matrix dimension (number of columns and rows in 
              matrices R)
 P            Number of columns in right-side matrices B
 L            Number of matrices
 R[L][SR]     Sequence of L upper triangular complex matrices R
 A[L][SA]     Sequence of L complex matrices A
 B[L][SB]     Sequence original right-side matrices B
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format
qB,qY         Fixed point representation of matrices B and y
Output:
y[L][SY]      Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
3. Matrix sizes M,N,P must be positive
4. M and N must be multiples of 4 
5. M>=N
---------------------------------------------------------------------------*/

void cholfwdmxnxpn ( void    * pScr,
                     complex_fract16 * restrict _y,
               const complex_fract16 * restrict _R, 
               const complex_fract16 * restrict _D,
               const complex_fract16 * restrict _A, 
               const complex_fract16 * restrict _B, 
                     int qB, int qY,
                     int M , int N ,
                     int P , int L )
{
          int16_t * restrict y=(      int16_t *)_y;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict A=(const int16_t *)_A;
    const int16_t * restrict B=(const int16_t *)_B;
    int32_t* Z=(int32_t* )pScr;
    int qYB=qY-qB;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT(M%4==0 && N%4==0 && M>0 && N>0 && L>0 && M>=N);
    if (L<=0 || N<=0 || P<=0)  return;

    // compute A'*B
    // computeAB(Z,A,B, qYB, M,N,P,L);

    if (P==1)
    {
        int SR=2*getSpace((N*(N+1))>>1);
        int SY=2*getSpace(N);
        int SD=2*getSpace(N);
        computeABmxnx1(Z,A,B, qYB, M,N,L);
        cholnFwdrec(y,R,D,Z,N,L,SR,SD,SY,2*N);
    }
    else
    {
        computeABmxnxp(Z,A,B, qYB, M,N,P,L);
        fwdnxp(y,R,D,Z,N,P,L);
    }
} /* cholfwdmxnxpn() */

/* scratch allocation function */
size_t cholfwdmxnxpn_getScratchSize (int M,int N, int P,int L)
{
    size_t Z_size;
    NASSERT(M%4==0 && N%4==0 && M>0 && N>0 && L>0 && M>=N);
    M=XT_MAX(0,M);
    N=XT_MAX(0,N);
    P=XT_MAX(0,P);
    L=XT_MAX(0,L);
    Z_size= (2*L*N*P*sizeof(int32_t));
    return (Z_size + 2 * BBE_SIMD_WIDTH);
} /* cholfwdmxnxpn_getScratchSize() */

#else
DISCARD_FUN(void,cholfwdmxnxpn,( void    * pScr,
                                 complex_fract16 * restrict y,
                           const complex_fract16 * restrict R,
                           const complex_fract16 * restrict D,
                           const complex_fract16 * restrict A,
                           const complex_fract16 * restrict B,
                                 int qB, int qY,
                                 int M , int N ,
                                 int P , int L ) )
size_t cholfwdmxnxpn_getScratchSize  (int M,int N, int P,int L) { (void)M,(void)N,(void)P,(void)L; return 0; }
#endif

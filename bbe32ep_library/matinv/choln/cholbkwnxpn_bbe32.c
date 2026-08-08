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
  Cholesky backward recursion for block ordered matrices
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

/*
    backward recursion: P!=1
*/
static void cholnBkwnxp( int16_t * restrict x, 
                   const int16_t * restrict Rt,
                   const int16_t * restrict D,
                   const int16_t * restrict y, 
                         int qXYA,
                         int N, int P, int L)
#if 0
{
    int m,k,p;
    int l;
    int SX=2*getSpace(N*P);
    int SD=2*getSpace(N);

    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    for (k=N-1; k>=0; k--)
    {
        for(l=0; l<L; l++)
        {
            const int16_t* pRt=Rt+l*N*(N-1)+(N-k-1)*(N-k-2);
              // calculate y(m,:)-R(m,:)*X, 1xP
            for(p=0; p<P; p++)
            {
                int16_t r_re,r_im;
                int16_t x_re,x_im;
                int32_t A_re,A_im;
                int64_t B_re,B_im;
                A_re=(y[l*SX+2*k*P+p*2+0]);
                A_im=(y[l*SX+2*k*P+p*2+1]);
                B_re=L_shl_l(A_re,qXYA); // qY->qX+qR
                B_im=L_shl_l(A_im,qXYA);  

                for (m=0; m<N-k-1; m++)
                {
                    x_re=x[l*SX+2*(k+1+m)*P+2*p+0];
                    x_im=x[l*SX+2*(k+1+m)*P+2*p+1];
                    r_re=pRt[2*m+0];
                    r_im=pRt[2*m+1];
                    B_re-=L_mul_ss(x_re, r_re)-L_mul_ss(x_im, r_im);    // ->qX+qR
                    B_im-=L_mul_ss(x_re, r_im)+L_mul_ss(x_im, r_re); 
                }
                // NOTE: having this 32x16 multiple is critical !
                A_re=mul32x16su((int32_t)B_re,D[l*SD+2*k+0]); // representation (qA+qY),qD -> qA+qY+qD-31
                A_im=mul32x16su((int32_t)B_im,D[l*SD+2*k+0]); // 
                x[l*SX+2*k*P+2*p+0]=packr1x40(A_re,D[l*SD+2*k+1]);
                x[l*SX+2*k*P+2*p+1]=packr1x40(A_im,D[l*SD+2*k+1]);
            }
        }
    }
}
#else
{
    int k;
    int l;
    int m;
    int p;
    int delta= 0;
    xb_vecNx16 tl,th,dd,d0,d1,b0,b1;

    xb_vecNx40 B0,B1;
    vsaN q;
    valign X_align, Y_align, Z_align;

    const int16_t   * restrict pYrd;
    const xb_vecNx16* restrict pXrd;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    const xb_vecNx16* restrict pXrd_tmp;
    const xb_vecNx16* restrict pYrd_tmp;
          xb_vecNx16* restrict pXwr_tmp;

    xb_vecNx16 coef;
    xb_vecNx16 xx0,xx1;
    xb_vecNx16 rr;

    const int SX= 2*getSpace(N*P);
    const int SD= 2*getSpace(N);

    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));

    // load shift amount
    coef = 1;
    coef = BBE_SLSNX16(coef,qXYA);

    const vsaN sh16= BBE_MOVVSA32(16);
    Z_align= BBE_ZALIGN();

    if (P<=8)
    {
        P = P<<2;
        __Pragma("loop_count min=1");
        for (k=N-1; k>=0; k--)
        {
            pXrd = (const xb_vecNx16*)((int8_t*)x+k*P+P);
            pXwr = (      xb_vecNx16*)((int8_t*)x+k*P);
            pYrd = (const int16_t   *)((int8_t*)y+k*P);
            //pDrd = (const xb_vecNx16*)&D[k+k];
            pDrd = (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
            pRrd = (const xb_vecNx16*)(Rt+(N-k-1)*(N-k-2));
            __Pragma("loop_count min=1");
            for(l=0; l<L; l++)
            {
                // calculate y(m,:)-R(m,:)*X, 1xP
                // load 16 bit complex value (16 bit re and 16 bit im)
                // load B and shl with saturation
                pXrd_tmp= pXrd;
                pYrd_tmp= (const xb_vecNx16*)pYrd;
                pYrd= (const int16_t*)XT_ADDX2(SX,(uintptr_t)pYrd);
                Y_align= BBE_LA_PP(pYrd_tmp);
                X_align= BBE_LA_PP(pXrd_tmp);
                BBE_LAVNX16_XP(b0,Y_align,pYrd_tmp,P);
                B0= BBE_MULNX16(b0,coef);
                // load D
                BBE_LPNX16_XP(dd,pDrd,2*SD);
                // replicate dd[0]
                d0 = BBE_REPNX16(dd,0);
                // make q
                d1= BBE_REPNX16(dd,1);
                q = BBE_MOVVSV(d1,0);
                __Pragma("no_unroll");
                for (m=0; m<N-1-k; m++)
                {
                    BBE_LAVNX16_XP(xx0,X_align,pXrd_tmp,P);
                    BBE_LPNX16_IP(rr, pRrd, 2*2);

                    //B-= rr*xx
                    rr= BBE_REPNX16C(rr,0);
                    BBE_MULSNX16C(B0,xx0,rr);
                }
                // pRrd++
                pRrd=(const xb_vecNx16*)XT_ADDX2(N*(N-1)-2*(N-k-1),(uintptr_t)pRrd);
                // pXrd++
                pXrd=(const xb_vecNx16*)XT_ADDX2(SX,(uintptr_t)pXrd);
                // 32x16 multiply
                tl = BBE_PACKLNX40(B0);
                th = BBE_PACKVNX40(B0,sh16);
                B0 = BBE_MULUUNX16(d0,tl);
                B0 = BBE_SRAINX40(B0,16);

                // result rounding
                B1= BBE_MULUSRNX16(d0,th,q);
                B0 = BBE_ADDNX40(B0,B1);
                // store res
                b0 = BBE_PACKVNX40(B0,q);
                pXwr_tmp= pXwr;
                BBE_SAVNX16_XP(b0,Z_align,pXwr_tmp,P);
                BBE_SAPOS_FP(Z_align,pXwr_tmp);
                pXwr= (xb_vecNx16*)XT_ADDX2(SX, (uintptr_t)pXwr);
            }
        }
    }
    else
    {
        __Pragma("loop_count min=1");
        for (k=N-1; k>=0; k--)
        {
            delta= XT_XOR(delta, delta);
            __Pragma("loop_count min=1");
            for(p= P; p>0; p-=16)
            {
                pXrd = (const xb_vecNx16*)XT_ADDX2(delta,(uintptr_t)&x[2*k*P+2*P]);
                pXwr = (xb_vecNx16*)XT_ADDX2(delta,(uintptr_t)&x[2*k*P]);
                pYrd = (const int16_t*)XT_ADDX2(delta,(uintptr_t)&y[2*k*P]);
                //pDrd = (const xb_vecNx16*)&D[k+k];
                pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
                pRrd = (const xb_vecNx16*)(Rt+(N-k-1)*(N-k-2));
                delta= XT_ADDI(delta, 2*16);
                __Pragma("loop_count min=1")
                for(l=0; l<L; l++)
                {
                    // calculate y(m,:)-R(m,:)*X, 1xP
                    // load 16 bit complex value (16 bit re and 16 bit im)
                    // load B and shl with saturation
                    pYrd_tmp= (const xb_vecNx16*)pYrd;
                    pYrd= (const int16_t*)XT_ADDX2(SX,(uintptr_t)pYrd);
                    Y_align= BBE_LA_PP(pYrd_tmp);
                    BBE_LAVNX16_XP(b0,Y_align,pYrd_tmp,p*4);
                    BBE_LAVNX16_XP(b1,Y_align,pYrd_tmp,p*4-2*BBE_SIMD_WIDTH);
                    B0= BBE_MULNX16(b0,coef);
                    B1= BBE_MULNX16(b1,coef);
                    // load D
                    BBE_LPNX16_XP(dd,pDrd,2*SD);
                    // replicate dd[0]
                    d0 = BBE_REPNX16(dd,0);
                    // make q
                    d1= BBE_REPNX16(dd,1);
                    q = BBE_MOVVSV(d1,0);
                    for (m=0; m<N-1-k; m++)
                    {
                        pXrd_tmp= pXrd;
                        X_align= BBE_LA_PP(pXrd_tmp);
                        BBE_LAVNX16_XP(xx0,X_align,pXrd_tmp,p*4);
                        BBE_LAVNX16_XP(xx1,X_align,pXrd_tmp,p*4-2*BBE_SIMD_WIDTH);
                        pXrd= (const xb_vecNx16*)XT_ADDX4(P, (uintptr_t)pXrd);

                        BBE_LPNX16_IP(rr, pRrd, 2*2);
                        rr= BBE_REPNX16C(rr,0);
                        //B-= rr*xx
                        BBE_MULSNX16C(B0,xx0,rr);
                        BBE_MULSNX16C(B1,xx1,rr);
                    }
                    // pRrd++
                    pRrd=(const xb_vecNx16*)XT_ADDX2(N*(N-1)-2*(N-k-1),(uintptr_t)pRrd);
                    // pXrd++
                    pXrd=(const xb_vecNx16*)XT_ADDX2(SX - (N-k-1)*2*P,(uintptr_t)pXrd);
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

                    // result rounding
                    B0 = BBE_RNDADJNX40(B0,q);
                    B1 = BBE_RNDADJNX40(B1,q);
                    // store res
                    b0 = BBE_PACKVNX40(B0,q);
                    b1 = BBE_PACKVNX40(B1,q);
                    pXwr_tmp= pXwr;
                    BBE_SAVNX16_XP(b0,Z_align,pXwr_tmp,4*p);
                    BBE_SAVNX16_XP(b1,Z_align,pXwr_tmp,4*p-2*BBE_SIMD_WIDTH);
                    BBE_SAPOS_FP(Z_align,pXwr_tmp);
                    pXwr= (xb_vecNx16*)XT_ADDX2(SX, (uintptr_t)pXwr);
                }
            }
        }
    }
}
#endif

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Matrix sizes SR,SD,SY,SX are selected as usual for complex block ordered 
matrix sequencies, i.e. total size is rounded up to the closest bigger 
multiple of BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the 
closest bigger multiple of degree of 2. 
SR=size(((N+1)*N)/2)
SD=size(N)
SY=size(N*P)
SX=size(N*P)
Scratch size in bytes is defined by cholbkwxxx_getScratchSize()

Input:
 N            Matrix dimension (number of columns and rows in 
              matrices R)
 P            Number of columns in right-side matrices B
 L            Number of matrices
 R[L][SR]     Sequence of L upper triangular complex matrices R
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format
 y[L][SY]     Sequence of intermediate decision matrices y
 qA,qX,qY     Fixed point representation of matrices A(or R which 
              is the same), x and y
Output:         
x[L][SX]      sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Number of matrices L must be positive
3. Matrix sizes M,N,P must be positive
4. M and N must be multiples of 4
5. qX+qA-qY must be <=16 
---------------------------------------------------------------------------*/

void cholbkwnxpn ( void    * pScr,
              complex_fract16* restrict _x, 
        const complex_fract16* restrict _R,
        const complex_fract16* restrict _D,
        const complex_fract16* restrict _y, 
                   int qA, int qY, int qX,
                   int N,int P, int L)
{
          int16_t * restrict x=(      int16_t *)_x;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict y=(const int16_t *)_y;
    int16_t* Rt=(int16_t*)pScr;
    int qXYA=qX-qY+qA;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    cholnTransformR(Rt,R,N,L);
    if (P!=1)
    {
        cholnBkwnxp(x,Rt,D,y,qXYA,N,P,L);
    }
    else
    {
        cholnBkwnx1(x,Rt,D,y,qXYA,N,L);
    }
} /* cholbkwnxpn() */

size_t cholbkwnxpn_getScratchSize (int N,int P,int L)
{
    NASSERT(L>0);
    return ((L*N*(N - 1) + BBE_SIMD_WIDTH)*sizeof(int16_t));
} /* cholbkwnxpn_getScratchSize() */

#else
DISCARD_FUN(void,cholbkwnxpn,(
            void *pScr,
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int N,int P, int L))
size_t cholbkwnxpn_getScratchSize(int N,int P,int L) { (void)N,(void)P,(void)L; return 0; }
#endif

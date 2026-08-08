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
    BBE32 code for Cholesky backward recursion, block format
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

#if 0
//#include "baseop.h"
#define MIN(x,y) ((x)<(y)?(x):(y))

static int32_t mul32x16su (int32_t x, uint16_t y)
{
  return (int32_t)((((int64_t)x)*y)>>16);
}

static  int16_t packr1x40(int64_t x,int rsh)
{
    int64_t rnd;
    rnd= (rsh<0) ? 0: (1<<rsh)>>1;
    if (rsh<0) x=(x<<-rsh);
    else       x=(x+rnd)>>rsh;
    if(x>MAX_INT16) x=MAX_INT16;
    if(x<MIN_INT16) x=MIN_INT16;
    return (int16_t)x;
}
#endif

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

/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
/*
    backward recursion: P==1
*/
void cholnBkwnx1( int16_t* restrict x, 
            const int16_t* restrict Rt,
            const int16_t* restrict D,
            const int16_t* restrict y, 
                  int qXYA,
                  int N,int L)
{
    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;
    valign R_align;
    valign X_align;
    xb_vecNx16 coef;
    xb_vecNx16 xx0,xx1,xx2,xx3;
    xb_vecNx16 rr0,rr1,rr2,rr3;
    const vsaN sh16= BBE_MOVVSA32(16);

    const int16_t   * restrict pYrd;
    const xb_vecNx16* restrict pXrd;
    const xb_vecNx16* restrict pXrd0;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;
    int m,k,M;
    int l;
    int SX=2*getSpace(N);
    int SD=2*getSpace(N);
#if 0
    int16_t r_re,r_im;
    int16_t x_re,x_im;
    int32_t A_re,A_im;
    int64_t B_re,B_im;
#endif
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    coef = BBE_SLSNX16(1,qXYA);

    k=N-1;
    M=0;    /* number of iterations in the innermost loop */
    for (; k>=XT_MAX(N-8-1,0); k--,M++)
    {
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd0= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            pXrd=pXrd0;
            BBE_LPNX16_XP(b,pXrd0,2*SX);
            pRrd = (const xb_vecNx16*)(Rt+l*N*(N-1)+(N-k-1)*(N-k-2));
            // calculate y(m,:)-R(m,:)*X, 1xP
            BBE_LPNX16_XP(b,pYrd,SX*2);// load 16 bit complex value (16 bit re and 16 bit im)
            B=BBE_MULNX16(b,coef);
            // load xx (load only M elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx0,X_align,pXrd,M*2*2);
            // load rr (load only M elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr0,R_align,pRrd,M*2*2);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx0,rr0);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // 32x16 multiply
            tl = BBE_PACKLNX40(B);
            th = BBE_PACKVNX40(B,sh16);
            BBE_LPNX16_XP(dd,pDrd,2*SD);// load D
            d0 = BBE_REPNX16(dd,0);// replicate dd[0]
            B = BBE_MULUUNX16(d0,tl);
            B = BBE_SRAINX40(B,16);
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
    for (; k>=XT_MAX(N-16-1,0); k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const int16_t* pRt=Rt+l*N*(N-1)+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); B_re = qXYA>0 ? B_re<<qXYA:B_re>>-qXYA;
            B_im=(y[l*SX+2*k+1]); B_im = qXYA>0 ? B_im<<qXYA:B_im>>-qXYA;
            NASSERT(M>4 && M<=8);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=L_mul_ss(x_re, r_re)-L_mul_ss(x_im, r_im);    // ->qX+qR
                B_im-=L_mul_ss(x_re, r_im)+L_mul_ss(x_im, r_re); 
            }
            A_re=mul32x16su((int32_t)B_re,D[l*SD+2*k+0]); 
            A_im=mul32x16su((int32_t)B_im,D[l*SD+2*k+0]); 
            x[l*SX+2*k+0]=packr1x40(A_re,D[l*SD+2*k+1]);
            x[l*SX+2*k+1]=packr1x40(A_im,D[l*SD+2*k+1]);
        }
#else
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd0= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            pXrd=pXrd0;
            BBE_LPNX16_XP(b,pXrd0,2*SX);
            pRrd = (const xb_vecNx16*)(Rt+l*N*(N-1)+(N-k-1)*(N-k-2));
            // calculate y(m,:)-R(m,:)*X, 1xP
            BBE_LPNX16_XP(b,pYrd,SX*2);// load 16 bit complex value (16 bit re and 16 bit im)
            B=BBE_MULNX16(b,coef);
            // load xx (load only M elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LANX16_IP (xx0,X_align,pXrd);
            BBE_LAVNX16_XP(xx1,X_align,pXrd,M*2*2-2*BBE_SIMD_WIDTH);
            // load rr (load only M elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LANX16_IP (rr0,R_align,pRrd);
            BBE_LAVNX16_XP(rr1,R_align,pRrd,M*2*2-2*BBE_SIMD_WIDTH);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx0,rr0);
            BBE_MULSNX16C(B,xx1,rr1);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // 32x16 multiply
            tl = BBE_PACKLNX40(B);
            th = BBE_PACKVNX40(B,sh16);
            BBE_LPNX16_XP(dd,pDrd,2*SD);// load D
            d0 = BBE_REPNX16(dd,0);// replicate dd[0]
            B = BBE_MULUUNX16(d0,tl);
            B = BBE_SRAINX40(B,16);
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
#endif
    }
    for (; k>=XT_MAX(N-32-1,0); k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const int16_t* pRt=Rt+l*N*(N-1)+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); B_re = qXYA>0 ? B_re<<qXYA:B_re>>-qXYA;
            B_im=(y[l*SX+2*k+1]); B_im = qXYA>0 ? B_im<<qXYA:B_im>>-qXYA;
            NASSERT(M>8 && M<=16);
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=L_mul_ss(x_re, r_re)-L_mul_ss(x_im, r_im);    // ->qX+qR
                B_im-=L_mul_ss(x_re, r_im)+L_mul_ss(x_im, r_re); 
            }
            A_re=mul32x16su((int32_t)B_re,D[l*SD+2*k+0]); 
            A_im=mul32x16su((int32_t)B_im,D[l*SD+2*k+0]); 
            x[l*SX+2*k+0]=packr1x40(A_re,D[l*SD+2*k+1]);
            x[l*SX+2*k+1]=packr1x40(A_im,D[l*SD+2*k+1]);
        }
#else
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd0=(const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            pXrd=pXrd0;
            BBE_LPNX16_XP(b,pXrd0,2*SX);
            pRrd = (const xb_vecNx16*)(Rt+l*N*(N-1)+(N-k-1)*(N-k-2));
            // calculate y(m,:)-R(m,:)*X, 1xP
            BBE_LPNX16_XP(b,pYrd,SX*2);// load 16 bit complex value (16 bit re and 16 bit im)
            B=BBE_MULNX16(b,coef);
            // load xx (load only M elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LANX16_IP (xx0,X_align,pXrd);
            BBE_LANX16_IP (xx1,X_align,pXrd);
            BBE_LAVNX16_XP(xx2,X_align,pXrd,M*2*2-4*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(xx3,X_align,pXrd,M*2*2-6*BBE_SIMD_WIDTH);
            // load rr (load only M elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LANX16_IP (rr0,R_align,pRrd);
            BBE_LANX16_IP (rr1,R_align,pRrd);
            BBE_LAVNX16_XP(rr2,R_align,pRrd,M*2*2-4*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(rr3,R_align,pRrd,M*2*2-6*BBE_SIMD_WIDTH);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx0,rr0);
            BBE_MULSNX16C(B,xx1,rr1);
            BBE_MULSNX16C(B,xx2,rr2);
            BBE_MULSNX16C(B,xx3,rr3);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // 32x16 multiply
            tl = BBE_PACKLNX40(B);
            th = BBE_PACKVNX40(B,sh16);
            BBE_LPNX16_XP(dd,pDrd,2*SD);// load D
            d0 = BBE_REPNX16(dd,0);// replicate dd[0]
            B = BBE_MULUUNX16(d0,tl);
            B = BBE_SRAINX40(B,16);
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
#endif
    }
    for (; k>=0; k--,M++)
    {
#if 0
        for(l=0; l<L; l++)
        {
            const int16_t* pRt=Rt+l*N*(N-1)+(N-k-1)*(N-k-2);
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re=(y[l*SX+2*k+0]); B_re = qXYA>0 ? B_re<<qXYA:B_re>>-qXYA;
            B_im=(y[l*SX+2*k+1]); B_im = qXYA>0 ? B_im<<qXYA:B_im>>-qXYA;
            for (m=0; m<M; m++)
            {
                x_re=x[l*SX+2*(k+1+m)+0];
                x_im=x[l*SX+2*(k+1+m)+1];
                r_re=pRt[2*m+0];
                r_im=pRt[2*m+1];
                B_re-=L_mul_ss(x_re, r_re)-L_mul_ss(x_im, r_im);    // ->qX+qR
                B_im-=L_mul_ss(x_re, r_im)+L_mul_ss(x_im, r_re); 
            }
            A_re=mul32x16su((int32_t)B_re,D[l*SD+2*k+0]); 
            A_im=mul32x16su((int32_t)B_im,D[l*SD+2*k+0]); 
            x[l*SX+2*k+0]=packr1x40(A_re,D[l*SD+2*k+1]);
            x[l*SX+2*k+1]=packr1x40(A_im,D[l*SD+2*k+1]);
        }
#else
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd0=(const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            int nbytes;
            pXrd=pXrd0;
            BBE_LPNX16_XP(b,pXrd0,2*SX);
            pRrd = (const xb_vecNx16*)(Rt+l*N*(N-1)+(N-k-1)*(N-k-2));
            // calculate y(m,:)-R(m,:)*X, 1xP
            BBE_LPNX16_XP(b,pYrd,SX*2);// load 16 bit complex value (16 bit re and 16 bit im)
            B=BBE_MULNX16(b,coef);
            // load xx (load only M elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LANX16_IP (xx0,X_align,pXrd);
            BBE_LANX16_IP (xx1,X_align,pXrd);
            BBE_LANX16_IP (xx2,X_align,pXrd);
            BBE_LANX16_IP (xx3,X_align,pXrd);
            // load rr (load only M elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LANX16_IP (rr0,R_align,pRrd);
            BBE_LANX16_IP (rr1,R_align,pRrd);
            BBE_LANX16_IP (rr2,R_align,pRrd);
            BBE_LANX16_IP (rr3,R_align,pRrd);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx0,rr0);
            BBE_MULSNX16C(B,xx1,rr1);
            BBE_MULSNX16C(B,xx2,rr2);
            BBE_MULSNX16C(B,xx3,rr3);
            nbytes=M*4-8*BBE_SIMD_WIDTH;
            for (m=0; m<M-32; m++)
            {
                BBE_LAVNX16_XP (xx0,X_align,pXrd,nbytes);
                BBE_LAVNX16_XP (rr0,R_align,pRrd,nbytes);
                BBE_MULSNX16C(B,xx0,rr0);
                nbytes-=2*BBE_SIMD_WIDTH;
            }
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // 32x16 multiply
            tl = BBE_PACKLNX40(B);
            th = BBE_PACKVNX40(B,sh16);
            BBE_LPNX16_XP(dd,pDrd,2*SD);// load D
            d0 = BBE_REPNX16(dd,0);// replicate dd[0]
            B = BBE_MULUUNX16(d0,tl);
            B = BBE_SRAINX40(B,16);
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
#endif
    }
}
#endif

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
    BBE32 code for Cholesky forward recursion, block format
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
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

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
/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
// N>24 && N<=32
void cholnFwdrec32(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
#if 1
{
    int l;
    int n;
    const int32_t* restrict pZrd;
    const xb_vecNx16* restrict pYrd;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;
    int16_t* restrict pYwr;
    const xb_vecNx16 *pYrd0;
    xb_vecNx16 dd,d0,rr0,rr1,rr2,rr3,yy0,yy1,yy2,yy3,tl,th,b_res,b;
    xb_vecNx40 B;
    xb_c40 r_summ;
    valign align,ay;
    vsaN q;
    const vsaN sh16= BBE_MOVVSA32(16);
    NASSERT(N>24 && N<=32);
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));

    cholnFwdrec24(y,R,D,Z,24,L,SR,SD,SY,SZ);
    for(n=24; n<N; n++)
    {
        pYrd = (const xb_vecNx16*)y;
        //pYwr = &y[n+n];
        pYwr= (int16_t*)XT_ADDX4(n,(uintptr_t)y);
        //pZrd = &Z[n+n];
        pZrd= (const int32_t*)XT_ADDX8(n,(uintptr_t)Z);
        pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
        //pDrd = (const xb_vecNx16*)&D[n+n];
        pDrd= (const xb_vecNx16*)XT_ADDX4(n,(uintptr_t)D);
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            // calculate A(:,n)'*B-Rn'*Y, 1xP
            // load 32 bit complex value (32 bit re and 32 bit im)
            // load B
            b = BBE_LV4X16_I(pZrd,0);
            B=BBE_MOVSWVL(b);
            // pZrd+=N*2
            pZrd=(const int32_t*)XT_ADDX4(SZ,(uintptr_t)pZrd);
            // load rr (load only n elements. others = 0)
            align = BBE_LA_PP(pRrd);
            BBE_LANX16_IP (rr0,align,pRrd);
            BBE_LANX16_IP (rr1,align,pRrd);
            BBE_LANX16_IP (rr2,align,pRrd);
            BBE_LAVNX16_XP(rr3,align,pRrd,n*4-6*BBE_SIMD_WIDTH);
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
            pYrd0=pYrd;
            ay=BBE_LA_PP(pYrd);
            BBE_LANX16_IP(yy0,ay,pYrd);
            BBE_LANX16_IP(yy1,ay,pYrd);
            BBE_LANX16_IP(yy2,ay,pYrd);
            BBE_LANX16_IP(yy3,ay,pYrd);
            pYrd=(const xb_vecNx16 *)(((uintptr_t)pYrd0)+SY*2);
            //B-= rr*yy
            BBE_MULSNX16J(B,yy0,rr0);
            BBE_MULSNX16J(B,yy1,rr1);
            BBE_MULSNX16J(B,yy2,rr2);
            BBE_MULSNX16J(B,yy3,rr3);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // get low 16 bits of B
            tl = BBE_PACKLNX40(B);
            // get high 16 bits of B
            th = BBE_PACKVNX40(B,sh16);
            // load D, pDrd+= SD*2
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
            BBE_SPNX16_XP(b_res,pYwr,2*SY);
        }
        //__Pragma("no_reorder")
    }
}
#else
{
    int l;
    int n,m;
    const int16_t *pR;
    //NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Z,2*BBE_SIMD_WIDTH);

    NASSERT(N<=32);

    for(n=0; n<N; n++)
    {
        for(l=0; l<L; l++)
        {
            int64_t B_re,B_im;
            pR=R+(n*(n+1))+l*SR;
            // calculate A(:,n)'*B-Rn'*Y, 1xP
            B_re=Z[SZ*l+2*n+0];
            B_im=Z[SZ*l+2*n+1];
            for (m=0; m<n; m++)   
            {
                int16_t r_re,r_im;
                int16_t y_re,y_im;
                r_re=pR[m*2+0];
                r_im=pR[m*2+1];
                y_re=y[l*SY+m*2+0];
                y_im=y[l*SY+m*2+1];
                B_re-=L_mul_ss(y_re,r_re)+L_mul_ss(y_im,r_im); // representation qA+qY
                B_im-=L_mul_ss(y_im,r_re)-L_mul_ss(y_re,r_im);
            }
            B_re=mul32x16su((int32_t)B_re,D[l*SD+2*n+0]); // representation (qA+qY),qD -> qA+qY+qD-31
            B_im=mul32x16su((int32_t)B_im,D[l*SD+2*n+0]); // 
            y[l*SY+m*2+0]=packr1x40(B_re,D[l*SD+2*n+1]);
            y[l*SY+m*2+1]=packr1x40(B_im,D[l*SD+2*n+1]);
        }
    }
}
#endif

#endif

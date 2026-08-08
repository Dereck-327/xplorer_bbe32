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
  NatureDSP_Baseband library. Cholesky backward recursion for block ordered matrices:
    These functions make backward recursion stage of pseudo-inversion. They use
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
#define GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(size) ((size)+(2*BBE_SIMD_WIDTH)-1)&(~(2*BBE_SIMD_WIDTH-1))

/* reversing R matrices for easier readings by rows (diagonal elements are omitted):
   original R    transformed R
   0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
     2 4 7 b
       5 8 c
         9 d
           e

   Input:
   R[L][SR]        L input matrices
   Rt[L*N*(N-1)]   stream of L trasposed matrices
*/
static void transformR8x8(int16_t* restrict  Rt,const int16_t  * restrict R,int L)
#if 0
{
    int k,l;
    /* 
        matlab code for generation of permutation indexes:
        N=8;
        a=zeros(N,N);
        m=1;
        for n=1:N
            a(1:n,n)=(m:m+n-1);
            m=m+n;
        end
        % read by rows from the end except for the diagonal elements
        b=[];
        k=1;m=1;
        for n=N-1:-1:1
            b(m:m+k-1)=a(n,n+1:end);
            m=m+k;
            k=k+1;
        end
        fprintf(1,'%d ',b-1);
    */
    static const int perm[28]={34,26,33,19,25,32,13,18,24,31,8,12,17,23,30,4,7,11,16,22,29,1,3,6,10,15,21,28};
    for(l=0; l<L; l++)
    {
        for (k=0; k<28; k++)
        {
            int ix;
            ix=perm[k];
            Rt[2*k+0]=R[2*ix+0];
            Rt[2*k+1]=R[2*ix+1];
        }
        R+=80;
        Rt+=56;
    }
}
#else
{
    /*
    	packed (2x8 bit --> 1x16 bit) selection indexes of R and intermediate selection
    */
    static
    const int16_t ALIGN(32) mask2[] = { 4,      261,    20,     277,    2,      259,    1536,   1793,   18,     275,    0,      257,    6656,   6913,   1024,   1281,
                                        0,      257,    3586,   3843,   4,      261,    22,     279,    4616,   4873,   7690,   7947,   3100,   3357,   30,     287,
                                        0,      257,    512,    769,    5120,   5377,   5640,   5897,   2048,   2305,   2560,   2817,   3072,   3329,   7704,   7961,
                                        14,     271,    22,     279,    0,      257,    3072,   3329,   6656,   6913,   2,      259,    6,      263,    12,     269,
                                        1024,   1281,   3586,   3843,   6676,   6933,   22,     279,    24,     281,    10,     267,    12,     269,    14,     271,
                                        0,      1,      2,      3,      4,      5,      24,     25,     0,      1,      0,      1,      0,      1,      0,      1};

    xb_vecNx16 RR0, RR1, RR2, RR3, RR4;
    xb_vecNx16 SEL0_1, SEL2_3, SEL4_5, SEL6_7, SEL8_9, SEL10_11;
    xb_vecNx16 RRT;
    vboolN _true=BBE_LTRNI(0);

    xb_vecNx16 RTMP0, RTMP1;
    valign align;

    vselN sel0, sel1, sel2, sel3, sel4, sel5, sel6, sel7, sel8, sel9, sel10;
    int l;
    const xb_vecNx16 *restrict ppR = (const xb_vecNx16*)R;
    xb_vecNx16 *restrict pRt = (xb_vecNx16*)Rt;
    SEL0_1=     BBE_LVNX16_I((const xb_vecNx16*)mask2, 0*2*BBE_SIMD_WIDTH);
    SEL2_3=     BBE_LVNX16_I((const xb_vecNx16*)mask2, 1*2*BBE_SIMD_WIDTH);
    SEL4_5=     BBE_LVNX16_I((const xb_vecNx16*)mask2, 2*2*BBE_SIMD_WIDTH);
    SEL6_7=     BBE_LVNX16_I((const xb_vecNx16*)mask2, 3*2*BBE_SIMD_WIDTH);
    SEL8_9=     BBE_LVNX16_I((const xb_vecNx16*)mask2, 4*2*BBE_SIMD_WIDTH);
    SEL10_11=   BBE_LVNX16_I((const xb_vecNx16*)mask2, 5*2*BBE_SIMD_WIDTH);

    sel0 = BBE_MOVVSELNX16(SEL0_1,0);
    sel1 = BBE_MOVVSELNX16(SEL0_1,8);

    align = BBE_ZALIGN();
    __Pragma("loop_count min=1");
    for(l=0; l<L; l++)
    {
        // load entire R
        BBE_LVNX16_IP(RR0, ppR,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(RR1, ppR,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(RR2, ppR,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(RR3, ppR,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(RR4, ppR,2*BBE_SIMD_WIDTH);
        // permutation
        // get RRT0
        //RR4 - 0..15, RR3 - 16..31
        RTMP0 = BBE_SELNX16(RR3, RR4, sel0);
        RTMP1 = BBE_SELNX16(RR1, RR2, sel1);
        sel2 = BBE_MOVVSELNX16(SEL2_3,0);
        RRT = BBE_SELNX16(RTMP1, RTMP0, sel2);
        BBE_SAVNX16_XP(RRT,align,pRt,2*BBE_SIMD_WIDTH);

        // get RRT1
        sel3 = BBE_MOVVSELNX16(SEL2_3,8);
        RTMP0 = BBE_SELNX16(RR2, RR3, sel3);
        sel4 = BBE_MOVVSELNX16(SEL4_5,0);
        RTMP1 = BBE_SELNX16(RR0, RR1, sel4);
        sel5=BBE_MOVVSELNX16(SEL4_5,8);
        RRT = BBE_SELNX16(RTMP1, RTMP0, sel5);
        BBE_SAVNX16_XP(RRT,align,pRt,2*BBE_SIMD_WIDTH);

        SEL2_3=BBE_MOVNX16T(SEL2_3,SEL2_3,_true);
        SEL4_5=BBE_MOVNX16T(SEL4_5,SEL4_5,_true);

        // get RRT2
        sel6=BBE_MOVVSELNX16(SEL6_7,0);
        RTMP0 = BBE_SELNX16(RR1, RR0, sel6);
        sel7=BBE_MOVVSELNX16(SEL6_7,8);
        RTMP1 = BBE_SELNX16(RR3, RR2, sel7);
        sel8=BBE_MOVVSELNX16(SEL8_9,0);
        RRT = BBE_SELNX16(RTMP1, RTMP0, sel8);
        BBE_SAVNX16_XP(RRT,align,pRt,2*BBE_SIMD_WIDTH);

        SEL6_7=BBE_MOVNX16T(SEL6_7,SEL6_7,_true);
        SEL8_9=BBE_MOVNX16T(SEL8_9,SEL8_9,_true);

        // get RRT3
        sel9=BBE_MOVVSELNX16(SEL8_9,8);
        RTMP0 = BBE_SELNX16(RR2, RR1, sel9);
        sel10=BBE_MOVVSELNX16(SEL10_11,0);
        RRT = BBE_SELNX16(RR3, RTMP0, sel10);
        BBE_SAVNX16_XP(RRT,align,pRt,2*(BBE_SIMD_WIDTH/2));

        SEL8_9=BBE_MOVNX16T(SEL8_9,SEL8_9,_true);
        SEL10_11=BBE_MOVNX16T(SEL10_11,SEL10_11,_true);

    }
    BBE_SAPOS_FP(align,pRt);
}
#endif

/*
    Convert Diagonal D from block to streaming format
    
    Input:
    D_in            D[L][8][2] diagonal elements in block format
    D_out           D[8][L][2] diagonal element in streaming format
    L               number of matrices
    restrictions:
    L is multiple of 8
*/
static void transformD8( int16_t* restrict D_out,
                   const int16_t* restrict D_in,
                   const int L)
{
    #define N 8
    #define SD 16

    NASSERT((L%8)==0);

    int l;
    const int stride= (L >> 3)*2*BBE_SIMD_WIDTH;
    const xb_vecNx16* restrict pD_in_rd= (const xb_vecNx16*)D_in;
    xb_vecNx16* restrict pD_out_wr= (xb_vecNx16*)D_out;
    xb_vecNx16 vX0, vX1, vX2, vX3, vX4, vX5, vX6, vX7;

    __Pragma("loop_count min=1")
    for (l=0; l<L; l+=8)
    {
        // read 8 D*8*2 elements
        BBE_LVNX16_IP(vX0, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX1, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX2, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX3, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX4, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX5, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX6, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX7, pD_in_rd, 2*BBE_SIMD_WIDTH);

        // complex transpose
        BBE_DSELNX16I(vX1,vX0,vX1,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX3,vX2,vX3,vX2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX5,vX4,vX5,vX4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX6,vX7,vX6,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vX2,vX0,vX2,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX3,vX1,vX3,vX1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX6,vX4,vX6,vX4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX5,vX7,vX5,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vX4,vX0,vX4,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX5,vX1,vX5,vX1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX6,vX2,vX6,vX2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX3,vX7,vX3,BBE_DSELI_DEINTERLEAVE_2);
        // now transposed matrix in vX0, vX1, ..., vX7
        
        // store res
        BBE_SVNX16_XP(vX0,pD_out_wr,stride);
        BBE_SVNX16_XP(vX1,pD_out_wr,stride);
        BBE_SVNX16_XP(vX2,pD_out_wr,stride);
        BBE_SVNX16_XP(vX3,pD_out_wr,stride);
        BBE_SVNX16_XP(vX4,pD_out_wr,stride);
        BBE_SVNX16_XP(vX5,pD_out_wr,stride);
        BBE_SVNX16_XP(vX6,pD_out_wr,stride);
        BBE_SVNX16_XP(vX7,pD_out_wr,-7*stride+2*BBE_SIMD_WIDTH);
    }
    #undef N
    #undef SD
}

/*
    backward recursion: P==1
    Input:
    Rt[L][8*8][2]       L transformed R matrices
    D[L][SD][2]         reciprocal of main diagonal (mantissa, exponent) 
                        in the special format
    y[L][SY][2]         sequence of intermediate decision matrices y
    qA,qX,qY            fixed point representation of matrices A(or R which 
                        is the same), x and y
    L                   number of matrices
    Output:
    x[L][SX][2]         sequence of decision matrix x

    Restrictions: L= 1
*/
static void bkw8x1(int16_t* restrict x, 
                    const int16_t* restrict Rt,
                    const int16_t* restrict D,
                    const int16_t* restrict y, 
                    int qXYA)
{
    #define N 8
    #define SX 16
    #define SD 16

    int k;
    
    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;
    valign R_align;
    valign X_align;

    const int16_t* restrict pYrd;
    const xb_vecNx16* restrict pXrd;
    xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    xb_vecNx16 coef;
    xb_vecNx16 xx;
    xb_vecNx16 rr;

    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));

    // load shift amount
    coef = 1;
    coef = BBE_SLSNX16(coef,qXYA);

    const vsaN sh16= BBE_MOVVSA32(16);
    for (k=N-1; k>=0; k--)
    {
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        pRrd = (const xb_vecNx16*)(Rt+(N-k-1)*(N-k-2));
        // calculate y(m,:)-R(m,:)*X, 1xP

        // load 16 bit complex value (16 bit re and 16 bit im)
        // load B and shl with saturation
        BBE_LPNX16_IP(b,pYrd,SX*2);
        B=BBE_MULNX16(b,coef);
        // load xx (load only (N-k-1) elements. others = 0)
        X_align = BBE_LA_PP(pXrd);
        BBE_LAVNX16_XP(xx,X_align,pXrd,(N-k-1)*2*2);
        
        // load rr (load only (N-k-1) elements. others = 0)
        R_align = BBE_LA_PP(pRrd);
        BBE_LAVNX16_XP(rr,R_align,pRrd,(N-k-1)*2*2);

        //B-= rr*xx
        BBE_MULSNX16C(B,xx,rr);
        // reduced add
        r_summ = BBE_RADDNX40C(B);
        // type conversion
        B=BBE_MOVNX40_FROMC40(r_summ);
        // get low 16 bits of B
        tl = BBE_PACKLNX40(B);
        // get high 16 bits of B
        th = BBE_PACKVNX40(B,sh16);
        // load D
        BBE_LPNX16_IP(dd,pDrd,2*SD);
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
        BBE_SPNX16_IP(b_res,pXwr,2*SX);
    }
    #undef SX
    #undef N
    #undef SD
}

/*
    backward recursion: P==1
    Input:
    Rt[L][8*8][2]       L transformed R matrices
    D[L][SD][2]         reciprocal of main diagonal (mantissa, exponent) 
                        in the special format
    y[L][SY][2]         sequence of intermediate decision matrices y
    qA,qX,qY            fixed point representation of matrices A(or R which 
                        is the same), x and y
    L                   number of matrices
    Output:
    x[L][SX][2]         sequence of decision matrix x
    Restrictions:
    L is multiple of 2
*/
static void bkw8x1_2(int16_t* restrict x, 
                        const int16_t* restrict Rt,
                        const int16_t* restrict D,
                        const int16_t* restrict y, 
                        int qXYA,
                        int L)
{
    #define N 8
    #define SX 16
    #define SD 16

    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT((L%2)==0);

    int k;
    int l;

    static
    const int16_t ALIGN(32) Sel[16]= {((1 << 8) + 0),((1 << 8) + 0),0,0,	0,0,0,0,	((17 << 8) + 16),((17 << 8) + 16),0,0,	0,0,0,0};
    xb_vecNx16 tl,th,dd0, dd1,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;
    valign X_align;
    valign R_align;

    const int16_t* restrict pYrd;
    const xb_vecNx16* restrict pXrd;
    xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    xb_vecNx16 coef;
    xb_vecNx16 xx;
    xb_vecNx16 rr;
    xb_vecNx16 vTmp0, vTmp1;

    vTmp0= BBE_LVNX16_I((const xb_vecNx16*)Sel, 0);
    const vselN sel_d0= BBE_MOVVSELNX16(vTmp0,0);
    const vselN sel_q= BBE_MOVVSELNX16(vTmp0,8);

    // load shift amount
    coef = 1;
    coef = BBE_SLSNX16(coef,qXYA);
    const vsaN sh16= BBE_MOVVSA32(16);

    for (k=N-1; k>=0; k--)
    {
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        pRrd = (const xb_vecNx16*)(Rt+(N-k-1)*(N-k-2));

        __Pragma("loop_count min=1");
        for(l=0; l<L; l+=2)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP

            // load 16 bit complex value (16 bit re and 16 bit im)
            // load B and shl with saturation
            BBE_LPNX16_IP(b,pYrd,SX*2);
            B=BBE_MULNX16(b,coef);
            // load xx
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx,X_align,pXrd,2*SX);
            // load rr (load only N-k-1 elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-k-1)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2((N*(N-1))-((N-k-1)*2),(uintptr_t)pRrd);
            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp0= BBE_MOVVWL(B);
            //--------------------------------------------------------------------------
            BBE_LPNX16_IP(b,pYrd,SX*2);
            B=BBE_MULNX16(b,coef);
            // load xx (load only N-k-1 elements. others = 0)
            BBE_LAVNX16_XP(xx,X_align,pXrd,(N-k-1)*2*2);
            pXrd=(const xb_vecNx16*)XT_ADDX2(SX-((N-k-1)*2),(uintptr_t)pXrd);
            // load rr (load only N-k-1 elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-k-1)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2((N*(N-1))-((N-k-1)*2),(uintptr_t)pRrd);
            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp1= BBE_MOVVWL(B);
            //--------------------------------------------------------
            B= BBE_MOVSWV(vTmp1,vTmp0);
            // get low 16 bits of B
            tl = BBE_PACKLNX40(B);
            // get high 16 bits of B
            th = BBE_PACKVNX40(B,sh16);
            // load D, pDrd+= SD*2
            BBE_LPNX16_IP(dd0,pDrd,2*SD);
            BBE_LPNX16_IP(dd1,pDrd,2*SD);
            // replicate dd[0]
            d0= BBE_SELNX16(dd1, dd0, sel_d0);
            // make q
            vTmp0= BBE_SELNX16(dd1, dd0, sel_q);
            q = BBE_MOVVSV(vTmp0,0);
            // mul low (unsigned*unsigned)
            B = BBE_MULUUNX16(d0,tl);
            // shift low left
            B = BBE_SRAINX40(B,16);
            // Bres+= mul high (unsigned*signed)
            BBE_MULUSANX16(B,d0,th);

            // result rounding
            B = BBE_RNDADJNX40(B,q);
            // store res
            b_res = BBE_PACKVNX40(B,q);
            BBE_SPNX16_IP(b_res,pXwr,2*SX);
            b_res= BBE_REPNX16C(b_res,4);
            BBE_SPNX16_IP(b_res,pXwr,2*SX);
        }
    }
    #undef SD
    #undef SX
    #undef N
}

/*
    backward recursion: P==1
    Input:
    Rt[L][8*8][2]       L transformed R matrices
    D[L][SD][2]         reciprocal of main diagonal (mantissa, exponent) 
                        in the special format
    y[L][SY][2]         sequence of intermediate decision matrices y
    qA,qX,qY            fixed point representation of matrices A(or R which 
                        is the same), x and y
    L                   number of matrices
    Output:
    x[L][SX][2]         sequence of decision matrix x
    Restrictions:
    L is multiple of 8
*/
static void bkw8x1_8(int16_t* restrict x, 
                        const int16_t* restrict Rt,
                        const int16_t* restrict Dtr,
                        const int16_t* restrict y, 
                        int qXYA,
                        int L,
                        void* pScr)
{
    #define N 8
    #define SX 16
    #define SD 16

    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Dtr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT((L%8)==0);

    int l;
    int k;
    const int16_t* restrict pYrd;
    const xb_vecNx16* restrict pXrd;
    xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd= (const xb_vecNx16*)&Dtr[7*L*2];
    xb_vecNx16* restrict pYrdwr= (xb_vecNx16*)pScr;

    xb_vecNx16 d0,rr,xx,tl,th,b,vTmp0, vTmp1;
    xb_vecNx16 coef;

    xb_vecNx40 B;
    xb_c40 r_summ;
    valign R_align;
    valign X_align;
    valign Yrdwr_align;
    vsaN q;

    // load shift amount
    coef = 1;
    coef = BBE_SLSNX16(coef,qXYA);
    const vsaN sh16= BBE_MOVVSA32(16);

    Yrdwr_align= BBE_ZALIGN();

    for (k=N-1; k>=0; k--)
    {
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        pXrd= (const xb_vecNx16*)XT_ADDI((uintptr_t)pXwr,2*2);
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        pRrd = (const xb_vecNx16*)(Rt+(N-k-1)*(N-k-2));
        X_align = BBE_LA_PP(pXrd);
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=1 factor=2")
        for(l=0; l<L; l+=2)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP

            // load 16 bit complex value (16 bit re and 16 bit im)
            // load B and shl with saturation
            BBE_LPNX16_IP(b,pYrd,SX*2);
            B=BBE_MULNX16(b,coef);
            // load xx
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx,X_align,pXrd,2*SX);
            // load rr (load only n elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-k-1)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2((N*(N-1))-((N-k-1)*2),(uintptr_t)pRrd);
            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp0= BBE_MOVVWL(B);
            //1-------------------------------------------------------------------------
            BBE_LPNX16_IP(b,pYrd,SX*2);
            B=BBE_MULNX16(b,coef);
            // load xx (load only N-k-1 elements. others = 0)
            BBE_LAVNX16_XP(xx,X_align,pXrd,(N-k-1)*2*2);
            pXrd=(const xb_vecNx16*)XT_ADDX2(SX-((N-k-1)*2),(uintptr_t)pXrd);
            // load rr (load only n elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-k-1)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2((N*(N-1))-((N-k-1)*2),(uintptr_t)pRrd);
            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp1= BBE_MOVVWL(B);

            vTmp0= BBE_SELNX16I(vTmp1,vTmp0,BBE_SELI_PACK_4);
            BBE_SAVNX16_XP(vTmp0,Yrdwr_align,pYrdwr,2*8);
        }
        BBE_SAPOS_FP(Yrdwr_align,pYrdwr);
        pYrdwr= (xb_vecNx16*)pScr;
        __Pragma("loop_count min=1")
        for(l=0; l<L; l+= 8)
        {
            BBE_LVNX16_IP(vTmp0,pYrdwr,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(vTmp1,pYrdwr,2*BBE_SIMD_WIDTH);
            B=BBE_MOVSWV(vTmp1, vTmp0);
            // get low 16 bits of B
            tl = BBE_PACKLNX40(B);
            // get high 16 bits of B
            th = BBE_PACKVNX40(B,sh16);
            // load D
            BBE_LVNX16_IP(d0,pDrd,2*BBE_SIMD_WIDTH);
            vTmp0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_ODD);
            d0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_EVEN);
            q = BBE_MOVVSV(vTmp0,0);
            // mul low (unsigned*unsigned)
            B = BBE_MULUUNX16(d0,tl);
            // shift low left
            B = BBE_SRAINX40(B,16);
            // Bres+= mul high (unsigned*signed)
            BBE_MULUSANX16(B,d0,th);
            // result rounding
            B = BBE_RNDADJNX40(B,q);
            // store res
            vTmp0 = BBE_PACKVNX40(B,q);
            BBE_SPNX16_IP(vTmp0,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,1);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,2);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,3);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,4);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,5);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,6);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,7);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
        }
        // set pDtr to prev element
        pDrd= (const xb_vecNx16*)XT_ADDX2(-2*L*2,(uintptr_t)pDrd);
        //__Pragma("no_reorder")
    }
    #undef N
    #undef SX
    #undef SD
}

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

void cholbkw8x1n ( void    * pScr,
          complex_fract16* restrict _x, 
    const complex_fract16* restrict _R,
    const complex_fract16* restrict _D,
    const complex_fract16* restrict _y, 
                   int qA, int qY, int qX,
                   int L )
{
          int16_t * restrict x=(      int16_t *)_x;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict y=(const int16_t *)_y;
    #define N 8
    #define SX 16
    #define SD 16
    #define SRt (3*BBE_SIMD_WIDTH + BBE_SIMD_WIDTH/2)
    int L1, L2, L8;
    int16_t* Rt=(int16_t*)pScr;
    int qXYA=qX-qY+qA;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    transformR8x8(Rt,R,L);
    
    // calculate values that less or equal than L and multiple of 8 (L8), multiple of 2(L2), multiple of 1 (L1)
    L8= L - (L&7);
    L2= (L - L8) - (L&1);
    L1= L - L8 - L2;

    if (L8)
    {
        size_t Rt_size  = GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(((L)*(N)*((N)-1))*sizeof(int16_t));
        size_t D_tr_size= GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(2*L*N*sizeof(int16_t));
        int16_t* D_tr= (int16_t*)((char*)pScr + Rt_size);
        pScr= ((char*)pScr + Rt_size + D_tr_size);
        transformD8(D_tr,D,L8);
        bkw8x1_8(x,Rt,D_tr,y,qXYA,L8,pScr);
    }
    if (L2)
    {
              int16_t* x2= x + L8*SX;
        const int16_t* Rt2= Rt+ L8*SRt;
        const int16_t* D2= D + L8*SD;
        const int16_t* y2= y + L8*SX;
        bkw8x1_2(x2,Rt2,D2,y2,qXYA,L2);
    }
    if (L1)
    {
              int16_t* x1= x + (L8+L2)*SX;
        const int16_t* Rt1= Rt+ (L8+L2)*SRt;
        const int16_t* D1= D + (L8+L2)*SD;
        const int16_t* y1= y + (L8+L2)*SX;
        bkw8x1(x1,Rt1,D1,y1,qXYA);
    }
    
    #undef N
    #undef SX
    #undef SD
    #undef SRt
} /* cholbkw8x1n() */

size_t cholbkw8x1n_getScratchSize (int N,int P,int L)
{
    size_t Rt_size;
    size_t D_tr_size = 0;
    size_t Intrmd_size = 0;
    NASSERT(N==8 && P==1 && L>0);
    Rt_size = GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(((L)*(N)*((N)-1))*sizeof(int16_t));
    if (L>7)
    {
        D_tr_size= (2*L*N*sizeof(int16_t));
        Intrmd_size= GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(L*2*sizeof(int32_t));
    }
    return (Rt_size + D_tr_size + Intrmd_size);
} /* cholbkw8x1n_getScratchSize() */
#else
DISCARD_FUN(void,cholbkw8x1n,(
            void *pScr,
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L))
size_t cholbkw8x1n_getScratchSize(int N,int P,int L) { (void)N,(void)P,(void)L; return 0; }

#endif

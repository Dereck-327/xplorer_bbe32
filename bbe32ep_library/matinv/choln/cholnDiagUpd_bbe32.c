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
    BBE32 code, part of Cholesky decomposition, block format
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
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))


static int32_t mul16x16su (int16_t x, uint16_t y)
{
  return ((int32_t)x)*y;
}

static void invSqrt_bbe32ep(int16_t* D,int64_t A)
{
  int d;
  
  xb_vecNx40 x_vec;
  xb_vecNx40 a_vec;
  xb_vecNx16   b_vec, xn_vec, cn_vec;
  vsaN c_vec,_24;

  _24=BBE_MOVVSA32(24);
  // Load the input number, Qx
  x_vec = BBE_MOVWA32((int32_t)A);
  x_vec =BBE_SLSINX40(x_vec,5);
  c_vec=BBE_NSAENX40(x_vec);
  x_vec=BBE_SLLNX40(x_vec,c_vec);
  BBE_RSQRTLUNX40_0(a_vec,b_vec, cn_vec, x_vec );
  BBE_MULUUSNX16( a_vec, cn_vec,  b_vec);
  a_vec=BBE_SRAINX40(a_vec,23);
  xn_vec=BBE_PACKLNX40(a_vec);
  D[0]= BBE_MOVAV16(xn_vec);
  xn_vec=BBE_MOVVVS(c_vec);
  d = BBE_MOVAV16(xn_vec);
  D[1] = 16-(d>>1);
}
#endif
/*-------------------------------------------------
   update N-th diagonal element
   Input:
   Z[L][N+1][2]  convolutions in N-th column
   Input/output:
   y             pointer to the begining of column 
                 in matrix R[L][SR] (N+1 elements 
                 is written)
   Output:
   D[L][SD]      reciprocals of main diagonal 
                 (pointer to the N-th element
-------------------------------------------------*/
void cholnDiagUpd(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD)
#if 0
{
    int64_t B_re;
    int sh; 
    int l,m;
    Z+=2*N; // points to diagonal element
    for(l=0; l<L; l++)
    {
        B_re=Z[2*l*(N+1)+0];
        for (m=0; m<N; m++)   
        {
            int16_t r_re,r_im;
            r_re=y[2*m+0]; r_im=y[2*m+1];
            B_re-=L_mul_ss(r_re,r_re)+L_mul_ss(r_im,r_im);
        }
        invSqrt_bbe32ep(D,B_re);
        sh=D[1]+1;
        B_re= (sh>0) ? B_re>>sh : B_re<<(-sh);
        B_re=mul16x16su((int16_t)(B_re),D[0]);
        B_re=(B_re+0x4000)>>15;
        B_re=MIN(MAX(B_re,MIN_INT16),MAX_INT16);
        y[2*N+0]=(int16_t)B_re;
        y[2*N+1]=0;
        y+=SR;
        D+=SD;
    }
}
#else
{
    int l, n;

    xb_vecNx16 b, rr, d0, d1;
    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN vsaSh;
    valign align;

    const int32_t* restrict pZrd = &Z[2*N]; // points to diagonal element;
    int16_t *restrict pYwr = &y[2*N];
    const xb_vecNx16* restrict pRrd = (const xb_vecNx16*)y;
    __Pragma("loop_count min=1");
    for(l=0; l<L; l++)
    {
        // load 32 bit real value
        // load B
        BBE_LPNX16_XP(b,pZrd,2*(N+1)*4);
        B=BBE_MOVSWVL(b);
        // load rr (load only N elements. others = 0)
        align = BBE_LA_PP(pRrd);
        __Pragma("loop_count min=1");
        for (n= N*4; n>0; n-= sizeof(xb_vecNx16))
        {
            BBE_LAVNX16_XP(rr,align,pRrd,n);
            //B-= rr*rr
            BBE_MULSNX16J(B,rr,rr);
        }
        // pRrd+=SR-N*2
        pRrd=(const xb_vecNx16*)XT_ADDX2(SR-N*2,(uintptr_t)pRrd);
        // reduced add
        r_summ = BBE_RADDNX40C(B);
        // type conversion
        B=BBE_MOVNX40_FROMC40(r_summ);
        // store B
        b=BBE_MOVVWL(B);
        // invSqrt
        {
            xb_vecNx16   b_vec, cn_vec;
            B = BBE_ADDNX40(B,B);
            vsaSh=BBE_NSAENX40(B);
            B=BBE_SLLNX40(B,vsaSh);
            BBE_RSQRTLUNX40_0(B,b_vec, cn_vec, B);
            BBE_MULUUSNX16(B, cn_vec,  b_vec);
            B=BBE_SRAINX40(B,23);
            d0=BBE_PACKLNX40(B);
            vsaSh = BBE_SUBSR1SAVSN(18,vsaSh);
        }
        // store D[0]
        BBE_SSNX16_IP(d0,D,2);
        // store D[1]
        d1 = BBE_MOVVVS(vsaSh);
        BBE_SSNX16_XP(d1,D,2*SD - 2);
        // sh+= 1
        vsaSh = BBE_ADDSAVSN(1,vsaSh);
        // restore B
        B=BBE_MOVSWVL(b);
        // pack B
        b = BBE_PACKVNX40(B,vsaSh);
        //(unsigned*signed)
        B = BBE_MULUSNX16(d0,b);
        b = BBE_PACKQNX40(B);
        // store res
        BBE_SPNX16_XP(b,pYwr,SR*2);
    }
}
#endif
#endif

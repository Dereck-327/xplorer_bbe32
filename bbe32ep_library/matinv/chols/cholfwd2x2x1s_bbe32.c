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
  NatureDSP_Baseband library. Cholesky forward recursion for pseudo-inversion API (complex data)
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

#if !(HAVE_MULUUNX16PACKH && 1)
DISCARD_FUN(void, cholfwd2x2x1s,(
                  complex_fract16* restrict _y,
            const complex_fract16* restrict _R, 
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _A, 
            const complex_fract16* restrict _B, 
            int qA,int qB,int qY,
            int L))
#else

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
M               Matrix dimension (number of rows in matrices A)
N               Matrix dimension (number of columns and rows in matrices 
                R)
P               Number of columns in right-side matrices B
L               Number of matrices
R[N*N][L]       Cholesky upper triangular matrices R
A[M*N][L]       Original left-side matrices A
B[M*P][L]       Original right-side matrices B
D[L/4][N][8]    Reciprocal of main diagonal (mantissa, exponent) in the 
                special format
qA,qB,qY        Fixed point representation of matrices A (or R which is 
                the same), B and y
Output:
y[N*P][L]       Decision matrix y

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8
3. Matrix sizes (M,N) must be greater than 1, P must be >=1
4. M >= N
---------------------------------------------------------------------------*/

void cholfwd2x2x1s (
                  complex_fract16* restrict _y,
            const complex_fract16* restrict _R, 
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _A, 
            const complex_fract16* restrict _B, 
            int qA,int qB,int qY,
            int L)
{
          int16_t* restrict y=(      int16_t*)_y;
    const int16_t* restrict R=(const int16_t*)_R; 
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict A=(const int16_t*)_A; 
    const int16_t* restrict B=(const int16_t*)_B; 

/*
  Reference Matlab code
  function [Y]=cholfwd(R,A,B)
  sz=size(A); M=sz(1); N=sz(2); 
  sz=size(B); P=sz(2); 
  D=real(1./diag(R));
  AB=A'*B;
  Y=zeros(N,P);
  for n=1:N
    Rn=R(:,n); 
    Bn=AB(n,:);
    y=(Bn-Rn'*Y)*D(n);
    Y(n,:)=y;
  end
*/
  const xb_vecNx16      *          pa0;
  const xb_vecNx16      *          pa1;
  const xb_vecNx16      * restrict pr0;
  const xb_vecNx16      * restrict pb0;
  const xb_vecNx16      * restrict pb1;
  const xb_vecNx16      * restrict pd;
        xb_vecNx16      * restrict py0;
        xb_vecNx16      * restrict py1;
  vsaN shft,sh16,sh;
  xb_vecNx16 t0;

  int l,_4L;
  // check alignment
  NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(A,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
  NASSERT(L>1);
  NASSERT((L%(BBE_SIMD_WIDTH/2))==0);
  t0 = BBE_MOVVA16(qY-qB);
  shft = BBE_MOVVSV(t0,0);
  t0 = BBE_MOVVINT16(16);
  sh16 = BBE_MOVVSV(t0,0);
  _4L=L*4;
  pd  = (xb_vecNx16   *) D;
  pa0 = (xb_vecNx16   *) A;
  pa1 = (xb_vecNx16   *) (((uintptr_t)A) + _4L);
  pr0 = (xb_vecNx16   *) (((uintptr_t)R) + _4L);
  pb0 = (xb_vecNx16   *) B;
  pb1 = (xb_vecNx16   *) (((uintptr_t)B) + _4L);
  py0 = (xb_vecNx16   *) (y);
  py1 = (xb_vecNx16   *) (((uintptr_t)y) + _4L);
  for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
  {
    xb_vecNx16    a0, a1, yx, b0,b1,r0,d0,d1,t0,t1;
    xb_vecNx40    B10;
    //n==0
    // take A'
    BBE_LVNX16_IP(a0,pa0,2*BBE_SIMD_WIDTH);
    a1 = BBE_LVNX16_X(pa0,(2*2*2*L-2*BBE_SIMD_WIDTH));
    // take R(:,n)
    // calculate A(:,n)'*B-Rn'*Y, 1xP
    BBE_LVNX16_IP(b0,pb0,2*BBE_SIMD_WIDTH);
    //Q(qA+qB)<-Q.qA*Q.qB
    B10=BBE_MULNX16J(b0,a0);
    BBE_LVNX16_IP(b1,pb1,2*BBE_SIMD_WIDTH);
    //Q(qA+qB)<-Q.qA*Q.qB
    BBE_MULANX16J(B10,b1,a1);
    //Q(qA+qY-16)<-Q(qA+qB)+qY-16-qB
    B10 = BBE_SLANX40(B10,shft);
    // k=0..7, n=0
    // y_re=Y[(k+(n)*L)*2+0]==0;
    // y_im=Y[(k+(n)*L)*2+1]==0;
    //C_re+=mul_ls(y_re,r_re)+mul_ls(y_im,r_im)==0; 
    //C_im+=mul_ls(y_im,r_re)-mul_ls(y_re,r_im)==0;

    BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);//Q.(31 - qA - exp)
    BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);//exp
    /*
      B_re=(int32_t)(((int64_t)(A_re)*D[2*k+0])>>16); 
      B_im=(int32_t)(((int64_t)(A_im)*D[2*k+1])>>16); 
      if(sh>0) { B_re<<=sh; B_im<<=sh; }
      else     { B_re>>=-sh;B_im>>=-sh;}
      B_re=(B_re+0x4000)>>15;
      B_im=(B_im+0x4000)>>15;
      B_re=MAX(MIN_INT16,MIN(MAX_INT16,B_re));
      B_im=MAX(MIN_INT16,MIN(MAX_INT16,B_im));
      Y[(k+(n)*L)*2+0]=(int16_t)(B_re);
      Y[(k+(n)*L)*2+1]=(int16_t)(B_im);
      */
#if 0
    t0=BBE_PACKLNX40(B10);
    t1=BBE_PACKVNX40(B10,sh16);
    B10=BBE_MULUUNX16(t0,d0);
    B10=BBE_SRAINX40(B10,16);
    BBE_MULUSANX16(B10,d0,t1);
#elif 1
    {
        xb_vecNx16U h,zero=0;
        t0=BBE_PACKLNX40(B10);
        t1=BBE_PACKVNX40(B10,sh16);
        B10=BBE_MULUSNX16(d0,t1);
        h=BBE_MULUUNX16PACKH(d0,t0);
        BBE_ADDWUANX16(B10,h,zero);
    }
#else
    {
        xb_mvecNx32 M10;
        t0=BBE_MOVSVWXL(B10);
        t1=BBE_MOVSVWXH(B10);
        M10=BBE_MULUSMNX16(d0,t1);
        B10=BBE_MULUUNX16 (d0,t0);
        BBE_SRAIWADDMNX40 (B10,M10,16);
    }
#endif
    sh=BBE_MOVVSV(d1,0);
    B10=BBE_SLSNX40(B10,sh);
    yx=BBE_PACKPNX40(B10);
    BBE_SVNX16_IP(yx,py0,2*BBE_SIMD_WIDTH);
    //n=1
    // take A(:,n)' 
    BBE_LVNX16_IP(a0,pa1,2*BBE_SIMD_WIDTH);
    a1 = BBE_LVNX16_X(pa1,(2*2*2*L-2*BBE_SIMD_WIDTH));
    // calculate A(:,n)'*B-Rn'*Y, 1xP
    //Q(qA+qB)<-Q.qA*Q.qB
    B10=BBE_MULNX16J(b0,a0);
    //Q(qA+qB)<-Q.qA*Q.qB
    BBE_MULANX16J(B10,b1,a1);

    //Q(qA+qY-16)<-Q(qA+qB)+qY-16-qB
    B10 = BBE_SLANX40(B10,shft);

    BBE_LVNX16_IP(r0,pr0,2*BBE_SIMD_WIDTH);
    //Q(qA+qY-16)<-Q(qA+qY-16)+Q.(qY-16)*Q.qA
    BBE_MULSNX16J(B10,yx,r0);

    BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);//Q.(31 - qA - exp)
    BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);//exp
#if 0
    t0=BBE_PACKLNX40(B10);
    t1=BBE_PACKVNX40(B10,sh16);
    B10=BBE_MULUUNX16(t0,d0);
    B10=BBE_SRAINX40(B10,16);
    BBE_MULUSANX16(B10,d0,t1);
#elif 1
    {
        xb_vecNx16U h,zero=0;
        t0=BBE_PACKLNX40(B10);
        t1=BBE_PACKVNX40(B10,sh16);
        B10=BBE_MULUSNX16(d0,t1);
        h=BBE_MULUUNX16PACKH(d0,t0);
        BBE_ADDWUANX16(B10,h,zero);
    }
#else
    {
        xb_mvecNx32 M10;
        t0=BBE_MOVSVWXL(B10);
        t1=BBE_MOVSVWXH(B10);
        M10=BBE_MULUSMNX16(d0,t1);
        B10=BBE_MULUUNX16 (d0,t0);
        BBE_SRAIWADDMNX40 (B10,M10,16);
    }
#endif
    sh=BBE_MOVVSV(d1,0);
    B10=BBE_SLSNX40(B10,sh);
    yx=BBE_PACKPNX40(B10);
    BBE_SVNX16_IP(yx,py1,2*BBE_SIMD_WIDTH);
  }
} /* cholfwd2x2x1s() */
#endif

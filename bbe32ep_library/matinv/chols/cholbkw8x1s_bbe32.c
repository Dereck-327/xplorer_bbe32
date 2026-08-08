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
  NatureDSP_Baseband library. Cholesky backward recursion for pseudo-inversion API (complex data)
    These functions make backward recursion stage of pseudo-inversion. They
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if !(HAVE_VSAMATH && 1)
DISCARD_FUN(void, cholbkw8x1s, (
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L))
#else

#define BKWITERATION(n,x,pd,py,pr,c0,_4L,_28L,inc)                     \
{                                                                      \
    vsaN sh;                                                           \
    xb_vecNx16 t0,y0,d0,d1,r0;                                         \
    xb_vecNx40 A10;                                                    \
    d1 = BBE_LVNX16_I(pd, ((7-n)*2+1)*(2*BBE_SIMD_WIDTH));             \
    d0 = BBE_LVNX16_I(pd, ((7-n)*2+0)*(2*BBE_SIMD_WIDTH));             \
    sh = BBE_MOVVSV(d1,0); sh=BBE_SUBSAVSN(11,sh);                  \
    if(n==7) {BBE_LVNX16_XP(y0,py,_28L);}                                \
    else     {BBE_LVNX16_XP(y0,py,_4L);}                                 \
    A10=BBE_MULUSRNX16(c0,y0,sh);                                      \
    if(n>=1) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##0,r0);};  \
    if(n>=2) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##1,r0);};  \
    if(n>=3) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##2,r0);};  \
    if(n>=4) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##3,r0);};  \
    if(n>=5) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##4,r0);};  \
    if(n>=6) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##5,r0);};  \
    if(n>=7) {BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x##6,r0);};  \
    t0 = BBE_PACKVNX40(A10,sh);                                        \
    A10=BBE_MULUSNX16(d0,t0);                                          \
    x##n = BBE_PACKQNX40(A10);                                         \
    BBE_LVNX16_XP(t0,pr,inc);                                          \
    inc-=_4L;                                                          \
}

/*-------------------------------------------------------------------------
   These functions make backward recursion stage of pseudo-inversion. They
   use Cholesky decomposition of original matrices and results of forward 
   recursion. 
   NOTE:
   Data layout for matrices is selected as for other matrices written in a 
   streaming order. 

   Input:
   N                Matrix dimension (number of columns and rows in 
                    matrices R)
   P                Number of columns in right-side matrices B
   L                Number of matrices
   R[N*N][L]        Cholesky upper triangular matrices R
   D[L/4][N][8]     reciprocal of main diagonal (mantissa, exponent) in the 
                    special format
   y[N*P][L]        Results of forward recursion stage
   qA,qX,qY         fixed point representation of matrices A(or R which is 
                    the same), x and y
   Output:
   x[N*P][L]        Decision matrix x


   Restrictions:
   1. All matrices must not overlap and be aligned on 32-byte boundary 
   2. Number of matrices L must be a multiple of 8
   3. Matrix sizes (M,N) must be greater than 1, P must be >=1
   4. qX+qA-qY must be <=16 
---------------------------------------------------------------------------*/

void cholbkw8x1s (
                  complex_fract16* restrict _x, 
            const complex_fract16* restrict _R,
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _y, 
            int qA, int qY, int qX,
            int L)
{
          int16_t* restrict x=(      int16_t*)_x;
    const int16_t* restrict R=(const int16_t*)_R;
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict y=(const int16_t*)_y;
/*
    Reference code:
    function [X]=cholbkw(R,y, extraBits, isFixedPoint, isShow)
    sz=size(y); N=sz(1);P=sz(2);
    X=zeros(N,P);
    D=real(1./diag(R));
    q=ceil(log2(min(abs(D))));
    q=15-q;
    q=q-extraBits;
    for m=N:-1:1
        Rm=R(m,:); 
        ym=y(m,:);
        x=(ym-Rm*X)*D(m);
        X(m,:)=x;
    end
*/
    //int arr[1];
    vsaN sh;                                                          
    xb_vecNx16 t0,y0,d0,d1,r0;                                        
    xb_vecNx40 A10;                                                   
    int inc;
    int _4L; // -4*L
    int _32L; // -32*L
    int _28L; // 28*L
    vsaN  shft;
    const xb_vecNx16      * restrict pr   = (xb_vecNx16   *) R ;
    const xb_vecNx16      * restrict pd   = (xb_vecNx16   *) D ;
    const xb_vecNx16      * restrict py   = (xb_vecNx16   *) y ;
          xb_vecNx16        * restrict px0  = (xb_vecNx16   *) x ;
    xb_vecNx16    c0,x0,x1,x2,x3,x4,x5,x6,x7;
    int l;

    // check alignment
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);

    shft = BBE_MOVVSA32(qX+qA-qY);
    c0   = BBE_MOVVINT16(1);
    c0   = BBE_SLLNX16(c0,shft);
    
    //D+=32*(8-1);
    y+=2*L*(8-1);
    x+=2*L*(8-1);
    R+=2*(8*8-1)*L;
    pd  = (xb_vecNx16 *) D ;
    l=L>>3;
    _4L=-(L<<2);
    //arr[0]=_4L;
    _32L=_4L<<3;
    px0 = (xb_vecNx16 *) (x);
    py  = (const xb_vecNx16 *) (y);
    inc=_32L;
    _28L=_4L-_32L;
    pr  = (xb_vecNx16 *) (R);
    BKWITERATION(0,x,pd,py,pr,c0,_4L,_28L,inc);
    BKWITERATION(1,x,pd,py,pr,c0,_4L,_28L,inc);
    for(l=(L>>(LOG2_BBE_SIMD_WIDTH-1))-1; l>0; l--)
    {
        BKWITERATION(2,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(3,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(4,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(5,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(6,x,pd,py,pr,c0,_4L,_28L,inc);
        BBE_LVNX16_IP(d0,pd,  2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(d1,pd,15*2*BBE_SIMD_WIDTH);
        sh = BBE_MOVVSV(d1,0); sh=BBE_SUBSAVSN(11,sh);
        BBE_LVNX16_XP(y0,py,_28L);                               
        A10=BBE_MULUSRNX16(c0,y0,sh);                                    
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x0,r0);  
        BBE_SVNX16_XP(x0,px0,_4L); 
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x1,r0);  
        BBE_SVNX16_XP(x1,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x2,r0);  
        BBE_SVNX16_XP(x2,px0,_4L); 
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x3,r0);  
        BBE_SVNX16_XP(x3,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x4,r0);  
        BBE_SVNX16_XP(x4,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x5,r0);  
        BBE_SVNX16_XP(x5,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x6,r0);  
        BBE_SVNX16_XP(x6,px0,_4L);  
        t0 = BBE_PACKVNX40(A10,sh);                                       
        A10=BBE_MULUSNX16(d0,t0);                                    
        x7 = BBE_PACKQNX40(A10);                                   
        BBE_LVNX16_XP(t0,pr,inc);                                         
        BBE_SVNX16_XP(x7,px0,_28L);
        // go to the next 8 matrices
#if defined(__cplusplus) || defined(COMPILER_GNU)
       {
            const xb_vecNx16* temp=(const xb_vecNx16*)R;
            BBE_LVNX16_IP(r0,temp,(2*BBE_SIMD_WIDTH));
            R=(const int16_t*)temp;
       }
#else   // please keep this - the equivalent code for C++ works slower!!!
        BBE_LVNX16_IP(r0,(const xb_vecNx16*)R,(2*BBE_SIMD_WIDTH));
#endif
        BBE_LVNX16_IP(r0,px0,(2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(r0,py ,(2*BBE_SIMD_WIDTH));
        inc-=_28L;                                                         
        pr  = (xb_vecNx16 *) (R);
        BKWITERATION(0,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(1,x,pd,py,pr,c0,_4L,_28L,inc);
    }
    // epilogue
    {
        BKWITERATION(2,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(3,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(4,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(5,x,pd,py,pr,c0,_4L,_28L,inc);
        BKWITERATION(6,x,pd,py,pr,c0,_4L,_28L,inc);
        BBE_LVNX16_IP(d0,pd,  2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(d1,pd,15*2*BBE_SIMD_WIDTH);
        sh = BBE_MOVVSV(d1,0); sh=BBE_SUBSAVSN(11,sh);
        BBE_LVNX16_XP(y0,py,_28L);                               
        A10=BBE_MULUSRNX16(c0,y0,sh);                                    
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x0,r0);  
        BBE_SVNX16_XP(x0,px0,_4L); 
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x1,r0);  
        BBE_SVNX16_XP(x1,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x2,r0);  
        BBE_SVNX16_XP(x2,px0,_4L); 
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x3,r0);  
        BBE_SVNX16_XP(x3,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x4,r0);  
        BBE_SVNX16_XP(x4,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x5,r0);  
        BBE_SVNX16_XP(x5,px0,_4L);  
        BBE_LVNX16_XP(r0,pr ,_4L);BBE_MULSNX16C(A10,x6,r0);  
        BBE_SVNX16_XP(x6,px0,_4L);  
        t0 = BBE_PACKVNX40(A10,sh);                                       
        A10=BBE_MULUSNX16(d0,t0);                                    
        x7 = BBE_PACKQNX40(A10);                                   
        BBE_LVNX16_XP(t0,pr,inc);                                         
        BBE_SVNX16_XP(x7,px0,_28L);
    }
} /* cholbkw8x1s() */
#endif

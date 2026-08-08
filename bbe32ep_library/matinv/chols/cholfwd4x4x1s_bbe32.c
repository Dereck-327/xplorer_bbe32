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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
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

#define FWDITERATION(n,y,pa,b,pd,pr,_4L,_12L,_16L,_44L,inc,shft)                     \
{                                                                                    \
    xb_vecNx16 d0,d1,t0,t1;                                                          \
    xb_vecNx40 B10;                                                                  \
    vsaN sh;                                                                         \
    BBE_LVNX16_XP(a0,pa, _16L); B10=BBE_MULNX16J (b##0,a0);                          \
    BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b##1,a0);                          \
    BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b##2,a0);                          \
    BBE_LVNX16_XP(a0,pa, _44L); BBE_MULANX16J(B10,b##3,a0);                          \
    B10 = BBE_SLANX40(B10,shft);                                                     \
    if(n>=1) {BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y##0,r0);}              \
    if(n>=2) {BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y##1,r0);}              \
    if(n>=3) {BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y##2,r0);}              \
    BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);                                           \
    BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);                                           \
    t0=BBE_PACKLNX40(B10);                                                           \
    t1=BBE_PACKVNX40(B10,sh16);                                                      \
    B10=BBE_MULUSNX16(t0,d0);                                                        \
    B10=BBE_SRAINX40(B10,16);                                                        \
    BBE_MULANX16(B10,d0,t1);                                                         \
    sh=BBE_MOVVSV(d1,0);                                                             \
    B10=BBE_SLLNX40(B10,sh);                                                         \
    y##n=BBE_PACKQNX40(B10);                                                         \
    BBE_LVNX16_XP(t0,pr,inc);                                                        \
    inc-=_16L;                                                                       \
}

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

void cholfwd4x4x1s (
                  complex_fract16* restrict _y,
            const complex_fract16* restrict _R, 
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _A, 
            const complex_fract16* restrict _B, 
            int qA,int qB,int qY,
            int L)
{
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
          int16_t* restrict y=(      int16_t*)_y;
    const int16_t* restrict R=(const int16_t*)_R; 
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict A=(const int16_t*)_A; 
    const int16_t* restrict B=(const int16_t*)_B; 
    int arr[1];
    int _16L; //16*L
    int inc;
    int _4L; //4*L
    int _44L; //-44*L
    int _12L; //-12*L
    vsaN shft,sh16;
    xb_vecNx16 t0;
    xb_vecNx16   a0,b0,b1,b2,b3,r0,y0,y1,y2,y3;
    int l;
    const xb_vecNx16      * restrict pa  = (const xb_vecNx16   *) A ;
    const xb_vecNx16      * restrict pr  = (const xb_vecNx16   *) R ;
    const xb_vecNx16      * restrict pb  = (const xb_vecNx16   *) B ;
    const xb_vecNx16      * restrict pd  = (const xb_vecNx16   *) D ;
          xb_vecNx16      * restrict py  = (xb_vecNx16 *) y ;

    // check alignment
  NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(A,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
  NASSERT(L>1);
  NASSERT((L%(BBE_SIMD_WIDTH/2))==0);

    _4L=L<<2;
    _16L=L<<4;
    _44L=L*(-44);
    _12L=_4L-_16L;
    t0 = BBE_MOVVA16(qY-qB);
    shft = BBE_MOVVSV(t0,0);
    t0 = BBE_MOVVINT16(16);
    sh16 = BBE_MOVVSV(t0,0);

    pd  = (const xb_vecNx16*)D;
    pb  = (const xb_vecNx16*)B;
    py  = (      xb_vecNx16*)y; 
    pa  = (const xb_vecNx16*)A;
    pr   =(const xb_vecNx16*)R;
    l=L>>(LOG2_BBE_SIMD_WIDTH-1);
    arr[0]= _4L;
    do
    {
        xb_vecNx16 d0,d1,t0,t1;      
        xb_vecNx40 B10;      
        vsaN sh;     
        inc=XT_L32I_N(arr,0);
        BBE_LVNX16_XP(b0,pb, _4L);
        BBE_LVNX16_XP(b1,pb, _4L); 
        BBE_LVNX16_XP(b2,pb, _4L); 
        BBE_LVNX16_XP(b3,pb,_12L); 
                                                                        
        BBE_LVNX16_XP(a0,pa, _16L); B10=BBE_MULNX16J (b0,a0);     
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b1,a0);     
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b2,a0);     
        BBE_LVNX16_XP(a0,pa, _44L); BBE_MULANX16J(B10,b3,a0);     
        B10 = BBE_SLANX40(B10,shft);
        BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);       
        BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);       
        t0=BBE_PACKLNX40(B10);      
        t1=BBE_PACKVNX40(B10,sh16); 
        B10=BBE_MULUUNX16(t0,d0);   
        B10=BBE_SRAINX40(B10,16);   
        BBE_MULUSANX16(B10,d0,t1);    
        sh=BBE_MOVVSV(d1,0); 
        B10=BBE_SLSNX40(B10,sh);    
        y0=BBE_PACKPNX40(B10);    
        BBE_LVNX16_XP(t0,pr,inc);   
        inc-=_16L;   

        BBE_LVNX16_XP(a0,pa, _16L); B10=BBE_MULNX16J (b0,a0);
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b1,a0);
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b2,a0);
        BBE_LVNX16_XP(a0,pa, _44L); BBE_MULANX16J(B10,b3,a0);
        B10 = BBE_SLANX40(B10,shft); 
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y0,r0);
        BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);
        t0=BBE_PACKLNX40(B10);      
        t1=BBE_PACKVNX40(B10,sh16); 
        #ifdef COMPILER_XTENSA
        #pragma no_reorder
        #endif
        B10=BBE_MULUUNX16(t0,d0);   
        B10=BBE_SRAINX40(B10,16);   
        BBE_MULUSANX16(B10,d0,t1);    
        sh=BBE_MOVVSV(d1,0); 
        B10=BBE_SLSNX40(B10,sh);    
        y1=BBE_PACKPNX40(B10);    
        BBE_LVNX16_XP(t0,pr,inc);    
        inc-=_16L;    

        BBE_LVNX16_XP(a0,pa, _16L); B10=BBE_MULNX16J (b0,a0);  
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b1,a0);  
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b2,a0);  
        BBE_LVNX16_XP(a0,pa, _44L); BBE_MULANX16J(B10,b3,a0);  
        B10 = BBE_SLANX40(B10,shft);   
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y0,r0);
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y1,r0);
        BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH);  
        BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH);  
        t0=BBE_PACKLNX40(B10);      
        t1=BBE_PACKVNX40(B10,sh16); 
        B10=BBE_MULUUNX16(t0,d0);   
        B10=BBE_SRAINX40(B10,16);   
        BBE_MULUSANX16(B10,d0,t1);    
        sh=BBE_MOVVSV(d1,0); 
        B10=BBE_SLSNX40(B10,sh);    
        y2=BBE_PACKPNX40(B10);    
        BBE_LVNX16_XP(t0,pr,inc);      
        inc-=_16L;      

        BBE_LVNX16_XP(a0,pa, _16L); B10=BBE_MULNX16J (b0,a0); 
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b1,a0); 
        BBE_LVNX16_XP(a0,pa, _16L); BBE_MULANX16J(B10,b2,a0); 
        BBE_LVNX16_XP(a0,pa, _44L); BBE_MULANX16J(B10,b3,a0); 
        B10 = BBE_SLANX40(B10,shft);  
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y0,r0);
        BBE_SVNX16_XP(y0,py,_4L);
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y1,r0);
        BBE_SVNX16_XP(y1,py,_4L);
        BBE_LVNX16_XP(r0,pr ,_16L);  BBE_MULSNX16J(B10,y2,r0);
        BBE_SVNX16_XP(y2,py,_4L);
        BBE_LVNX16_IP(d0,pd,2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(d1,pd,2*BBE_SIMD_WIDTH); 
        t0=BBE_PACKLNX40(B10);      
        t1=BBE_PACKVNX40(B10,sh16); 
        B10=BBE_MULUUNX16(t0,d0);   
        B10=BBE_SRAINX40(B10,16);   
        BBE_MULUSANX16(B10,d0,t1);    
        sh=BBE_MOVVSV(d1,0); 
        B10=BBE_SLSNX40(B10,sh);    
        y3=BBE_PACKPNX40(B10);    
        inc-=_16L;      
        BBE_LVNX16_XP(t0,pr,inc);     

        BBE_SVNX16_XP(y3,py,_12L);
        // go to the next 8 matrices
        pa= (const xb_vecNx16*)(((uintptr_t)pa)-_16L); 
        BBE_LVNX16_IP(r0,pb,(2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(r0,py,(2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(r0,pa,(2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(r0,pr,(2*BBE_SIMD_WIDTH));
    }
    while(--l);
} /* cholfwd4x4x1s() */

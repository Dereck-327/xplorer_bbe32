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
DISCARD_FUN(void, cholbkw4x1s, (
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L))
#else

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

void cholbkw4x1s (
                  complex_fract16* restrict _x, 
            const complex_fract16* restrict _R,
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _y, 
            int qA, int qY, int qX,
            int L)
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
#if 0
{
  int l;
  vsaN  shft,sh;
  xb_vecNx40    A10;
  xb_vecNx16    _11,x0,x1,x2;    
  xb_vecNx16    d0_0, d0_1,c0;
  const xb_vecNx16      * restrict pr0  = (xb_vecNx16   *) R ;
  const xb_vecNx16      * restrict pr1  = (xb_vecNx16   *) R ;
  const xb_vecNx16      * restrict pr2  = (xb_vecNx16   *) R ;
  const xb_vecNx16      * restrict pd   = (xb_vecNx16   *) D ;
  const xb_vecNx16      * restrict py1  = (xb_vecNx16 *) y ;
  const xb_vecNx16      * restrict py3  = (xb_vecNx16 *) y ;
        xb_vecNx16      * restrict px1  = (xb_vecNx16   *) x ;
        xb_vecNx16      * restrict px3  = (xb_vecNx16   *) x ;
    // check alignment
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT (qX+qA-qY<=16);
    x0 = BBE_MOVVA16(qX+qA-qY);
    shft = BBE_MOVVSV(x0,0);
    c0 = BBE_MOVVINT16(1);
    c0 = BBE_SLLNX16(c0,shft);

    _11 = BBE_MOVVINT16(11);   
    py1 = (xb_vecNx16 *) (y+2*2*L);
    py3 = (xb_vecNx16 *) (y+0*2*L);
    px1 = (xb_vecNx16   *) (x+2*2*L);
    px3 = (xb_vecNx16   *) (x+0*2*L);
    pr0 = (xb_vecNx16   *) (R+ 22*L) ;
    pr1 = (xb_vecNx16   *) (R+ 12*L) ;
    pr2 = (xb_vecNx16   *) (R+  2*L) ;
    pd  = (xb_vecNx16   *) D+7;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 t0,y0,r0,r1,r2;
   
        //m=3
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        d0_1 = BBE_SUBNX16(_11,d0_1);
        y0=BBE_LVNX16_X(py1,2*2*L);
        sh=BBE_MOVVSV(d0_1,0);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x0 = BBE_PACKQNX40(A10);
        BBE_SVNX16_X(x0,px1,2*2*L); 

        //m=2
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        d0_1 = BBE_SUBNX16(_11,d0_1);
        sh=BBE_MOVVSV(d0_1,0);
        BBE_LVNX16_IP(y0,py1,2*BBE_SIMD_WIDTH);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        BBE_LVNX16_IP(r0,pr0,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x1 = BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(x1,px1,2*BBE_SIMD_WIDTH);     

        //m=1
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        d0_1 = BBE_SUBNX16(_11,d0_1);
        sh=BBE_MOVVSV(d0_1,0);
        y0=BBE_LVNX16_X(py3,2*2*L);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        r1 = BBE_LVNX16_X(pr1,2*2*L);
        BBE_LVNX16_IP(r0,pr1,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r1);
        BBE_MULSNX16C(A10,x1,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x2 = BBE_PACKQNX40(A10);
        BBE_SVNX16_X(x2,px3,2*2*L);   

       // m=0
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_XP(d0_0,pd, 15*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        d0_1 = BBE_SUBNX16(_11,d0_1);
        sh=BBE_MOVVSV(d0_1,0);
        BBE_LVNX16_IP(y0,py3,2*BBE_SIMD_WIDTH);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        r2 = BBE_LVNX16_X(pr2,4*2*L);
        r1 = BBE_LVNX16_X(pr2,2*2*L);
        BBE_LVNX16_IP(r0,pr2,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r2);
        BBE_MULSNX16C(A10,x1,r1);
        BBE_MULSNX16C(A10,x2,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x0 = BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(x0,px3,2*BBE_SIMD_WIDTH);   
    }
} /* cholbkw4x1s() */
#else
// -----------------------------------------------------------------------
// the best variant ....
// -----------------------------------------------------------------------
{
          int16_t* restrict x=(      int16_t*)_x;
    const int16_t* restrict R=(const int16_t*)_R;
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict y=(const int16_t*)_y;
  int l;
  vsaN  shft,sh;
  xb_vecNx40    A10;
  xb_vecNx16    x0,x1,x2,x3;    
  xb_vecNx16 t0,y0,r0,r1,r2;
  xb_vecNx16    d0_0, d0_1,c0;
  const xb_vecNx16      * restrict pr0  = (const xb_vecNx16 *) R ;
  const xb_vecNx16      * restrict pr1  = (const xb_vecNx16 *) R ;
  const xb_vecNx16      * restrict pr2  = (const xb_vecNx16 *) R ;
  const xb_vecNx16      * restrict pd   = (const xb_vecNx16 *) D ;
  const xb_vecNx16      * restrict py1  = (const xb_vecNx16 *) y ;
  const xb_vecNx16      * restrict py3  = (const xb_vecNx16 *) y ;
        xb_vecNx16      * restrict px1  = (xb_vecNx16   *) x ;
        xb_vecNx16      * restrict px3  = (xb_vecNx16   *) x ;
    // check alignment
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT (qX+qA-qY<=16);
    x0 = BBE_MOVVA16(qX+qA-qY);
    shft = BBE_MOVVSV(x0,0);
    c0 = BBE_MOVVINT16(1);
    c0 = BBE_SLLNX16(c0,shft);

    py1 = (const xb_vecNx16 *) (y+2*2*L);
    py3 = (const xb_vecNx16 *) (y+0*2*L);
    px1 = (      xb_vecNx16 *) (x+2*2*L);
    px3 = (      xb_vecNx16 *) (x+0*2*L);
    pr0 = (const xb_vecNx16 *) (R+ 22*L) ;
    pr1 = (const xb_vecNx16 *) (R+ 12*L) ;
    pr2 = (const xb_vecNx16 *) (R+  2*L) ;
    pd  = (const xb_vecNx16 *) D+7;

    // prologue of manual sw pipelining
    BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
    BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
    sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
    y0=BBE_LVNX16_X(py1,2*2*L);
    A10=BBE_MULUSRNX16(c0,y0,sh);
    t0 = BBE_PACKVNX40(A10,sh);
    A10=BBE_MULUSNX16(d0_0,t0);
    x0 = BBE_PACKQNX40(A10);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1))-1; l++)
    {
        //m=3
        BBE_SVNX16_X(x0,px1,2*2*L); 

        //m=2
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
        BBE_LVNX16_IP(y0,py1,2*BBE_SIMD_WIDTH);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        BBE_LVNX16_IP(r0,pr0,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x1 = BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(x1,px1,2*BBE_SIMD_WIDTH);     

        //m=1
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
        y0=BBE_LVNX16_X(py3,2*2*L);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        r1 = BBE_LVNX16_X(pr1,2*2*L);
        BBE_LVNX16_IP(r0,pr1,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r1);
        BBE_MULSNX16C(A10,x1,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x2 = BBE_PACKQNX40(A10);
        BBE_SVNX16_X(x2,px3,2*2*L);   

       // m=0
        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_XP(d0_0,pd, 15*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
        BBE_LVNX16_IP(y0,py3,2*BBE_SIMD_WIDTH);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        r2 = BBE_LVNX16_X(pr2,4*2*L);
        r1 = BBE_LVNX16_X(pr2,2*2*L);
        BBE_LVNX16_IP(r0,pr2,2*BBE_SIMD_WIDTH);
        BBE_MULSNX16C(A10,x0,r2);
        BBE_MULSNX16C(A10,x1,r1);
        BBE_MULSNX16C(A10,x2,r0);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x3 = BBE_PACKQNX40(A10);

        BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
        BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
        sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
        y0=BBE_LVNX16_X(py1,2*2*L);
        A10=BBE_MULUSRNX16(c0,y0,sh);
        t0 = BBE_PACKVNX40(A10,sh);
        A10=BBE_MULUSNX16(d0_0,t0);
        x0 = BBE_PACKQNX40(A10);

        BBE_SVNX16_IP(x3,px3,2*BBE_SIMD_WIDTH);   
    }
    // epilogue of manual sw pipelining
    //m=3
    BBE_SVNX16_X(x0,px1,2*2*L); 
   //m=2
    BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
    BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
    sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
    BBE_LVNX16_IP(y0,py1,2*BBE_SIMD_WIDTH);
    A10=BBE_MULUSRNX16(c0,y0,sh);
    BBE_LVNX16_IP(r0,pr0,2*BBE_SIMD_WIDTH);
    BBE_MULSNX16C(A10,x0,r0);
    t0 = BBE_PACKVNX40(A10,sh);
    A10=BBE_MULUSNX16(d0_0,t0);
    x1 = BBE_PACKQNX40(A10);
    BBE_SVNX16_IP(x1,px1,2*BBE_SIMD_WIDTH);     

    //m=1
    BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
    BBE_LVNX16_IP(d0_0,pd, -1*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
    sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
    y0=BBE_LVNX16_X(py3,2*2*L);
    A10=BBE_MULUSRNX16(c0,y0,sh);
    r1 = BBE_LVNX16_X(pr1,2*2*L);
    BBE_LVNX16_IP(r0,pr1,2*BBE_SIMD_WIDTH);
    BBE_MULSNX16C(A10,x0,r1);
    BBE_MULSNX16C(A10,x1,r0);
    t0 = BBE_PACKVNX40(A10,sh);
    A10=BBE_MULUSNX16(d0_0,t0);
    x2 = BBE_PACKQNX40(A10);
    BBE_SVNX16_X(x2,px3,2*2*L);   

   // m=0
    BBE_LVNX16_IP(d0_1,pd, -1*2*BBE_SIMD_WIDTH);// exp
    BBE_LVNX16_XP(d0_0,pd, 15*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
    sh=BBE_MOVVSV(d0_1,0); sh=BBE_SUBSAVSN(11,sh);
    BBE_LVNX16_IP(y0,py3,2*BBE_SIMD_WIDTH);
    A10=BBE_MULUSRNX16(c0,y0,sh);
    r2 = BBE_LVNX16_X(pr2,4*2*L);
    r1 = BBE_LVNX16_X(pr2,2*2*L);
    BBE_LVNX16_IP(r0,pr2,2*BBE_SIMD_WIDTH);
    BBE_MULSNX16C(A10,x0,r2);
    BBE_MULSNX16C(A10,x1,r1);
    BBE_MULSNX16C(A10,x2,r0);
    t0 = BBE_PACKVNX40(A10,sh);
    A10=BBE_MULUSNX16(d0_0,t0);
    x3 = BBE_PACKQNX40(A10);
    BBE_SVNX16_IP(x3,px3,2*BBE_SIMD_WIDTH);   
} /* cholbkw4x1s() */
#endif

#endif

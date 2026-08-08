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
DISCARD_FUN(void, cholbkw2x1s, (
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

void cholbkw2x1s (
                  complex_fract16* restrict _x, 
            const complex_fract16* restrict _R,
            const complex_fract16* restrict _D,
            const complex_fract16* restrict _y, 
            int qA, int qY, int qX,
            int L)
{
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
          int16_t* restrict x=(      int16_t*)_x;
    const int16_t* restrict R=(const int16_t*)_R;
    const int16_t* restrict D=(const int16_t*)_D;
    const int16_t* restrict y=(const int16_t*)_y;
    int l,_4L;
    vsaN  shft,sh;
    xb_vecNx40    A10;
    xb_vecNx16    d0_0, d0_1, d1_0, d1_1;  
    xb_vecNx16    r0, x0, x1, y0, y1, t0, c0;
    const xb_vecNx16      * restrict pr;
    const xb_vecNx16      * restrict pd;
    const xb_vecNx16      * restrict py;
          xb_vecNx16      * restrict pxw;
    const xb_vecNx16      * restrict pxr; // read-only pointer for x

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
    _4L=L<<2;

    if (L<=BBE_SIMD_WIDTH*2)
    {   
        /*--------------------------------------
        code for smaller L
        --------------------------------------*/
        py  = (const xb_vecNx16   *) y;
        pxw = (      xb_vecNx16   *) x;
        pr  = (const xb_vecNx16   *) (((uintptr_t)R)+_4L);
        pd  = (const xb_vecNx16   *) D;
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            BBE_LVNX16_IP(d1_0,pd,2*BBE_SIMD_WIDTH);// exp
            BBE_LVNX16_IP(d1_1,pd,2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
            BBE_LVNX16_IP(d0_0,pd,2*BBE_SIMD_WIDTH);// exp
            BBE_LVNX16_IP(d0_1,pd,2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
            y1=BBE_LVNX16_X( py ,_4L);
            BBE_LVNX16_IP(y0,py ,2*BBE_SIMD_WIDTH);
            sh=BBE_MOVVSV(d0_1,0);
            sh=BBE_SUBSAVSN(11,sh);
            A10=BBE_MULUSRNX16(c0,y1,sh);
            t0 = BBE_PACKVNX40(A10,sh);
            A10=BBE_MULUSNX16(d0_0,t0);
            x1 = BBE_PACKQNX40(A10);
            sh=BBE_MOVVSV(d1_1,0);
            sh=BBE_SUBSAVSN(11,sh);
            A10=BBE_MULUSRNX16(c0,y0,sh);
            BBE_LVNX16_XP(r0,pr,2*BBE_SIMD_WIDTH);
            //Q(qX+qA)<=Q.qX * Q.qA
            BBE_MULSNX16C(A10,x1,r0);
            t0 = BBE_PACKVNX40(A10,sh);
            A10=BBE_MULUSNX16(d1_0,t0);
            x0 = BBE_PACKQNX40(A10);
            BBE_SVNX16_X (x1,pxw,_4L); 
            BBE_SVNX16_IP(x0,pxw,2*BBE_SIMD_WIDTH); 
        }
    }
    else
    {
        /*--------------------------------------
        code for bigger L
        --------------------------------------*/
        pxr =
        pxw = (      xb_vecNx16   *) (((uintptr_t)x)+_4L);
        py  = (      xb_vecNx16   *) (((uintptr_t)y)+_4L);
        pr  = (const xb_vecNx16   *) (((uintptr_t)R)+_4L);
        pd  = (const xb_vecNx16   *) (((uintptr_t)D)+4*BBE_SIMD_WIDTH);
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            BBE_LVNX16_IP(d0_0,pd,  2*BBE_SIMD_WIDTH);// exp
            BBE_LVNX16_IP(d0_1,pd,3*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
            BBE_LVNX16_IP(y1,py ,2*BBE_SIMD_WIDTH);
            sh=BBE_MOVVSV(d0_1,0);
            sh=BBE_SUBSAVSN(11,sh);
            A10=BBE_MULUSRNX16(c0,y1,sh);
            t0 = BBE_PACKVNX40(A10,sh);
            A10=BBE_MULUSNX16(d0_0,t0);
            x1 = BBE_PACKQNX40(A10);
            BBE_SVNX16_IP(x1,pxw,2*BBE_SIMD_WIDTH); 
        }
        pxw = (      xb_vecNx16   *) (((uintptr_t)pxr)-_4L);
        py  = (      xb_vecNx16   *) y;
        pd  = (const xb_vecNx16   *) D;
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            BBE_LVNX16_IP(d1_0,pd,  2*BBE_SIMD_WIDTH);// exp
            BBE_LVNX16_IP(d1_1,pd,3*2*BBE_SIMD_WIDTH);// Q.(31 - qA - exp)
            BBE_LVNX16_IP(y0,py ,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(x1,pxr,2*BBE_SIMD_WIDTH);
            sh=BBE_MOVVSV(d1_1,0);
            sh=BBE_SUBSAVSN(11,sh);
            A10=BBE_MULUSRNX16(c0,y0,sh);
            BBE_LVNX16_XP(r0,pr,2*BBE_SIMD_WIDTH);
            //Q(qX+qA)<=Q.qX * Q.qA
            BBE_MULSNX16C(A10,x1,r0);
            t0 = BBE_PACKVNX40(A10,sh);
            A10=BBE_MULUSNX16(d1_0,t0);
            x0 = BBE_PACKQNX40(A10);
            BBE_SVNX16_IP(x0,pxw,2*BBE_SIMD_WIDTH); 
        }
    }
} /* cholbkw2x1s() */
#endif

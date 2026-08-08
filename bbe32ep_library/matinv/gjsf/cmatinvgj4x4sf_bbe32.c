/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
    Complex Matrix Gauss-Jordan inversion, floating point complex data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include <math.h>
#include <float.h>
#include "cmatinvgjnxnsf_common.h"
/*
    Reference Matlab code:
    % Gauss-Jordan inversion with full pivoting
    function A=gj_inversion1(A)
    sz=size(A);
    N=sz(1);
    swap_row=[];
    swap_col=[];
    for k=1:N
        % pivoting
         T=abs(A);
        T(swap_col,:)=0;
        T(:,swap_col)=0;
        [amax,imax]=max(reshape(T,1,N*N));
        pivot_col=floor((imax-1)/N)+1;
        pivot_row=rem(imax-1,N)+1;
        % swap rows
        t=A(pivot_row,:); A(pivot_row,:)=A(pivot_col,:); A(pivot_col,:)=t;
        swap_row(k)=pivot_row;
        swap_col(k)=pivot_col;
        % process pivot row
        norm=1./A(pivot_col,pivot_col);
        A(pivot_col,pivot_col)=1;
        A(pivot_col,:)=A(pivot_col,:)*norm;
        % elimination
        for i=1:N
            if(i==pivot_col) continue; end
            t=A(i,pivot_col);
            A(i,pivot_col)=0;
            A(i,:)=A(i,:)-A(pivot_col,:)*t;
        end
    end
    % final permutation of columns
    for k = N:-1:1
        t = A(:, swap_row(k));
        A(:, swap_row(k))=A(:, swap_col(k));
        A(:, swap_col(k))=t;
    end

*/
#if HAVE_VFPU

#include <complex.h>
#if 0
static complex_float makecomplexf(float32_t re,float32_t im)
{
    union { complex_float c; struct { float32_t re,im;}s; } w;
    w.s.re=re;
    w.s.im=im;
    return w.c;
}

static complex_float subc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)-crealf(y);
    z.s.im=cimagf(x)-cimagf(y);
    return z.u;
}

static complex_float mulc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)*crealf(y) - cimagf(x)*cimagf(y);
    z.s.im=crealf(x)*cimagf(y) + cimagf(x)*crealf(y);
    return z.u;
}

static complex_float recipc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    float32_t d;
    d=crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
    d=1.0f/d;
    z.s.re= crealf(x)*d;
    z.s.im=-cimagf(x)*d;
    return z.u;
}
static float32_t sqrc(complex_float x)
{
    return crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
}
#endif

/*-------------------------------------------------------------------------
Inversion of Stream Ordered Matrices By Gauss-Jordan Algortihm

Description: perform in-place inversion of real/complex matrices by Gauss-
Jordan elimination method, with full pivoting. The algorithm is applied to
a sequence of input matrices stored in stream order.

Inversion result is not defined for a close to singular input matrix.

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [c]matinvgj<size>sf_getScratchSize(N,L)
Input:
  N          Matrix size
  L          Number of matrices
Input/Output:
  A[N*N][L]  Input matrices, inverted matrices on output
Restrictions:
  pScr,A     Must not overlap and must be aligned on 32-byte boundary 
  N          Must be greater than 1
  L          Must be a multiple of 8 for real-valued functions, or a multiple
             of 4 for complex-valued functions
---------------------------------------------------------------------------*/
void cmatinvgj4x4sf ( void* pScr, complex_float * restrict z, int L  )
#if 0
{
    int32_t *srow,*scol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);

    if(L<=0 ) return;

    srow=(int32_t*)(pScr);
    scol=srow+4*L;

    /*------------*/
    /* stage 1    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {   
        float32_t t,maxval;
        maxval=sqrc(z[l+L* 0]);          {           srow[0*L+l]=0; scol[0*L+l]=0; }
        t=sqrc(z[l+L* 1]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=1; }
        t=sqrc(z[l+L* 2]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=2; }
        t=sqrc(z[l+L* 3]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=3; }

        t=sqrc(z[l+L* 4]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=0; }
        t=sqrc(z[l+L* 5]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=1; }
        t=sqrc(z[l+L* 6]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=2; }
        t=sqrc(z[l+L* 7]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=3; }

        t=sqrc(z[l+L* 8]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=0; }
        t=sqrc(z[l+L* 9]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=1; }
        t=sqrc(z[l+L*10]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=2; }
        t=sqrc(z[l+L*11]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=3; }

        t=sqrc(z[l+L*12]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=0; }
        t=sqrc(z[l+L*13]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=1; }
        t=sqrc(z[l+L*14]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=2; }
        t=sqrc(z[l+L*15]);  if(t>maxval) {           srow[0*L+l]=3; scol[0*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        complex_float t;
        int prow=srow[0*L+l];
        int pcol=scol[0*L+l];
        /* swap rows / columns*/
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 0]; z[l+L* 0]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 1]; z[l+L* 1]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L* 2]; z[l+L* 2]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L* 3]; z[l+L* 3]=t;
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 0]; z[l+L* 0]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 4]; z[l+L* 4]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L* 8]; z[l+L* 8]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*12]; z[l+L*12]=t;
    }

    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L* 0]=recipc(z[l+L* 0]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L* 1] = mulc(z[l+L* 1] , z[l+L* 0]);
        z[l+L* 2] = mulc(z[l+L* 2] , z[l+L* 0]);
        z[l+L* 3] = mulc(z[l+L* 3] , z[l+L* 0]);
        z[l+L* 5] = subc(z[l+L* 5], mulc(z[l+L* 1] , z[l+L* 4])); 
        z[l+L* 6] = subc(z[l+L* 6], mulc(z[l+L* 2] , z[l+L* 4])); 
        z[l+L* 7] = subc(z[l+L* 7], mulc(z[l+L* 3] , z[l+L* 4])); 
        z[l+L* 9] = subc(z[l+L* 9], mulc(z[l+L* 1] , z[l+L* 8])); 
        z[l+L*10] = subc(z[l+L*10], mulc(z[l+L* 2] , z[l+L* 8])); 
        z[l+L*11] = subc(z[l+L*11], mulc(z[l+L* 3] , z[l+L* 8])); 
        z[l+L*13] = subc(z[l+L*13], mulc(z[l+L* 1] , z[l+L*12])); 
        z[l+L*14] = subc(z[l+L*14], mulc(z[l+L* 2] , z[l+L*12])); 
        z[l+L*15] = subc(z[l+L*15], mulc(z[l+L* 3] , z[l+L*12])); 
        z[l+L* 4] = subc(makecomplexf(0.f,0.f), mulc(z[l+L* 0] , z[l+L* 4])); 
        z[l+L* 8] = subc(makecomplexf(0.f,0.f), mulc(z[l+L* 0] , z[l+L* 8])); 
        z[l+L*12] = subc(makecomplexf(0.f,0.f), mulc(z[l+L* 0] , z[l+L*12])); 
    }
    /*------------*/
    /* stage 2    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {   
        float32_t t,maxval;
        maxval=sqrc(z[l+L*(1*4+1)]);          {           srow[1*L+l]=1; scol[1*L+l]=1; }
        t=sqrc(z[l+L*(1*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=1; scol[1*L+l]=2; }
        t=sqrc(z[l+L*(1*4+3)]);  if(t>maxval) { maxval=t; srow[1*L+l]=1; scol[1*L+l]=3; }

        t=sqrc(z[l+L*(2*4+1)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=1; }
        t=sqrc(z[l+L*(2*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=2; }
        t=sqrc(z[l+L*(2*4+3)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=3; }

        t=sqrc(z[l+L*(3*4+1)]);  if(t>maxval) { maxval=t; srow[1*L+l]=3; scol[1*L+l]=1; }
        t=sqrc(z[l+L*(3*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=3; scol[1*L+l]=2; }
        t=sqrc(z[l+L*(3*4+3)]);  if(t>maxval) {           srow[1*L+l]=3; scol[1*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        complex_float t;
        int prow=srow[1*L+l];
        int pcol=scol[1*L+l];
        /* swap rows / columns*/
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 4]; z[l+L* 4]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 5]; z[l+L* 5]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L* 6]; z[l+L* 6]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L* 7]; z[l+L* 7]=t;
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 1]; z[l+L* 1]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 5]; z[l+L* 5]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L* 9]; z[l+L* 9]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*13]; z[l+L*13]=t;
    }

    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*5]=recipc(z[l+L*5]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L* 4] = mulc(z[l+L *4] , z[l+L* 5]);
        z[l+L* 6] = mulc(z[l+L *6] , z[l+L* 5]);
        z[l+L* 7] = mulc(z[l+L *7] , z[l+L* 5]);
        z[l+L* 0] = subc(z[l+L* 0],mulc(z[l+L* 4] , z[l+L* 1])); 
        z[l+L* 2] = subc(z[l+L* 2],mulc(z[l+L* 6] , z[l+L* 1])); 
        z[l+L* 3] = subc(z[l+L* 3],mulc(z[l+L* 7] , z[l+L* 1])); 
        z[l+L* 8] = subc(z[l+L* 8],mulc(z[l+L* 4] , z[l+L* 9])); 
        z[l+L*10] = subc(z[l+L*10],mulc(z[l+L* 6] , z[l+L* 9])); 
        z[l+L*11] = subc(z[l+L*11],mulc(z[l+L* 7] , z[l+L* 9])); 
        z[l+L*12] = subc(z[l+L*12],mulc(z[l+L* 4] , z[l+L*13])); 
        z[l+L*14] = subc(z[l+L*14],mulc(z[l+L* 6] , z[l+L*13])); 
        z[l+L*15] = subc(z[l+L*15],mulc(z[l+L* 7] , z[l+L*13])); 
        z[l+L* 1] = subc(makecomplexf(0.f,0.f),mulc(z[l+L* 5] , z[l+L* 1])); 
        z[l+L* 9] = subc(makecomplexf(0.f,0.f),mulc(z[l+L* 5] , z[l+L* 9])); 
        z[l+L*13] = subc(makecomplexf(0.f,0.f),mulc(z[l+L* 5] , z[l+L*13])); 
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {   
        float32_t t,maxval;
        maxval=sqrc(z[l+L*(2*4+2)]);          {           srow[2*L+l]=2; scol[2*L+l]=2; }
        t=sqrc(z[l+L*(2*4+3)]);  if(t>maxval) { maxval=t; srow[2*L+l]=2; scol[2*L+l]=3; }
        t=sqrc(z[l+L*(3*4+2)]);  if(t>maxval) { maxval=t; srow[2*L+l]=3; scol[2*L+l]=2; }
        t=sqrc(z[l+L*(3*4+3)]);  if(t>maxval) {           srow[2*L+l]=3; scol[2*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        complex_float t;
        int prow=srow[2*L+l];
        int pcol=scol[2*L+l];
        /* swap rows / columns*/
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 8]; z[l+L* 8]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 9]; z[l+L* 9]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L*10]; z[l+L*10]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L*11]; z[l+L*11]=t;
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 2]; z[l+L* 2]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 6]; z[l+L* 6]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L*10]; z[l+L*10]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*14]; z[l+L*14]=t;
    }

    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*10]=recipc(z[l+L*10]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L* 8] = mulc(z[l+L* 8] , z[l+L*10]);
        z[l+L* 9] = mulc(z[l+L* 9] , z[l+L*10]);
        z[l+L*11] = mulc(z[l+L*11] , z[l+L*10]);
        z[l+L* 0] = subc(z[l+L* 0] , mulc(z[l+L* 8] , z[l+L* 2])); 
        z[l+L* 1] = subc(z[l+L* 1] , mulc(z[l+L* 9] , z[l+L* 2])); 
        z[l+L* 3] = subc(z[l+L* 3] , mulc(z[l+L*11] , z[l+L* 2])); 
        z[l+L* 4] = subc(z[l+L* 4] , mulc(z[l+L* 8] , z[l+L* 6])); 
        z[l+L* 5] = subc(z[l+L* 5] , mulc(z[l+L* 9] , z[l+L* 6])); 
        z[l+L* 7] = subc(z[l+L* 7] , mulc(z[l+L*11] , z[l+L* 6])); 
        z[l+L*12] = subc(z[l+L*12] , mulc(z[l+L* 8] , z[l+L*14])); 
        z[l+L*13] = subc(z[l+L*13] , mulc(z[l+L* 9] , z[l+L*14])); 
        z[l+L*15] = subc(z[l+L*15] , mulc(z[l+L*11] , z[l+L*14])); 
        z[l+L* 2] =  subc(makecomplexf(0.f,0.f),mulc(z[l+L*10] , z[l+L* 2])); 
        z[l+L* 6] =  subc(makecomplexf(0.f,0.f),mulc(z[l+L*10] , z[l+L* 6])); 
        z[l+L*14] =  subc(makecomplexf(0.f,0.f),mulc(z[l+L*10] , z[l+L*14])); 
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*15]=recipc(z[l+L*15]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L*12] = mulc(z[l+L*12] , z[l+L*15]);
        z[l+L*13] = mulc(z[l+L*13] , z[l+L*15]);
        z[l+L*14] = mulc(z[l+L*14] , z[l+L*15]);
        z[l+L* 0] = subc(z[l+L* 0],mulc(z[l+L*12] , z[l+L* 3])); 
        z[l+L* 1] = subc(z[l+L* 1],mulc(z[l+L*13] , z[l+L* 3])); 
        z[l+L* 2] = subc(z[l+L* 2],mulc(z[l+L*14] , z[l+L* 3])); 
        z[l+L* 4] = subc(z[l+L* 4],mulc(z[l+L*12] , z[l+L* 7])); 
        z[l+L* 5] = subc(z[l+L* 5],mulc(z[l+L*13] , z[l+L* 7])); 
        z[l+L* 6] = subc(z[l+L* 6],mulc(z[l+L*14] , z[l+L* 7])); 
        z[l+L* 8] = subc(z[l+L* 8],mulc(z[l+L*12] , z[l+L*11])); 
        z[l+L* 9] = subc(z[l+L* 9],mulc(z[l+L*13] , z[l+L*11])); 
        z[l+L*10] = subc(z[l+L*10],mulc(z[l+L*14] , z[l+L*11])); 
        z[l+L* 3] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*15] , z[l+L* 3])); 
        z[l+L* 7] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*15] , z[l+L* 7])); 
        z[l+L*11] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*15] , z[l+L*11])); 
    }
    /* final reverse permulation of columns  */
    for (l=0; l<L; l++)
    {
        complex_float t;
        int prow,pcol;
        /* swap rows / columns*/
        prow=scol[2*L+l];
        pcol=srow[2*L+l];
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 2]; z[l+L* 2]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 6]; z[l+L* 6]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L*10]; z[l+L*10]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*14]; z[l+L*14]=t;
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 8]; z[l+L* 8]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 9]; z[l+L* 9]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L*10]; z[l+L*10]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L*11]; z[l+L*11]=t;

        /* swap rows / columns*/
        prow=scol[1*L+l];
        pcol=srow[1*L+l];
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 1]; z[l+L* 1]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 5]; z[l+L* 5]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L* 9]; z[l+L* 9]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*13]; z[l+L*13]=t;
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 4]; z[l+L* 4]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 5]; z[l+L* 5]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L* 6]; z[l+L* 6]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L* 7]; z[l+L* 7]=t;

        prow=scol[0*L+l];
        pcol=srow[0*L+l];
        t=z[l+L*( 0+pcol)]; z[l+L*( 0+pcol)]=z[l+L* 0]; z[l+L* 0]=t;
        t=z[l+L*( 4+pcol)]; z[l+L*( 4+pcol)]=z[l+L* 4]; z[l+L* 4]=t;
        t=z[l+L*( 8+pcol)]; z[l+L*( 8+pcol)]=z[l+L* 8]; z[l+L* 8]=t;
        t=z[l+L*(12+pcol)]; z[l+L*(12+pcol)]=z[l+L*12]; z[l+L*12]=t;
        t=z[l+L*(prow*4+0)]; z[l+L*(prow*4+0)]=z[l+L* 0]; z[l+L* 0]=t;
        t=z[l+L*(prow*4+1)]; z[l+L*(prow*4+1)]=z[l+L* 1]; z[l+L* 1]=t;
        t=z[l+L*(prow*4+2)]; z[l+L*(prow*4+2)]=z[l+L* 2]; z[l+L* 2]=t;
        t=z[l+L*(prow*4+3)]; z[l+L*(prow*4+3)]=z[l+L* 3]; z[l+L* 3]=t;
    }
}
#else
{
    xb_vecNx16 *pRowCol;
    const xb_vecN_4xcf32 * restrict pZrd;
          xb_vecN_4xcf32 * restrict pZwr;
    int32_t *srow;
    vboolN_2 * browcol;  /* pairs of row/col boolean flags for permutations, [12*L>>((LOG2_BBE_SIMD_WIDTH-2))/]*/
    vboolN_2 * restrict bRowCol0;
    vboolN_2 * restrict bRowCol1;
    vboolN_2 * restrict bRowCol2;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);

    if(L<=0 ) return;

    srow=(int32_t*)(pScr);
    browcol=(vboolN_2*)(srow+4*L);

    /*------------*/
    /* stage 1    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used 
        NOTE: search loop is separated by 2 parts for better scheduling 
    */
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pRowCol=(xb_vecNx16*)srow;
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 t,maxval;
        xb_vecN_4xcf32 a,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11;
        vboolN_2 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11;
        xb_vecN_2xc16 xcol=0,xrow=0;
        BBE_LVN_4XCF32_XP(z0 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z5 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z10,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,2*BBE_SIMD_WIDTH-11*L*sizeof(complex_float));

        a=BBE_MULMN_4XCF32(z0 ,z0 ,0,0); BBE_MULMASN_4XCF32(a,z0 ,z0 ,0,15); maxval=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        a=BBE_MULMN_4XCF32(z1 ,z1 ,0,0); BBE_MULMASN_4XCF32(a,z1 ,z1 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b1 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b1 ); 
        a=BBE_MULMN_4XCF32(z2 ,z2 ,0,0); BBE_MULMASN_4XCF32(a,z2 ,z2 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b2 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b2 ); 
        a=BBE_MULMN_4XCF32(z3 ,z3 ,0,0); BBE_MULMASN_4XCF32(a,z3 ,z3 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b3 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b3 ); 
        a=BBE_MULMN_4XCF32(z4 ,z4 ,0,0); BBE_MULMASN_4XCF32(a,z4 ,z4 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b4 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b4 ); 
        a=BBE_MULMN_4XCF32(z5 ,z5 ,0,0); BBE_MULMASN_4XCF32(a,z5 ,z5 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b5 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b5 ); 
        a=BBE_MULMN_4XCF32(z6 ,z6 ,0,0); BBE_MULMASN_4XCF32(a,z6 ,z6 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b6 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b6 ); 
        a=BBE_MULMN_4XCF32(z7 ,z7 ,0,0); BBE_MULMASN_4XCF32(a,z7 ,z7 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b7 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b7 ); xrow=BBE_MOVN_2XC16T(1,xrow,b4|b5|b6|b7 ); 
        a=BBE_MULMN_4XCF32(z8 ,z8 ,0,0); BBE_MULMASN_4XCF32(a,z8 ,z8 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b8 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b8 ); 
        a=BBE_MULMN_4XCF32(z9 ,z9 ,0,0); BBE_MULMASN_4XCF32(a,z9 ,z9 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b9 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b9 ); 
        a=BBE_MULMN_4XCF32(z10,z10,0,0); BBE_MULMASN_4XCF32(a,z10,z10,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b10=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b10); 
        a=BBE_MULMN_4XCF32(z11,z11,0,0); BBE_MULMASN_4XCF32(a,z11,z11,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b11); xrow=BBE_MOVN_2XC16T(2,xrow,b8|b9|b10|b11); 

        BBE_SVNX16_IP(BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(xcol,xrow,BBE_SELI_INTERLEAVE_2_EVEN)),pRowCol,2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(BBE_MOVNX16_FROMN_2XF32(maxval),pRowCol,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+12*L);
    pRowCol=(xb_vecNx16*)srow;
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecNx16 tmp;
        xb_vecN_2xf32 t,maxval;
        xb_vecN_4xcf32 a,z12,z13,z14,z15;
        vboolN_2 b12,b13,b14,b15;
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xc16 xcol=0,xrow=0;
        BBE_LVN_4XCF32_XP(z12,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));

        BBE_LVNX16_IP(tmp,pRowCol,2*BBE_SIMD_WIDTH);
        xrow=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(tmp,tmp,BBE_SELI_INTERLEAVE_2_EVEN));
        xcol=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(tmp,tmp,BBE_SELI_INTERLEAVE_2_ODD));
        BBE_LVNX16_IP(tmp,pRowCol,2*BBE_SIMD_WIDTH);
        maxval=BBE_MOVN_2XF32_FROMNX16(tmp);

        a=BBE_MULMN_4XCF32(z12,z12,0,0); BBE_MULMASN_4XCF32(a,z12,z12,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b12=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b12); 
        a=BBE_MULMN_4XCF32(z13,z13,0,0); BBE_MULMASN_4XCF32(a,z13,z13,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b13=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b13); 
        a=BBE_MULMN_4XCF32(z14,z14,0,0); BBE_MULMASN_4XCF32(a,z14,z14,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b14); 
        a=BBE_MULMN_4XCF32(z15,z15,0,0); BBE_MULMASN_4XCF32(a,z15,z15,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b15=BBE_OGTN_2XF32(t,maxval);                                     xcol=BBE_MOVN_2XC16T(3,xcol,b15); xrow=BBE_MOVN_2XC16T(3,xrow,b12|b13|b14|b15); 
        brow1=BBE_EQN_2XC16(xrow,1);
        brow2=BBE_EQN_2XC16(xrow,2);
        brow3=BBE_EQN_2XC16(xrow,3);
        bcol1=BBE_EQN_2XC16(xcol,1);
        bcol2=BBE_EQN_2XC16(xcol,2);
        bcol3=BBE_EQN_2XC16(xcol,3);
        BBE_SBN_2_IP(brow1,bRowCol0,sizeof(vboolN_2));
        BBE_SBN_2_IP(brow2,bRowCol0,sizeof(vboolN_2));
        BBE_SBN_2_IP(brow3,bRowCol0,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol1,bRowCol0,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol2,bRowCol0,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol3,bRowCol0,sizeof(vboolN_2));
    }

    // swap rows/columns
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z);
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xf32 t,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol0,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z0 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-15*L*sizeof(complex_float));

        t=z0 ; z0 =BBE_MOVN_2XF32T(z4 ,z0 ,brow1); z4 =BBE_MOVN_2XF32T(t,z4 ,brow1);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z8 ,z0 ,brow2); z8 =BBE_MOVN_2XF32T(t,z8 ,brow2);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z12,z0 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z5 ,z1 ,brow1); z5 =BBE_MOVN_2XF32T(t,z5 ,brow1);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z9 ,z1 ,brow2); z9 =BBE_MOVN_2XF32T(t,z9 ,brow2);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z13,z1 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z6 ,z2 ,brow1); z6 =BBE_MOVN_2XF32T(t,z6 ,brow1);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z10,z2 ,brow2); z10=BBE_MOVN_2XF32T(t,z10,brow2);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z14,z2 ,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z7 ,z3 ,brow1); z7 =BBE_MOVN_2XF32T(t,z7 ,brow1);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z11,z3 ,brow2); z11=BBE_MOVN_2XF32T(t,z11,brow2);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z15,z3 ,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        t=z0 ; z0 =BBE_MOVN_2XF32T(z1, z0 ,bcol1); z1 =BBE_MOVN_2XF32T(t,z1 ,bcol1);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z2 ,z0 ,bcol2); z2 =BBE_MOVN_2XF32T(t,z2 ,bcol2);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z3 ,z0 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z5 ,z4 ,bcol1); z5 =BBE_MOVN_2XF32T(t,z5 ,bcol1);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z6 ,z4 ,bcol2); z6 =BBE_MOVN_2XF32T(t,z6 ,bcol2);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z7 ,z4 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z9 ,z8 ,bcol1); z9 =BBE_MOVN_2XF32T(t,z9 ,bcol1);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z10,z8 ,bcol2); z10=BBE_MOVN_2XF32T(t,z10,bcol2);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z11,z8 ,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z12; z12=BBE_MOVN_2XF32T(z13,z12,bcol1); z13=BBE_MOVN_2XF32T(t,z13,bcol1);
        t=z12; z12=BBE_MOVN_2XF32T(z14,z12,bcol2); z14=BBE_MOVN_2XF32T(t,z14,bcol2);
        t=z12; z12=BBE_MOVN_2XF32T(z15,z12,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        BBE_SVN_2XF32_XP(z0 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-15*L*sizeof(complex_float));
    }
    /* reciprocal of main diagonal */
    cmatinvgjnxnsf_recip(z,L);
    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 1*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_4XCF32_XP(z0 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z5 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z10,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z12,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-15*L*sizeof(complex_float));

        z1 =BBE_MULN_4XCF32( z1 ,z0 );  
        z2 =BBE_MULN_4XCF32( z2 ,z0 );  
        z3 =BBE_MULN_4XCF32( z3 ,z0 );

        BBE_MULSN_4XCF32(z5 ,z1 ,z4 ); 
        BBE_MULSN_4XCF32(z6 ,z2 ,z4 ); 
        BBE_MULSN_4XCF32(z7 ,z3 ,z4 ); 

        BBE_MULSN_4XCF32(z9 ,z1 ,z8 ); 
        BBE_MULSN_4XCF32(z10,z2 ,z8 ); 
        BBE_MULSN_4XCF32(z11,z3 ,z8 ); 

        BBE_MULSN_4XCF32(z13,z1 ,z12); 
        BBE_MULSN_4XCF32(z14,z2 ,z12); 
        BBE_MULSN_4XCF32(z15,z3 ,z12); 

        BBE_SVN_4XCF32_XP(z1 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z2 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z3 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z5 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z6 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z7 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z9 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z10,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z11,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z13,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z14,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-14*L*sizeof(complex_float));
    }
    pZrd   =(const xb_vecN_4xcf32*)(z+ 0*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 4*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 t,z0,z4,z8,z12;
        BBE_LVN_4XCF32_XP(z0 ,pZrd,  4*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,  4*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,  4*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z12,pZrd,2*BBE_SIMD_WIDTH-12*L*sizeof(complex_float));
        t=z4 ; z4 =BBE_MULMN_4XCF32 (z0,t, 3, 4); BBE_MULMASN_4XCF32 (z4 , z0,t, 2, 11);
        t=z8 ; z8 =BBE_MULMN_4XCF32 (z0,t, 3, 4); BBE_MULMASN_4XCF32 (z8 , z0,t, 2, 11);
        t=z12; z12=BBE_MULMN_4XCF32 (z0,t, 3, 4); BBE_MULMASN_4XCF32 (z12, z0,t, 2, 11);
        BBE_SVN_4XCF32_XP(z4 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z8 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z12,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(complex_float));
    }

    /*------------*/
    /* stage 2    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    // search the pivot in the submatrix 3x3
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+5*L);
    bRowCol1=browcol+6*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 t,maxval;
        xb_vecN_4xcf32 a,z5,z6,z7,z9,z10,z11,z13,z14,z15;
        vboolN_2 b5,b6,b7,b9,b10,b11,b13,b14,b15;
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xc16 xcol=1,xrow=1;
        BBE_LVN_4XCF32_XP(z5 ,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd,2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z10,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-10*L*sizeof(complex_float));

        a=BBE_MULMN_4XCF32(z5 ,z5 ,0,0); BBE_MULMASN_4XCF32(a,z5 ,z5 ,0,15); maxval=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        a=BBE_MULMN_4XCF32(z6 ,z6 ,0,0); BBE_MULMASN_4XCF32(a,z6 ,z6 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b6 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b6 ); 
        a=BBE_MULMN_4XCF32(z7 ,z7 ,0,0); BBE_MULMASN_4XCF32(a,z7 ,z7 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b7 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b7 ); xrow=BBE_MOVN_2XC16T(1,xrow,b5|b6|b7 ); 
        a=BBE_MULMN_4XCF32(z9 ,z9 ,0,0); BBE_MULMASN_4XCF32(a,z9 ,z9 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b9 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b9 ); 
        a=BBE_MULMN_4XCF32(z10,z10,0,0); BBE_MULMASN_4XCF32(a,z10,z10,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b10=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b10); 
        a=BBE_MULMN_4XCF32(z11,z11,0,0); BBE_MULMASN_4XCF32(a,z11,z11,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b11); xrow=BBE_MOVN_2XC16T(2,xrow,b9|b10|b11); 
        a=BBE_MULMN_4XCF32(z13,z13,0,0); BBE_MULMASN_4XCF32(a,z13,z13,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b13=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b13); 
        a=BBE_MULMN_4XCF32(z14,z14,0,0); BBE_MULMASN_4XCF32(a,z14,z14,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b14); 
        a=BBE_MULMN_4XCF32(z15,z15,0,0); BBE_MULMASN_4XCF32(a,z15,z15,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b15=BBE_OGTN_2XF32(t,maxval);                                     xcol=BBE_MOVN_2XC16T(3,xcol,b15); xrow=BBE_MOVN_2XC16T(3,xrow,b13|b14|b15); 
        brow2=BBE_EQN_2XC16(xrow,2);
        brow3=BBE_EQN_2XC16(xrow,3);
        bcol2=BBE_EQN_2XC16(xcol,2);
        bcol3=BBE_EQN_2XC16(xcol,3);
        BBE_SBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));
    }
    // swap rows/columns
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+L);
    pZwr   =(      xb_vecN_4xcf32*)(z+L);
    bRowCol1=browcol+6*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xf32 t,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-14*L*sizeof(complex_float));

        t=z4 ; z4 =BBE_MOVN_2XF32T(z8 ,z4 ,brow2); z8 =BBE_MOVN_2XF32T(t,z8 ,brow2);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z12,z4 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z9 ,z5 ,brow2); z9 =BBE_MOVN_2XF32T(t,z9 ,brow2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z13,z5 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z10,z6 ,brow2); z10=BBE_MOVN_2XF32T(t,z10,brow2);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z14,z6 ,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z11,z7 ,brow2); z11=BBE_MOVN_2XF32T(t,z11,brow2);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z15,z7 ,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        t=z1 ; z1 =BBE_MOVN_2XF32T(z2 ,z1 ,bcol2); z2 =BBE_MOVN_2XF32T(t,z2 ,bcol2);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z3 ,z1 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z6 ,z5 ,bcol2); z6 =BBE_MOVN_2XF32T(t,z6 ,bcol2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z7 ,z5 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z10,z9 ,bcol2); z10=BBE_MOVN_2XF32T(t,z10,bcol2);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z11,z9 ,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z13; z13=BBE_MOVN_2XF32T(z14,z13,bcol2); z14=BBE_MOVN_2XF32T(t,z14,bcol2);
        t=z13; z13=BBE_MOVN_2XF32T(z15,z13,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        BBE_SVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-14*L*sizeof(complex_float));
    }
    /* reciprocal of main diagonal */
    cmatinvgjnxnsf_recip(z+L*5,L);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_4xcf32*)(z+ 5*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 4*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_4XCF32_XP(z5 ,pZrd, -1*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd, -6*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd, -1*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z0 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3 ,pZrd,  6*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd, -1*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z10,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd, -1*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z12,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-10*L*sizeof(complex_float));

        z4 =BBE_MULN_4XCF32( z4 ,z5 );  
        z6 =BBE_MULN_4XCF32( z6 ,z5 );  
        z7 =BBE_MULN_4XCF32( z7 ,z5 );  
        BBE_MULSN_4XCF32(z0 ,z4 ,z1 ); 
        BBE_MULSN_4XCF32(z2 ,z6 ,z1 ); 
        BBE_MULSN_4XCF32(z3 ,z7 ,z1 ); 
        BBE_MULSN_4XCF32(z8 ,z4 ,z9 ); 
        BBE_MULSN_4XCF32(z10,z6 ,z9 ); 
        BBE_MULSN_4XCF32(z11,z7 ,z9 ); 
        BBE_MULSN_4XCF32(z12,z4 ,z13); 
        BBE_MULSN_4XCF32(z14,z6 ,z13); 
        BBE_MULSN_4XCF32(z15,z7 ,z13); 

        BBE_SVN_4XCF32_XP(z4 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z6 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z7 ,pZwr, -7*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z0 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z2 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z3 ,pZwr,  5*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z8 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z10,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z11,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z12,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z14,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-11*L*sizeof(complex_float));
    }
    pZrd   =(const xb_vecN_4xcf32*)(z+ 5*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 1*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 t,z1,z9,z13,z5;
        BBE_LVN_4XCF32_XP(z5 ,pZrd, -4*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd,  8*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,  4*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(complex_float));
        t=z1 ; z1 =BBE_MULMN_4XCF32 (z5,t, 3, 4); BBE_MULMASN_4XCF32 (z1 , z5,t, 2, 11);
        t=z9 ; z9 =BBE_MULMN_4XCF32 (z5,t, 3, 4); BBE_MULMASN_4XCF32 (z9 , z5,t, 2, 11);
        t=z13; z13=BBE_MULMN_4XCF32 (z5,t, 3, 4); BBE_MULMASN_4XCF32 (z13, z5,t, 2, 11);
        BBE_SVN_4XCF32_XP(z1 ,pZwr, 8* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z9 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z13,pZwr,2*BBE_SIMD_WIDTH-12*L*sizeof(complex_float));
    }

    /*------------*/
    /* stage 3    */
    /*------------*/
    // search the pivot in the submatrix 2x2
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+10*L);
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 t,maxval;
        xb_vecN_4xcf32 a,z10,z11,z14,z15;
        vboolN_2 b11,b14,b15;
        vboolN_2 brow3,bcol3;
        BBE_LVN_4XCF32_XP(z10,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,3*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-5*L*sizeof(complex_float));

        a=BBE_MULMN_4XCF32(z10,z10,0,0); BBE_MULMASN_4XCF32(a,z10,z10,0,15); maxval=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        a=BBE_MULMN_4XCF32(z11,z11,0,0); BBE_MULMASN_4XCF32(a,z11,z11,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        a=BBE_MULMN_4XCF32(z14,z14,0,0); BBE_MULMASN_4XCF32(a,z14,z14,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        a=BBE_MULMN_4XCF32(z15,z15,0,0); BBE_MULMASN_4XCF32(a,z15,z15,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));
        b15=BBE_OGTN_2XF32(t,maxval);                                     
        brow3=b14;
        bcol3=b11&~b14;
        brow3|=b15;
        bcol3|=b15;

       BBE_SBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));
       BBE_SBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));
    }

    // swap rows/columns
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+2*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+2*L);
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow3,bcol3;
        xb_vecN_2xf32 t,z2,z3,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),3*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-13*L*sizeof(complex_float));

        t=z8 ; z8 =BBE_MOVN_2XF32T(z12,z8 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z13,z9 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z10; z10=BBE_MOVN_2XF32T(z14,z10,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z11; z11=BBE_MOVN_2XF32T(z15,z11,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        t=z2 ; z2 =BBE_MOVN_2XF32T(z3 ,z2 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z7 ,z6 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z10; z10=BBE_MOVN_2XF32T(z11,z10,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z14; z14=BBE_MOVN_2XF32T(z15,z14,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),3*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-13*L*sizeof(complex_float));
    }
    /* reciprocal of main diagonal */
    cmatinvgjnxnsf_recip(z+L*10,L);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_4xcf32*)(z+10*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 8*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_4XCF32_XP(z10,pZrd, -2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd, -9*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd, -2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z0 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3 ,pZrd,  3*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd, -2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z5 ,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd,  7*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd, -2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z12,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,  2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-5*L*sizeof(complex_float));

        z8 =BBE_MULN_4XCF32(z8 ,z10);  
        z9 =BBE_MULN_4XCF32(z9 ,z10);  
        z11=BBE_MULN_4XCF32(z11,z10);  
        BBE_MULSN_4XCF32(z0 ,z8 ,z2 ); 
        BBE_MULSN_4XCF32(z1 ,z9 ,z2 ); 
        BBE_MULSN_4XCF32(z3 ,z11,z2 ); 
        BBE_MULSN_4XCF32(z4 ,z8 ,z6 ); 
        BBE_MULSN_4XCF32(z5 ,z9 ,z6 ); 
        BBE_MULSN_4XCF32(z7 ,z11,z6 ); 
        BBE_MULSN_4XCF32(z12,z8 ,z14); 
        BBE_MULSN_4XCF32(z13,z9 ,z14); 
        BBE_MULSN_4XCF32(z15,z11,z14); 

        BBE_SVN_4XCF32_XP(z8 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z9 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z11,pZwr,-11*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z0 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z1 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z3 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z4 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z5 ,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z7 ,pZwr,  5*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z12,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z13,pZwr,  2*L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-7*L*sizeof(complex_float));
    }
    pZrd   =(const xb_vecN_4xcf32*)(z+10*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 2*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 t,z2,z6,z10,z14;
        BBE_LVN_4XCF32_XP(z10,pZrd, -8* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd,  4* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd,  8* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,2*BBE_SIMD_WIDTH-4*L*sizeof(complex_float));
        t=z2 ; z2 =BBE_MULMN_4XCF32 (z10,t, 3, 4); BBE_MULMASN_4XCF32 (z2 , z10,t, 2, 11);
        t=z6 ; z6 =BBE_MULMN_4XCF32 (z10,t, 3, 4); BBE_MULMASN_4XCF32 (z6 , z10,t, 2, 11);
        t=z14; z14=BBE_MULMN_4XCF32 (z10,t, 3, 4); BBE_MULMASN_4XCF32 (z14, z10,t, 2, 11);
        BBE_SVN_4XCF32_XP(z2 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z6 ,pZwr, 8* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z14,pZwr,2*BBE_SIMD_WIDTH-12*L*sizeof(complex_float));
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    /* reciprocal of main diagonal */
    cmatinvgjnxnsf_recip(z+L*15,L);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_4xcf32*)(z+15*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+12*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_4XCF32_XP(z15,pZrd,-3* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z12,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z13,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z14,pZrd,-11*L*sizeof(complex_float));
        z12=BBE_MULN_4XCF32(z12,z15);       
        z13=BBE_MULN_4XCF32(z13,z15);       
        z14=BBE_MULN_4XCF32(z14,z15);       
        BBE_SVN_4XCF32_XP(z12,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z13,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z14,pZwr,-14*L*sizeof(complex_float));

        BBE_LVN_4XCF32_XP(z3 ,pZrd,-3* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z0 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2 ,pZrd, 5* L*sizeof(complex_float));
        BBE_MULSN_4XCF32(z0 ,z12,z3 );      
        BBE_MULSN_4XCF32(z1 ,z13,z3 );      
        BBE_MULSN_4XCF32(z2 ,z14,z3 );      
        BBE_SVN_4XCF32_XP(z0 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z1 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z2 ,pZwr,  2*L*sizeof(complex_float));

        BBE_LVN_4XCF32_XP(z7 ,pZrd,-3* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z4 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z5 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z6 ,pZrd, 5* L*sizeof(complex_float));
        BBE_MULSN_4XCF32(z4 ,z12,z7 );     
        BBE_MULSN_4XCF32(z5 ,z13,z7 );     
        BBE_MULSN_4XCF32(z6 ,z14,z7 );     
        BBE_SVN_4XCF32_XP(z4 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z5 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z6 ,pZwr,  2*L*sizeof(complex_float));

        BBE_LVN_4XCF32_XP(z11,pZrd,-3* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z8 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z9 ,pZrd,    L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z10,pZrd,2*BBE_SIMD_WIDTH+5*L*sizeof(complex_float));
        BBE_MULSN_4XCF32(z8 ,z12,z11);     
        BBE_MULSN_4XCF32(z9 ,z13,z11);     
        BBE_MULSN_4XCF32(z10,z14,z11);     
        BBE_SVN_4XCF32_XP(z8 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z9 ,pZwr,    L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z10,pZwr,2*BBE_SIMD_WIDTH+2*L*sizeof(complex_float));
    }
    pZrd   =(const xb_vecN_4xcf32*)(z+15*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+ 3*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 t,z3,z7,z11,z15;
        BBE_LVN_4XCF32_XP(z15,pZrd,-12*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3 ,pZrd,  4* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z7 ,pZrd,  4* L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z11,pZrd,2*BBE_SIMD_WIDTH+4*L*sizeof(complex_float));
        t=z3 ; z3 =BBE_MULMN_4XCF32 (z15,t, 3, 4); BBE_MULMASN_4XCF32 (z3 , z15,t, 2, 11);
        t=z7 ; z7 =BBE_MULMN_4XCF32 (z15,t, 3, 4); BBE_MULMASN_4XCF32 (z7 , z15,t, 2, 11);
        t=z11; z11=BBE_MULMN_4XCF32 (z15,t, 3, 4); BBE_MULMASN_4XCF32 (z11, z15,t, 2, 11);
        BBE_SVN_4XCF32_XP(z3 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z7 ,pZwr, 4* L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z11,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(complex_float));
    }
    /* final reverse permulation of columns  */
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+2*L);
    pZwr   =(      xb_vecN_4xcf32*)(z+2*L);
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow3,bcol3;
        xb_vecN_2xf32 t,z2,z3,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),3*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),  L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-13*L*sizeof(complex_float));

        t=z2 ; z2 =BBE_MOVN_2XF32T(z3 ,z2 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z7 ,z6 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z10; z10=BBE_MOVN_2XF32T(z11,z10,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z14; z14=BBE_MOVN_2XF32T(z15,z14,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        t=z8 ; z8 =BBE_MOVN_2XF32T(z12,z8 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z13,z9 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z10; z10=BBE_MOVN_2XF32T(z14,z10,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z11; z11=BBE_MOVN_2XF32T(z15,z11,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),3*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),  L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-13*L*sizeof(complex_float));
    }
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z+L);
    pZwr   =(      xb_vecN_4xcf32*)(z+L);
    bRowCol1=browcol+6*(L>>(LOG2_BBE_SIMD_WIDTH-2));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xf32 t,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-14*L*sizeof(complex_float));

        t=z4 ; z4 =BBE_MOVN_2XF32T(z8 ,z4 ,brow2); z8 =BBE_MOVN_2XF32T(t,z8 ,brow2);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z12,z4 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z9 ,z5 ,brow2); z9 =BBE_MOVN_2XF32T(t,z9 ,brow2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z13,z5 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z10,z6 ,brow2); z10=BBE_MOVN_2XF32T(t,z10,brow2);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z14,z6 ,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z11,z7 ,brow2); z11=BBE_MOVN_2XF32T(t,z11,brow2);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z15,z7 ,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        t=z1 ; z1 =BBE_MOVN_2XF32T(z2 ,z1 ,bcol2); z2 =BBE_MOVN_2XF32T(t,z2 ,bcol2);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z3 ,z1 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z6 ,z5 ,bcol2); z6 =BBE_MOVN_2XF32T(t,z6 ,bcol2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z7 ,z5 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z10,z9 ,bcol2); z10=BBE_MOVN_2XF32T(t,z10,bcol2);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z11,z9 ,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z13; z13=BBE_MOVN_2XF32T(z14,z13,bcol2); z14=BBE_MOVN_2XF32T(t,z14,bcol2);
        t=z13; z13=BBE_MOVN_2XF32T(z15,z13,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        BBE_SVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-14*L*sizeof(complex_float));
    }
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z);
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xf32 t,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LBN_2_IP(bcol1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol0,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z0 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZrd),-11*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZrd),-11*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZrd),-11*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZrd),  4*L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-15*L*sizeof(complex_float));

        t=z0 ; z0 =BBE_MOVN_2XF32T(z1, z0 ,bcol1); z1 =BBE_MOVN_2XF32T(t,z1 ,bcol1);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z2 ,z0 ,bcol2); z2 =BBE_MOVN_2XF32T(t,z2 ,bcol2);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z3 ,z0 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z5 ,z4 ,bcol1); z5 =BBE_MOVN_2XF32T(t,z5 ,bcol1);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z6 ,z4 ,bcol2); z6 =BBE_MOVN_2XF32T(t,z6 ,bcol2);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z7 ,z4 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z9 ,z8 ,bcol1); z9 =BBE_MOVN_2XF32T(t,z9 ,bcol1);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z10,z8 ,bcol2); z10=BBE_MOVN_2XF32T(t,z10,bcol2);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z11,z8 ,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z12; z12=BBE_MOVN_2XF32T(z13,z12,bcol1); z13=BBE_MOVN_2XF32T(t,z13,bcol1);
        t=z12; z12=BBE_MOVN_2XF32T(z14,z12,bcol2); z14=BBE_MOVN_2XF32T(t,z14,bcol2);
        t=z12; z12=BBE_MOVN_2XF32T(z15,z12,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z4 ,z0 ,brow1); z4 =BBE_MOVN_2XF32T(t,z4 ,brow1);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z8 ,z0 ,brow2); z8 =BBE_MOVN_2XF32T(t,z8 ,brow2);
        t=z0 ; z0 =BBE_MOVN_2XF32T(z12,z0 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z5 ,z1 ,brow1); z5 =BBE_MOVN_2XF32T(t,z5 ,brow1);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z9 ,z1 ,brow2); z9 =BBE_MOVN_2XF32T(t,z9 ,brow2);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z13,z1 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z6 ,z2 ,brow1); z6 =BBE_MOVN_2XF32T(t,z6 ,brow1);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z10,z2 ,brow2); z10=BBE_MOVN_2XF32T(t,z10,brow2);
        t=z2 ; z2 =BBE_MOVN_2XF32T(z14,z2 ,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z7 ,z3 ,brow1); z7 =BBE_MOVN_2XF32T(t,z7 ,brow1);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z11,z3 ,brow2); z11=BBE_MOVN_2XF32T(t,z11,brow2);
        t=z3 ; z3 =BBE_MOVN_2XF32T(z15,z3 ,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        BBE_SVN_2XF32_XP(z0 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z4 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z8 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z12,castxcc(xb_vecN_2xf32,pZwr),-11*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z1 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z5 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z9 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z13,castxcc(xb_vecN_2xf32,pZwr),-11*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z6 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z10,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z14,castxcc(xb_vecN_2xf32,pZwr),-11*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z7 ,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z11,castxcc(xb_vecN_2xf32,pZwr),  4*L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z15,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-15*L*sizeof(complex_float));
    }

}
#endif

size_t cmatinvgj4x4sf_getScratchSize (int N, int L  )
{
    NASSERT(N==4);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int32_t)*4*L+                              /* L row indices, L column indices, 2*L maxval */
           sizeof(vboolN_2)*12*(L>>(LOG2_BBE_SIMD_WIDTH-2)); /* boolean masks */
}

#else
DISCARD_FUN(void, cmatinvgj4x4sf ,( void* pScr, complex_float * restrict z, int L  ))

size_t cmatinvgj4x4sf_getScratchSize (int N, int L  )
{
    NASSERT(N==4);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

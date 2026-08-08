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
    Real Matrix Gauss-Jordan inversion, floating point real data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "matinvgjnxnsf_common.h"

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
void matinvgj4x4sf ( void* pScr, float32_t * restrict z, int L  )
#if 0
{
    int32_t *srow,*scol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);

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
        maxval=fabsf(z[l+L* 0]);          {           srow[0*L+l]=0; scol[0*L+l]=0; }
        t=fabsf(z[l+L* 1]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=1; }
        t=fabsf(z[l+L* 2]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=2; }
        t=fabsf(z[l+L* 3]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=3; }

        t=fabsf(z[l+L* 4]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=0; }
        t=fabsf(z[l+L* 5]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=1; }
        t=fabsf(z[l+L* 6]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=2; }
        t=fabsf(z[l+L* 7]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=3; }

        t=fabsf(z[l+L* 8]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=0; }
        t=fabsf(z[l+L* 9]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=1; }
        t=fabsf(z[l+L*10]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=2; }
        t=fabsf(z[l+L*11]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=3; }

        t=fabsf(z[l+L*12]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=0; }
        t=fabsf(z[l+L*13]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=1; }
        t=fabsf(z[l+L*14]);  if(t>maxval) { maxval=t; srow[0*L+l]=3; scol[0*L+l]=2; }
        t=fabsf(z[l+L*15]);  if(t>maxval) {           srow[0*L+l]=3; scol[0*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        float32_t t;
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
    for (l=0; l<L; l++) z[l+L* 0]=1.0f/z[l+L* 0];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L* 1] *= z[l+L* 0];
        z[l+L* 2] *= z[l+L* 0];
        z[l+L* 3] *= z[l+L* 0];
        z[l+L* 5] -= z[l+L* 1] * z[l+L* 4]; 
        z[l+L* 6] -= z[l+L* 2] * z[l+L* 4]; 
        z[l+L* 7] -= z[l+L* 3] * z[l+L* 4]; 
        z[l+L* 4] = -z[l+L* 0] * z[l+L* 4]; 
        z[l+L* 9] -= z[l+L* 1] * z[l+L* 8]; 
        z[l+L*10] -= z[l+L* 2] * z[l+L* 8]; 
        z[l+L*11] -= z[l+L* 3] * z[l+L* 8]; 
        z[l+L* 8] = -z[l+L* 0] * z[l+L* 8]; 
        z[l+L*13] -= z[l+L* 1] * z[l+L*12]; 
        z[l+L*14] -= z[l+L* 2] * z[l+L*12]; 
        z[l+L*15] -= z[l+L* 3] * z[l+L*12]; 
        z[l+L*12] = -z[l+L* 0] * z[l+L*12]; 
    }
    /*------------*/
    /* stage 2    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {   
        float32_t t,maxval;
        maxval=fabsf(z[l+L*(1*4+1)]);          {           srow[1*L+l]=1; scol[1*L+l]=1; }
        t=fabsf(z[l+L*(1*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=1; scol[1*L+l]=2; }
        t=fabsf(z[l+L*(1*4+3)]);  if(t>maxval) { maxval=t; srow[1*L+l]=1; scol[1*L+l]=3; }

        t=fabsf(z[l+L*(2*4+1)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=1; }
        t=fabsf(z[l+L*(2*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=2; }
        t=fabsf(z[l+L*(2*4+3)]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=3; }

        t=fabsf(z[l+L*(3*4+1)]);  if(t>maxval) { maxval=t; srow[1*L+l]=3; scol[1*L+l]=1; }
        t=fabsf(z[l+L*(3*4+2)]);  if(t>maxval) { maxval=t; srow[1*L+l]=3; scol[1*L+l]=2; }
        t=fabsf(z[l+L*(3*4+3)]);  if(t>maxval) {           srow[1*L+l]=3; scol[1*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        float32_t t;
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
    for (l=0; l<L; l++) z[l+L*5]=1.0f/z[l+L*5];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L*4] *= z[l+L*5];
        z[l+L*6] *= z[l+L*5];
        z[l+L*7] *= z[l+L*5];
        z[l+L* 0] -= z[l+L* 4] * z[l+L* 1]; 
        z[l+L* 2] -= z[l+L* 6] * z[l+L* 1]; 
        z[l+L* 3] -= z[l+L* 7] * z[l+L* 1]; 
        z[l+L* 1] = -z[l+L* 5] * z[l+L* 1]; 
        z[l+L* 8] -= z[l+L* 4] * z[l+L* 9]; 
        z[l+L*10] -= z[l+L* 6] * z[l+L* 9]; 
        z[l+L*11] -= z[l+L* 7] * z[l+L* 9]; 
        z[l+L* 9] = -z[l+L* 5] * z[l+L* 9]; 
        z[l+L*12] -= z[l+L* 4] * z[l+L*13]; 
        z[l+L*14] -= z[l+L* 6] * z[l+L*13]; 
        z[l+L*15] -= z[l+L* 7] * z[l+L*13]; 
        z[l+L*13] = -z[l+L* 5] * z[l+L*13]; 
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {   
        float32_t t,maxval;
        maxval=fabsf(z[l+L*(2*4+2)]);          {           srow[2*L+l]=2; scol[2*L+l]=2; }
        t=fabsf(z[l+L*(2*4+3)]);  if(t>maxval) { maxval=t; srow[2*L+l]=2; scol[2*L+l]=3; }
        t=fabsf(z[l+L*(3*4+2)]);  if(t>maxval) { maxval=t; srow[2*L+l]=3; scol[2*L+l]=2; }
        t=fabsf(z[l+L*(3*4+3)]);  if(t>maxval) {           srow[2*L+l]=3; scol[2*L+l]=3; }
    }
    for (l=0; l<L; l++)
    {
        float32_t t;
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
    for (l=0; l<L; l++) z[l+L*10]=1.0f/z[l+L*10];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L* 8] *= z[l+L*10];
        z[l+L* 9] *= z[l+L*10];
        z[l+L*11] *= z[l+L*10];
        z[l+L* 0] -= z[l+L* 8] * z[l+L* 2]; 
        z[l+L* 1] -= z[l+L* 9] * z[l+L* 2]; 
        z[l+L* 3] -= z[l+L*11] * z[l+L* 2]; 
        z[l+L* 2] = -z[l+L*10] * z[l+L* 2]; 
        z[l+L* 4] -= z[l+L* 8] * z[l+L* 6]; 
        z[l+L* 5] -= z[l+L* 9] * z[l+L* 6]; 
        z[l+L* 7] -= z[l+L*11] * z[l+L* 6]; 
        z[l+L* 6] = -z[l+L*10] * z[l+L* 6]; 
        z[l+L*12] -= z[l+L* 8] * z[l+L*14]; 
        z[l+L*13] -= z[l+L* 9] * z[l+L*14]; 
        z[l+L*15] -= z[l+L*11] * z[l+L*14]; 
        z[l+L*14] = -z[l+L*10] * z[l+L*14]; 
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*15]=1.0f/z[l+L*15];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        z[l+L*12] *= z[l+L*15];
        z[l+L*13] *= z[l+L*15];
        z[l+L*14] *= z[l+L*15];
        z[l+L* 0] -= z[l+L*12] * z[l+L* 3]; 
        z[l+L* 1] -= z[l+L*13] * z[l+L* 3]; 
        z[l+L* 2] -= z[l+L*14] * z[l+L* 3]; 
        z[l+L* 3] = -z[l+L*15] * z[l+L* 3]; 
        z[l+L* 4] -= z[l+L*12] * z[l+L* 7]; 
        z[l+L* 5] -= z[l+L*13] * z[l+L* 7]; 
        z[l+L* 6] -= z[l+L*14] * z[l+L* 7]; 
        z[l+L* 7] = -z[l+L*15] * z[l+L* 7]; 
        z[l+L* 8] -= z[l+L*12] * z[l+L*11]; 
        z[l+L* 9] -= z[l+L*13] * z[l+L*11]; 
        z[l+L*10] -= z[l+L*14] * z[l+L*11]; 
        z[l+L*11] = -z[l+L*15] * z[l+L*11]; 
    }
    /* final reverse permulation of columns  */
    for (l=0; l<L; l++)
    {
        float32_t t;
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
    const xb_vecN_2xf32 * restrict pZrd;
          xb_vecN_2xf32 * restrict pZwr;
    vboolN_2 * browcol;  /* pairs of row/col boolean flags for permutations, [12*(L>>(LOG2_BBE_SIMD_WIDTH-1))]*/
    vboolN_2 * restrict bRowCol0;
    vboolN_2 * restrict bRowCol1;
    vboolN_2 * restrict bRowCol2;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);

    if(L<=0 ) return;

    browcol=(vboolN_2*)pScr;

    /*------------*/
    /* stage 1    */
    /*------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    pZrd   =(const xb_vecN_2xf32*)(z);
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 t,maxval,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        vboolN_2 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15;
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xc16 xcol=0,xrow=0;
        BBE_LVN_2XF32_XP(z0 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));

        maxval=BBE_ABSN_2XF32(z0);
        t=BBE_ABSN_2XF32(z1 ); b1 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b1 ); 
        t=BBE_ABSN_2XF32(z2 ); b2 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b2 ); 
        t=BBE_ABSN_2XF32(z3 ); b3 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b3 ); 
        t=BBE_ABSN_2XF32(z4 ); b4 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b4 ); 
        t=BBE_ABSN_2XF32(z5 ); b5 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b5 ); 
        t=BBE_ABSN_2XF32(z6 ); b6 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b6 ); 
        t=BBE_ABSN_2XF32(z7 ); b7 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b7 ); xrow=BBE_MOVN_2XC16T(1,xrow,b4|b5|b6|b7 ); 
        t=BBE_ABSN_2XF32(z8 ); b8 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b8 ); 
        t=BBE_ABSN_2XF32(z9 ); b9 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b9 ); 
        t=BBE_ABSN_2XF32(z10); b10=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b10); 
        t=BBE_ABSN_2XF32(z11); b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b11); xrow=BBE_MOVN_2XC16T(2,xrow,b8|b9|b10|b11); 
        t=BBE_ABSN_2XF32(z12); b12=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b12); 
        t=BBE_ABSN_2XF32(z13); b13=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b13); 
        t=BBE_ABSN_2XF32(z14); b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b14); 
        t=BBE_ABSN_2XF32(z15); b15=BBE_OGTN_2XF32(t,maxval);                                     xcol=BBE_MOVN_2XC16T(3,xcol,b15); xrow=BBE_MOVN_2XC16T(3,xrow,b12|b13|b14|b15); 
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
    pZrd   =(const xb_vecN_2xf32*)(z);
    pZwr   =(      xb_vecN_2xf32*)(z);
    bRowCol0=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xf32 t,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol0,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z0 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));

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

        BBE_SVN_2XF32_XP(z0 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));
    }

    matinvgjnxnsf_recip(z,L);  /* reciprocal of main diagonal */

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_2XF32_XP(z0 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));

        z1=BBE_MULN_2XF32(z1,z0);
        z2=BBE_MULN_2XF32(z2,z0);
        z3=BBE_MULN_2XF32(z3,z0);
        BBE_MULSN_2XF32(z5,z1,z4);
        BBE_MULSN_2XF32(z6,z2,z4);
        BBE_MULSN_2XF32(z7,z3,z4);
        z4=BBE_MULMN_2XF32(z0,z4,3,12);
        BBE_MULSN_2XF32(z9 ,z1,z8);
        BBE_MULSN_2XF32(z10,z2,z8);
        BBE_MULSN_2XF32(z11,z3,z8);
        z8=BBE_MULMN_2XF32(z0,z8,3,12);
        BBE_MULSN_2XF32(z13,z1,z12);
        BBE_MULSN_2XF32(z14,z2,z12);
        BBE_MULSN_2XF32(z15,z3,z12);
        z12=BBE_MULMN_2XF32(z0,z12,3,12);

        BBE_SVN_2XF32_XP(z0 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));
    }
    /*------------*/
    /* stage 2    */
    /*------------*/
    // search the pivot in the submatrix 3x3
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+5*L);
    bRowCol1=browcol+6*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 t,maxval,z5,z6,z7,z9,z10,z11,z13,z14,z15;
        vboolN_2 b5,b6,b7,b9,b10,b11,b13,b14,b15;
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xc16 xcol=1,xrow=1;
        BBE_LVN_2XF32_XP(z5 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-10*L*sizeof(float32_t));

        maxval=BBE_ABSN_2XF32(z5);
        t=BBE_ABSN_2XF32(z6 ); b6 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b6 ); 
        t=BBE_ABSN_2XF32(z7 ); b7 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b7 ); xrow=BBE_MOVN_2XC16T(1,xrow,b5|b6|b7 ); 
        t=BBE_ABSN_2XF32(z9 ); b9 =BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b9 ); 
        t=BBE_ABSN_2XF32(z10); b10=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b10); 
        t=BBE_ABSN_2XF32(z11); b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(3,xcol,b11); xrow=BBE_MOVN_2XC16T(2,xrow,b9|b10|b11); 
        t=BBE_ABSN_2XF32(z13); b13=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b13); 
        t=BBE_ABSN_2XF32(z14); b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b14); 
        t=BBE_ABSN_2XF32(z15); b15=BBE_OGTN_2XF32(t,maxval);                                     xcol=BBE_MOVN_2XC16T(3,xcol,b15); xrow=BBE_MOVN_2XC16T(3,xrow,b13|b14|b15); 
        brow2=BBE_EQN_2XC16(xrow,2);
        brow3=BBE_EQN_2XC16(xrow,3);
        bcol2=BBE_EQN_2XC16(xcol,2);
        bcol3=BBE_EQN_2XC16(xcol,3);
        BBE_SBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));
    }
        /* swap rows / columns*/
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+L);
    pZwr   =(      xb_vecN_2xf32*)(z+L);
    bRowCol1=browcol+6*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xf32 t,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z1 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-14*L*sizeof(float32_t));

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

        BBE_SVN_2XF32_XP(z1 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-14*L*sizeof(float32_t));
    }
    
    matinvgjnxnsf_recip(z+5*L,L);  /* reciprocal of main diagonal */
    /* read/write order 5 4 6 7 1 0 2 3 9 8 10 11 13 12 14 15 for better scheduling */
    pZrd   =(const xb_vecN_2xf32*)(z+5*L);
    pZwr   =(      xb_vecN_2xf32*)(z+5*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_2XF32_XP(z5 ,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,-6*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z0 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd, 6*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-10*L*sizeof(float32_t));

        z4=BBE_MULN_2XF32(z4,z5);            
        z6=BBE_MULN_2XF32(z6,z5);            
        z7=BBE_MULN_2XF32(z7,z5);            
        BBE_MULSN_2XF32(z0,z4,z1);           
        BBE_MULSN_2XF32(z2,z6,z1);           
        BBE_MULSN_2XF32(z3,z7,z1);           
        z1=BBE_MULMN_2XF32(z5,z1,3,12);      
        BBE_MULSN_2XF32(z8 ,z4,z9);          
        BBE_MULSN_2XF32(z10,z6,z9);          
        BBE_MULSN_2XF32(z11,z7,z9);          
        z9=BBE_MULMN_2XF32(z5,z9,3,12);      
        BBE_MULSN_2XF32(z12,z4,z13);         
        BBE_MULSN_2XF32(z14,z6,z13);         
        BBE_MULSN_2XF32(z15,z7,z13);         
        z13=BBE_MULMN_2XF32(z5,z13,3,12);    

        BBE_SVN_2XF32_XP(z5 ,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,-6*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z0 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr, 6*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-10*L*sizeof(float32_t));
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    // search the pivot in the submatrix 2x2
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+10*L);
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 t,maxval,z10,z11,z14,z15;
        vboolN_2 b11,b14,b15;
        vboolN_2 brow3,bcol3;
        BBE_LVN_2XF32_XP(z10,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-5*L*sizeof(float32_t));

        maxval=BBE_ABSN_2XF32(z10);
        t=BBE_ABSN_2XF32(z11); b11=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        t=BBE_ABSN_2XF32(z14); b14=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        t=BBE_ABSN_2XF32(z15); b15=BBE_OGTN_2XF32(t,maxval);                                     
        brow3=b14;
        bcol3=b11&~b14;
        brow3|=b15;
        bcol3|=b15;

        BBE_SBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));
    }
    // swap rows/columns
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+2*L);
    pZwr   =(      xb_vecN_2xf32*)(z+2*L);
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow3,bcol3;
        xb_vecN_2xf32 t,z2,z3,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(z2 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,  L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-13*L*sizeof(float32_t));

        t=z8 ; z8 =BBE_MOVN_2XF32T(z12,z8 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z13,z9 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z10; z10=BBE_MOVN_2XF32T(z14,z10,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z11; z11=BBE_MOVN_2XF32T(z15,z11,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        t=z2 ; z2 =BBE_MOVN_2XF32T(z3 ,z2 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z7 ,z6 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z10; z10=BBE_MOVN_2XF32T(z11,z10,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z14; z14=BBE_MOVN_2XF32T(z15,z14,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);

        BBE_SVN_2XF32_XP(z2 ,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,  L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-13*L*sizeof(float32_t));
    }

    matinvgjnxnsf_recip(z+10*L,L);  /* reciprocal of main diagonal */
    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    /* read/write order 10 8 9 11 2 0 1 3 6 4 5 7 14 12 13 15 for better scheduling */
    pZrd   =(const xb_vecN_2xf32*)(z+10*L);
    pZwr   =(      xb_vecN_2xf32*)(z+10*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_2XF32_XP(z10,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,-9*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z0 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd, 3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd, 7*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-5*L*sizeof(float32_t));
        
        z8 =BBE_MULN_2XF32(z8 ,z10);        
        z9 =BBE_MULN_2XF32(z9 ,z10);        
        z11=BBE_MULN_2XF32(z11,z10);        
        BBE_MULSN_2XF32(z0 ,z8 ,z2 );       
        BBE_MULSN_2XF32(z1 ,z9 ,z2 );       
        BBE_MULSN_2XF32(z3 ,z11,z2 );       
        z2 =BBE_MULMN_2XF32(z10,z2 ,3,12);  
        BBE_MULSN_2XF32(z4 ,z8 ,z6 );       
        BBE_MULSN_2XF32(z5 ,z9 ,z6 );       
        BBE_MULSN_2XF32(z7 ,z11,z6 );       
        z6 =BBE_MULMN_2XF32(z10,z6 ,3,12);  
        BBE_MULSN_2XF32(z12,z8 ,z14);       
        BBE_MULSN_2XF32(z13,z9 ,z14);       
        BBE_MULSN_2XF32(z15,z11,z14);       
        z14=BBE_MULMN_2XF32(z10,z14,3,12);  

        BBE_SVN_2XF32_XP(z10,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,-9*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z0 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr, 3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr, 7*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-5*L*sizeof(float32_t));
    }
    /*------------*/
    /* stage 3    */
    /*------------*/
    matinvgjnxnsf_recip(z+15*L,L);  /* reciprocal of main diagonal */
    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    /* read/write order 15 12 13 14 3 0 1 2 7 4 5 6 11 8 9 10 for better scheduling */
    pZrd   =(const xb_vecN_2xf32*)(z+15*L);
    pZwr   =(      xb_vecN_2xf32*)(z+15*L);
    for(l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_2XF32_XP(z15,pZrd,-3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,-3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z0 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd, 5*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,-3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd, 5*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,-3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,2*BBE_SIMD_WIDTH+5*L*sizeof(float32_t));
        
        z12=BBE_MULN_2XF32(z12,z15);       
        z13=BBE_MULN_2XF32(z13,z15);       
        z14=BBE_MULN_2XF32(z14,z15);       
        BBE_MULSN_2XF32(z0 ,z12,z3 );      
        BBE_MULSN_2XF32(z1 ,z13,z3 );      
        BBE_MULSN_2XF32(z2 ,z14,z3 );      
        z3 =BBE_MULMN_2XF32(z15,z3 ,3,12); 
        BBE_MULSN_2XF32(z4 ,z12,z7 );      
        BBE_MULSN_2XF32(z5 ,z13,z7 );      
        BBE_MULSN_2XF32(z6 ,z14,z7 );      
        z7 =BBE_MULMN_2XF32(z15,z7 ,3,12); 
        BBE_MULSN_2XF32(z8 ,z12,z11);      
        BBE_MULSN_2XF32(z9 ,z13,z11);      
        BBE_MULSN_2XF32(z10,z14,z11);      
        z11=BBE_MULMN_2XF32(z15,z11,3,12); 

        BBE_SVN_2XF32_XP(z15,pZwr,-3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,-3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z0 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr, 5*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,-3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr, 5*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,-3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,2*BBE_SIMD_WIDTH+5*L*sizeof(float32_t));
    }
    /* final reverse permulation of columns/rows  */
    __Pragma("no_reorder")
    bRowCol0=browcol;
    bRowCol1=browcol+ 6*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    bRowCol2=browcol+10*(L>>(LOG2_BBE_SIMD_WIDTH-1));
    pZrd   =(const xb_vecN_2xf32*)(z+8*L);
    pZwr   =(      xb_vecN_2xf32*)(z+8*L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow2,brow3,bcol2,bcol3;
        xb_vecN_2xf32 t,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;

        BBE_LVN_2XF32_XP(z8 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd, -3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd, -3*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,2*BBE_SIMD_WIDTH+3*L*sizeof(float32_t));

        BBE_LBN_2_IP(bcol3,bRowCol2,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol2,sizeof(vboolN_2));

        t=z2 ; z2 =BBE_MOVN_2XF32T(z3 ,z2 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z7 ,z6 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z10; z10=BBE_MOVN_2XF32T(z11,z10,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z14; z14=BBE_MOVN_2XF32T(z15,z14,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);
        t=z8 ; z8 =BBE_MOVN_2XF32T(z12,z8 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z13,z9 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z10; z10=BBE_MOVN_2XF32T(z14,z10,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z11; z11=BBE_MOVN_2XF32T(z15,z11,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        BBE_LBN_2_IP(bcol2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol1,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol1,sizeof(vboolN_2));
        t=z1 ; z1 =BBE_MOVN_2XF32T(z2 ,z1 ,bcol2); z2 =BBE_MOVN_2XF32T(t,z2 ,bcol2);
        t=z1 ; z1 =BBE_MOVN_2XF32T(z3 ,z1 ,bcol3); z3 =BBE_MOVN_2XF32T(t,z3 ,bcol3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z6 ,z5 ,bcol2); z6 =BBE_MOVN_2XF32T(t,z6 ,bcol2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z7 ,z5 ,bcol3); z7 =BBE_MOVN_2XF32T(t,z7 ,bcol3);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z10,z9 ,bcol2); z10=BBE_MOVN_2XF32T(t,z10,bcol2);
        t=z9 ; z9 =BBE_MOVN_2XF32T(z11,z9 ,bcol3); z11=BBE_MOVN_2XF32T(t,z11,bcol3);
        t=z13; z13=BBE_MOVN_2XF32T(z14,z13,bcol2); z14=BBE_MOVN_2XF32T(t,z14,bcol2);
        t=z13; z13=BBE_MOVN_2XF32T(z15,z13,bcol3); z15=BBE_MOVN_2XF32T(t,z15,bcol3);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z8 ,z4 ,brow2); z8 =BBE_MOVN_2XF32T(t,z8 ,brow2);
        t=z4 ; z4 =BBE_MOVN_2XF32T(z12,z4 ,brow3); z12=BBE_MOVN_2XF32T(t,z12,brow3);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z9 ,z5 ,brow2); z9 =BBE_MOVN_2XF32T(t,z9 ,brow2);
        t=z5 ; z5 =BBE_MOVN_2XF32T(z13,z5 ,brow3); z13=BBE_MOVN_2XF32T(t,z13,brow3);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z10,z6 ,brow2); z10=BBE_MOVN_2XF32T(t,z10,brow2);
        t=z6 ; z6 =BBE_MOVN_2XF32T(z14,z6 ,brow3); z14=BBE_MOVN_2XF32T(t,z14,brow3);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z11,z7 ,brow2); z11=BBE_MOVN_2XF32T(t,z11,brow2);
        t=z7 ; z7 =BBE_MOVN_2XF32T(z15,z7 ,brow3); z15=BBE_MOVN_2XF32T(t,z15,brow3);

        BBE_SVN_2XF32_XP(z8 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr, -3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr, -3*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,2*BBE_SIMD_WIDTH+3*L*sizeof(float32_t));
    }
    __Pragma("no_reorder")

    pZrd   =(const xb_vecN_2xf32*)(z);
    pZwr   =(      xb_vecN_2xf32*)(z);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow1,brow2,brow3,bcol1,bcol2,bcol3;
        xb_vecN_2xf32 t,z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15;
        BBE_LVN_2XF32_XP(z0 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z12,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z9 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z13,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z10,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z14,pZrd,-11*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7 ,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z11,pZrd,  4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z15,pZrd,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));

        BBE_LBN_2_IP(bcol1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol3,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow1,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol0,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow3,bRowCol0,sizeof(vboolN_2));
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

        BBE_SVN_2XF32_XP(z0 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z12,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z9 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z13,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z10,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z14,pZwr,-11*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7 ,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z11,pZwr,  4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z15,pZwr,2*BBE_SIMD_WIDTH-15*L*sizeof(float32_t));
    }

}
#endif

size_t matinvgj4x4sf_getScratchSize (int N, int L  )
{
    NASSERT(N==4);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 );
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(vboolN_2)*12*(L>>(LOG2_BBE_SIMD_WIDTH-1));
}
#else
DISCARD_FUN(void, matinvgj4x4sf,( void* pScr, float32_t * restrict z, int L  ))

size_t matinvgj4x4sf_getScratchSize (int N, int L  )
{
    NASSERT(N==4);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

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

#if HAVE_VFPU

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
void matinvgj3x3sf ( void* pScr, float32_t * restrict z, int L  )
#if 0
{
    int32_t *srow,*scol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 );

    if(L<=0 ) return;

    srow=(int32_t*)(pScr);
    scol=srow+3*L;

    for (l=0; l<3*L; l++) srow[l]=scol[l]=0;
    /*----------------------*/
    /* 1-st stage           */
    /*----------------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {
        float32_t t,maxval;
        maxval=fabsf(z[l+L*(0*3+0)]);                      srow[0*L+l]=0; scol[0*L+l]=0;
        t=fabsf(z[l+L*(0*3+1)]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=1; }
        t=fabsf(z[l+L*(0*3+2)]);  if(t>maxval) { maxval=t; srow[0*L+l]=0; scol[0*L+l]=2; }
        t=fabsf(z[l+L*(1*3+0)]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=0; }
        t=fabsf(z[l+L*(1*3+1)]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=1; }
        t=fabsf(z[l+L*(1*3+2)]);  if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=2; }
        t=fabsf(z[l+L*(2*3+0)]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=0; }
        t=fabsf(z[l+L*(2*3+1)]);  if(t>maxval) { maxval=t; srow[0*L+l]=2; scol[0*L+l]=1; }
        t=fabsf(z[l+L*(2*3+2)]);  if(t>maxval) {           srow[0*L+l]=2; scol[0*L+l]=2; }
    }
    for (l=0; l<L; l++)
    {
        float32_t t;
        int prow=srow[0*L+l];
        int pcol=scol[0*L+l];
        /* swap rows */
        t=z[l+L*(prow*3+0)]; z[l+L*(prow*3+0)]=z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(prow*3+1)]; z[l+L*(prow*3+1)]=z[l+L*1]; z[l+L*1]=t;
        t=z[l+L*(prow*3+2)]; z[l+L*(prow*3+2)]=z[l+L*2]; z[l+L*2]=t;
        /* swap columns */
        t=z[l+L*(0+pcol)];z[l+L*(0+pcol)]=z[l+L*0];z[l+L*0]=t;
        t=z[l+L*(3+pcol)];z[l+L*(3+pcol)]=z[l+L*3];z[l+L*3]=t;
        t=z[l+L*(6+pcol)];z[l+L*(6+pcol)]=z[l+L*6];z[l+L*6]=t;
    }
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*0]=1.0f/z[l+L*0];

    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        float32_t t;
        t=z[l+L*0];
        z[l+L*1]*=t;
        z[l+L*2]*=t;
        t=z[l+L*3];
        z[l+L*3] = -z[l+L*0] * t; 
        z[l+L*4] -= z[l+L*1] * t; 
        z[l+L*5] -= z[l+L*2] * t; 
        t=z[l+L*6];
        z[l+L*6] = -z[l+L*0] * t; 
        z[l+L*7] -= z[l+L*1] * t; 
        z[l+L*8] -= z[l+L*2] * t; 
    }
    /*----------------------*/
    /* 2-nd stage           */
    /*----------------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {
        float32_t t,maxval;
        maxval=fabsf(z[l+L*4]);  srow[1*L+l]=1; scol[1*L+l]=1; 
        t=fabsf(z[l+L*5]);  if(t>maxval) { maxval=t; srow[1*L+l]=1; scol[1*L+l]=2; }
        t=fabsf(z[l+L*7]);  if(t>maxval) { maxval=t; srow[1*L+l]=2; scol[1*L+l]=1; }
        t=fabsf(z[l+L*8]);  if(t>maxval) { srow[1*L+l]=2; scol[1*L+l]=2; }
    }
    for (l=0; l<L; l++)
    {
        float32_t t;
        int prow=srow[1*L+l];
        int pcol=scol[1*L+l];
        /* swap rows */
        t=z[l+L*(prow*3+0)]; z[l+L*(prow*3+0)]=z[l+L*3]; z[l+L*3]=t;
        t=z[l+L*(prow*3+1)]; z[l+L*(prow*3+1)]=z[l+L*4]; z[l+L*4]=t;
        t=z[l+L*(prow*3+2)]; z[l+L*(prow*3+2)]=z[l+L*5]; z[l+L*5]=t;
        /* swap columns */
        t=z[l+L*(0+pcol)];z[l+L*(0+pcol)]=z[l+L*1];z[l+L*1]=t;
        t=z[l+L*(3+pcol)];z[l+L*(3+pcol)]=z[l+L*4];z[l+L*4]=t;
        t=z[l+L*(6+pcol)];z[l+L*(6+pcol)]=z[l+L*7];z[l+L*7]=t;
    }
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*4]=1.0f/z[l+L*4];

    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        float32_t t;
        t=z[l+L*4];
        z[l+L*3]*=t;
        z[l+L*5]*=t;

        t=z[l+L*1];
        z[l+L*0] -= z[l+L*3] * t; 
        z[l+L*1] = -z[l+L*4] * t; 
        z[l+L*2] -= z[l+L*5] * t; 

        t=z[l+L*7];
        z[l+L*6] -= z[l+L*3] * t; 
        z[l+L*7] = -z[l+L*4] * t; 
        z[l+L*8] -= z[l+L*5] * t; 
    }
    /*----------------------*/
    /* 3-rd stage           */
    /*----------------------*/
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*8]=1.0f/z[l+L*8];
    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++) 
    {
        float32_t t;
        t=z[l+L*8];
        z[l+L*6]*=t;
        z[l+L*7]*=t;
        t=z[l+L*2];
        z[l+L*0] -= z[l+L*6] * t; 
        z[l+L*1] -= z[l+L*7] * t; 
        z[l+L*2] = -z[l+L*8] * t; 
        t=z[l+L*5];
        z[l+L*3] -= z[l+L*6] * t; 
        z[l+L*4] -= z[l+L*7] * t; 
        z[l+L*5] = -z[l+L*8] * t; 
    }

    /* final reverse permulation of columns  */
    for (l=0; l<L; l++)
    {
        int prow,pcol;
        float32_t t;
        pcol=srow[1*L+l];
        prow=scol[1*L+l];
        /* swap rows */
        t=z[l+L*(prow*3+0)]; z[l+L*(prow*3+0)]=z[l+L*3]; z[l+L*3]=t;
        t=z[l+L*(prow*3+1)]; z[l+L*(prow*3+1)]=z[l+L*4]; z[l+L*4]=t;
        t=z[l+L*(prow*3+2)]; z[l+L*(prow*3+2)]=z[l+L*5]; z[l+L*5]=t;
        /* swap columns */
        t=z[l+L*(0+pcol)];z[l+L*(0+pcol)]=z[l+L*1];z[l+L*1]=t;
        t=z[l+L*(3+pcol)];z[l+L*(3+pcol)]=z[l+L*4];z[l+L*4]=t;
        t=z[l+L*(6+pcol)];z[l+L*(6+pcol)]=z[l+L*7];z[l+L*7]=t;
        pcol=srow[0*L+l];
        prow=scol[0*L+l];
        /* swap rows */
        t=z[l+L*(prow*3+0)]; z[l+L*(prow*3+0)]=z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(prow*3+1)]; z[l+L*(prow*3+1)]=z[l+L*1]; z[l+L*1]=t;
        t=z[l+L*(prow*3+2)]; z[l+L*(prow*3+2)]=z[l+L*2]; z[l+L*2]=t;
        /* swap columns */
        t=z[l+L*(0+pcol)];z[l+L*(0+pcol)]=z[l+L*0];z[l+L*0]=t;
        t=z[l+L*(3+pcol)];z[l+L*(3+pcol)]=z[l+L*3];z[l+L*3]=t;
        t=z[l+L*(6+pcol)];z[l+L*(6+pcol)]=z[l+L*6];z[l+L*6]=t;
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pZrd;
          xb_vecN_2xf32 * restrict pZwr;

    int32_t *srow,*scol;
    xb_vecN_2xc16 * restrict pRow;
    xb_vecN_2xc16 * restrict pCol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 );

    if(L<=0 ) return;

    srow=(int32_t*)(pScr);
    scol=srow+3*L;
    /*----------------------*/
    /* 1-st stage           */
    /*----------------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    pZrd   =(const xb_vecN_2xf32*)(z);
    pRow=(xb_vecN_2xc16*)srow;
    pCol=(xb_vecN_2xc16*)scol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,t,maxval;
        vboolN_2 b1,b2,b3,b4,b5,b6,b7,b8;
        xb_vecN_2xc16 xrow=0,xcol=0;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
        maxval=BBE_ABSN_2XF32(z0);
        t=BBE_ABSN_2XF32(z1); b1=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b1); 
        t=BBE_ABSN_2XF32(z2); b2=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b2); 
        t=BBE_ABSN_2XF32(z3); b3=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b3); 
        t=BBE_ABSN_2XF32(z4); b4=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b4); 
        t=BBE_ABSN_2XF32(z5); b5=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b5); xrow=BBE_MOVN_2XC16T(1,xrow,b3|b4|b5);
        t=BBE_ABSN_2XF32(z6); b6=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(0,xcol,b6); 
        t=BBE_ABSN_2XF32(z7); b7=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b7); 
        t=BBE_ABSN_2XF32(z8); b8=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b8); xrow=BBE_MOVN_2XC16T(2,xrow,b6|b7|b8);
        BBE_SVN_2XC16_IP(xrow,pRow,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XC16_IP(xcol,pCol,2*BBE_SIMD_WIDTH);
    }
    // first permutation
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z);
    pZwr   =(      xb_vecN_2xf32*)(z);
    pRow=(xb_vecN_2xc16*)srow;
    pCol=(xb_vecN_2xc16*)scol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8,t;
        xb_vecN_2xc16 xrow,xcol;
        vboolN_2 brow1,brow2,bcol1,bcol2;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
        BBE_LVN_2XC16_IP(xrow,pRow,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XC16_IP(xcol,pCol,2*BBE_SIMD_WIDTH);

        brow1=BBE_EQN_2XC16(xrow,1);
        brow2=BBE_EQN_2XC16(xrow,2);
        bcol1=BBE_EQN_2XC16(xcol,1);
        bcol2=BBE_EQN_2XC16(xcol,2);

        t=z0;  z0=BBE_MOVN_2XF32T(z3,z0,brow1); z3=BBE_MOVN_2XF32T(t,z3,brow1);
        t=z0;  z0=BBE_MOVN_2XF32T(z6,z0,brow2); z6=BBE_MOVN_2XF32T(t,z6,brow2);
        t=z1;  z1=BBE_MOVN_2XF32T(z4,z1,brow1); z4=BBE_MOVN_2XF32T(t,z4,brow1);
        t=z1;  z1=BBE_MOVN_2XF32T(z7,z1,brow2); z7=BBE_MOVN_2XF32T(t,z7,brow2);
        t=z2;  z2=BBE_MOVN_2XF32T(z5,z2,brow1); z5=BBE_MOVN_2XF32T(t,z5,brow1);
        t=z2;  z2=BBE_MOVN_2XF32T(z8,z2,brow2); z8=BBE_MOVN_2XF32T(t,z8,brow2);

        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol1); z1=BBE_MOVN_2XF32T(t,z1,bcol1);
        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,bcol2); z2=BBE_MOVN_2XF32T(t,z2,bcol2);
        t=z3;  z3=BBE_MOVN_2XF32T(z4,z3,bcol1); z4=BBE_MOVN_2XF32T(t,z4,bcol1);
        t=z3;  z3=BBE_MOVN_2XF32T(z5,z3,bcol2); z5=BBE_MOVN_2XF32T(t,z5,bcol2);
        t=z6;  z6=BBE_MOVN_2XF32T(z7,z6,bcol1); z7=BBE_MOVN_2XF32T(t,z7,bcol1);
        t=z6;  z6=BBE_MOVN_2XF32T(z8,z6,bcol2); z8=BBE_MOVN_2XF32T(t,z8,bcol2);

        BBE_SVN_2XF32_XP(z0,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
    }
    /* reciprocal of main diagonal */
    matinvgjnxnsf_recip(z,L);
    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
        z1=BBE_MULN_2XF32(z1,z0);
        z2=BBE_MULN_2XF32(z2,z0);
        BBE_MULSN_2XF32(z4,z1,z3);
        BBE_MULSN_2XF32(z5,z2,z3);
        z3=BBE_NEGN_2XF32(BBE_MULN_2XF32(z0,z3));
        BBE_MULSN_2XF32(z7,z1,z6);
        BBE_MULSN_2XF32(z8,z2,z6);
        z6=BBE_NEGN_2XF32(BBE_MULN_2XF32(z0,z6));
        BBE_SVN_2XF32_XP(z0,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
    }
    /*----------------------*/
    /* 2-nd stage           */
    /*----------------------*/
    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+L*4);
    pRow=(xb_vecN_2xc16*)(srow+L);
    pCol=(xb_vecN_2xc16*)(scol+L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z4,z5,z7,z8,t,maxval;
        vboolN_2 b1,b2,b3;
        xb_vecN_2xc16 xrow=1,xcol=1;
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-4*L*sizeof(float32_t));
        maxval=BBE_ABSN_2XF32(z4);
        t=BBE_ABSN_2XF32(z5); b1=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b1); 
        t=BBE_ABSN_2XF32(z7); b2=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(1,xcol,b2); 
        t=BBE_ABSN_2XF32(z8); b3=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); xcol=BBE_MOVN_2XC16T(2,xcol,b3); xrow=BBE_MOVN_2XC16T(2,xrow,b3|b2); 
        BBE_SVN_2XC16_IP(xrow,pRow,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XC16_IP(xcol,pCol,2*BBE_SIMD_WIDTH);
    }
    // permutation of columns/rows
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z+L);
    pZwr   =(      xb_vecN_2xf32*)(z+L);
    pRow=(xb_vecN_2xc16*)(srow+L);
    pCol=(xb_vecN_2xc16*)(scol+L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z1,z2,z3,z4,z5,z6,z7,z8,t;
        xb_vecN_2xc16 xrow,xcol;
        vboolN_2 brow2,bcol2;
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-7*L*sizeof(float32_t));
        BBE_LVN_2XC16_IP(xrow,pRow,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XC16_IP(xcol,pCol,2*BBE_SIMD_WIDTH);

        brow2=BBE_EQN_2XC16(xrow,2);
        bcol2=BBE_EQN_2XC16(xcol,2);

        t=z3;  z3=BBE_MOVN_2XF32T(z6,z3,brow2); z6=BBE_MOVN_2XF32T(t,z6,brow2);
        t=z4;  z4=BBE_MOVN_2XF32T(z7,z4,brow2); z7=BBE_MOVN_2XF32T(t,z7,brow2);
        t=z5;  z5=BBE_MOVN_2XF32T(z8,z5,brow2); z8=BBE_MOVN_2XF32T(t,z8,brow2);

        t=z1;  z1=BBE_MOVN_2XF32T(z2,z1,bcol2); z2=BBE_MOVN_2XF32T(t,z2,bcol2);
        t=z4;  z4=BBE_MOVN_2XF32T(z5,z4,bcol2); z5=BBE_MOVN_2XF32T(t,z5,bcol2);
        t=z7;  z7=BBE_MOVN_2XF32T(z8,z7,bcol2); z8=BBE_MOVN_2XF32T(t,z8,bcol2);

        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-7*L*sizeof(float32_t));
    }
    /* reciprocal of main diagonal */
    matinvgjnxnsf_recip(z+L*4,L);
    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    __Pragma("no_reorder")
#if 0
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
        z3=BBE_MULN_2XF32(z3,z4);                         
        z5=BBE_MULN_2XF32(z5,z4);                      
        BBE_MULSN_2XF32(z0,z3,z1);                     
        BBE_MULSN_2XF32(z2,z5,z1);                     
        z1=BBE_NEGN_2XF32(BBE_MULN_2XF32(z4,z1));
        BBE_MULSN_2XF32(z6,z3,z7);
        BBE_MULSN_2XF32(z8,z5,z7);
        z7=BBE_NEGN_2XF32(BBE_MULN_2XF32(z4,z7));
        BBE_SVN_2XF32_XP(z0,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
    }
#else
        // use read order 4 3 5 1 0 2 7 6 8 for better scheduling 
    pZrd   =(const xb_vecN_2xf32*)(z+4*L);
    pZwr   =(      xb_vecN_2xf32*)(z+4*L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8;
        BBE_LVN_2XF32_XP(z4,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,-4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z0,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd, 5*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,-1*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd, 2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-4*L*sizeof(float32_t));
        z3=BBE_MULN_2XF32(z3,z4);                         
        z5=BBE_MULN_2XF32(z5,z4);                      
        BBE_MULSN_2XF32(z0,z3,z1);                     
        BBE_MULSN_2XF32(z2,z5,z1);                     
        z1=BBE_NEGN_2XF32(BBE_MULN_2XF32(z4,z1));
        BBE_MULSN_2XF32(z6,z3,z7);
        BBE_MULSN_2XF32(z8,z5,z7);
        z7=BBE_NEGN_2XF32(BBE_MULN_2XF32(z4,z7));
        BBE_SVN_2XF32_XP(z4,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,-4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z0,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr, 5*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,-1*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-4*L*sizeof(float32_t));
    }
#endif

    /*----------------------*/
    /* 3-rd stage           */
    /*----------------------*/
    /* reciprocal of main diagonal */
    matinvgjnxnsf_recip(z+L*8,L);
    /* process pivot row & GS elimination : for all rows excluding the pivot one */
#if 0
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
        z6=BBE_MULN_2XF32(z6,z8);
        z7=BBE_MULN_2XF32(z7,z8);
        BBE_MULSN_2XF32(z0,z6,z2);
        BBE_MULSN_2XF32(z1,z7,z2);
        z2=BBE_NEGN_2XF32(BBE_MULN_2XF32(z8,z2));
        BBE_MULSN_2XF32(z3,z6,z5);
        BBE_MULSN_2XF32(z4,z7,z5);
        z5=BBE_NEGN_2XF32(BBE_MULN_2XF32(z8,z5));
        BBE_SVN_2XF32_XP(z0,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
    }
#else
    // another read/write oreder (8 6 7 2 0 1 5 3 4) for better scheduling
    pZrd   =(const xb_vecN_2xf32*)(z+8*L);
    pZwr   =(      xb_vecN_2xf32*)(z+8*L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,z4,z5,z6,z7,z8;
        BBE_LVN_2XF32_XP(z8,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,-5*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z0,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,-2*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,   L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,2*BBE_SIMD_WIDTH+4*L*sizeof(float32_t));
        z6=BBE_MULN_2XF32(z6,z8);            
        z7=BBE_MULN_2XF32(z7,z8);
        BBE_MULSN_2XF32(z0,z6,z2);
        BBE_MULSN_2XF32(z1,z7,z2);
        z2=BBE_NEGN_2XF32(BBE_MULN_2XF32(z8,z2));
        BBE_MULSN_2XF32(z3,z6,z5);
        BBE_MULSN_2XF32(z4,z7,z5);
        z5=BBE_NEGN_2XF32(BBE_MULN_2XF32(z8,z5));
        BBE_SVN_2XF32_XP(z8,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,-5*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z0,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,-2*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,   L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,2*BBE_SIMD_WIDTH+4*L*sizeof(float32_t));
    }
#endif

    /* final reverse permulation of columns/rows  */
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    pRow=(xb_vecN_2xc16*)(srow);
    pCol=(xb_vecN_2xc16*)(scol);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 t,z0,z1,z2,z3,z4,z5,z6,z7,z8;
        xb_vecN_2xc16 xrow,xcol;
        vboolN_2 brow1,brow2,bcol1,bcol2;
        BBE_LVN_2XF32_XP(z0,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z1,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z2,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z3,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z4,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z5,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z6,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z7,pZrd,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(z8,pZrd,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));

        xcol=BBE_LVN_2XC16_X(pRow,L*sizeof(int32_t));
        xrow=BBE_LVN_2XC16_X(pCol,L*sizeof(int32_t));

        brow2=BBE_EQN_2XC16(xrow,2);
        bcol2=BBE_EQN_2XC16(xcol,2);

        t=z1;  z1=BBE_MOVN_2XF32T(z2,z1,bcol2); z2=BBE_MOVN_2XF32T(t,z2,bcol2);
        t=z4;  z4=BBE_MOVN_2XF32T(z5,z4,bcol2); z5=BBE_MOVN_2XF32T(t,z5,bcol2);
        t=z7;  z7=BBE_MOVN_2XF32T(z8,z7,bcol2); z8=BBE_MOVN_2XF32T(t,z8,bcol2);
        t=z3;  z3=BBE_MOVN_2XF32T(z6,z3,brow2); z6=BBE_MOVN_2XF32T(t,z6,brow2);
        t=z4;  z4=BBE_MOVN_2XF32T(z7,z4,brow2); z7=BBE_MOVN_2XF32T(t,z7,brow2);
        t=z5;  z5=BBE_MOVN_2XF32T(z8,z5,brow2); z8=BBE_MOVN_2XF32T(t,z8,brow2);

        BBE_LVN_2XC16_IP(xcol,pRow,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XC16_IP(xrow,pCol,2*BBE_SIMD_WIDTH);

        brow1=BBE_EQN_2XC16(xrow,1);
        brow2=BBE_EQN_2XC16(xrow,2);
        bcol1=BBE_EQN_2XC16(xcol,1);
        bcol2=BBE_EQN_2XC16(xcol,2);

        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol1); z1=BBE_MOVN_2XF32T(t,z1,bcol1);
        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,bcol2); z2=BBE_MOVN_2XF32T(t,z2,bcol2);
        t=z3;  z3=BBE_MOVN_2XF32T(z4,z3,bcol1); z4=BBE_MOVN_2XF32T(t,z4,bcol1);
        t=z3;  z3=BBE_MOVN_2XF32T(z5,z3,bcol2); z5=BBE_MOVN_2XF32T(t,z5,bcol2);
        t=z6;  z6=BBE_MOVN_2XF32T(z7,z6,bcol1); z7=BBE_MOVN_2XF32T(t,z7,bcol1);
        t=z6;  z6=BBE_MOVN_2XF32T(z8,z6,bcol2); z8=BBE_MOVN_2XF32T(t,z8,bcol2);
        t=z0;  z0=BBE_MOVN_2XF32T(z3,z0,brow1); z3=BBE_MOVN_2XF32T(t,z3,brow1);
        t=z0;  z0=BBE_MOVN_2XF32T(z6,z0,brow2); z6=BBE_MOVN_2XF32T(t,z6,brow2);
        t=z1;  z1=BBE_MOVN_2XF32T(z4,z1,brow1); z4=BBE_MOVN_2XF32T(t,z4,brow1);
        t=z1;  z1=BBE_MOVN_2XF32T(z7,z1,brow2); z7=BBE_MOVN_2XF32T(t,z7,brow2);
        t=z2;  z2=BBE_MOVN_2XF32T(z5,z2,brow1); z5=BBE_MOVN_2XF32T(t,z5,brow1);
        t=z2;  z2=BBE_MOVN_2XF32T(z8,z2,brow2); z8=BBE_MOVN_2XF32T(t,z8,brow2);


        BBE_SVN_2XF32_XP(z0,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z1,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z2,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z3,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z4,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z5,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z6,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z7,pZwr,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(z8,pZwr,2*BBE_SIMD_WIDTH-8*L*sizeof(float32_t));
    }
}
#endif

size_t matinvgj3x3sf_getScratchSize (int N, int L  )
{
    NASSERT(N==3);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int32_t)*6*L;
}
#else
DISCARD_FUN(void, matinvgj3x3sf,( void* pScr, float32_t * restrict z, int L  ))

size_t matinvgj3x3sf_getScratchSize (int N, int L  )
{
    NASSERT(N==3);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

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
#include "common.h"

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

#undef BBE_MOVN_2XF32T
#define BBE_MOVN_2XF32T(x,y,b) (BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(x),BBE_MOVNX16_FROMN_2XF32(y),BBE_MOVN_FROMN_2(b))))

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
void matinvgj2x2sf ( void* pScr, float32_t * restrict z, int L  )
#if 0
{
    float32_t t;
    float32_t *maxVal;
    int32_t *srow,*scol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 );

    if(L<=0 ) return;

    maxVal=(float32_t*)pScr;
    srow=(int32_t*)(maxVal+L);
    scol=srow+2*L;

    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {
        maxVal[l]=fabsf(z[l+L*0]);  
        srow[0*L+l]=0; scol[0*L+l]=0; 
        t=fabsf(z[l+L*1]);  
        if(t>maxVal[l]) { maxVal[l]=t; scol[0*L+l]=1; }

        t=fabsf(z[l+L*2]);  
        if(t>maxVal[l]) { maxVal[l]=t; srow[0*L+l]=1; scol[0*L+l]=0; }

        t=fabsf(z[l+L*3]);  
        if(t>maxVal[l]) { maxVal[l]=t; srow[0*L+l]=1; scol[0*L+l]=1; }
    }
    /* swap rows/columns */
    for (l=0; l<L; l++)
    {
        int prow=srow[0*L+l];
        int pcol=scol[0*L+l];
        t=z[l+L*(prow*2+0)];  z[l+L*(prow*2+0)]=z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(prow*2+1)];  z[l+L*(prow*2+1)]=z[l+L*1]; z[l+L*1]=t;
        t=z[l+L*(0+pcol)];    z[l+L*(0+pcol)]  =z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(2+pcol)];    z[l+L*(2+pcol)]  =z[l+L*2]; z[l+L*2]=t;
    }
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) maxVal[l]=1.0f/z[l+L*0];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++)
    {
        float32_t t;
        z[l+L*0] =maxVal[l];
        z[l+L*1]*=maxVal[l];
        t=z[l+L*2];
        z[l+L*2] = -z[l+L*0] * t; 
        z[l+L*3] -= z[l+L*1] * t; 
    }
    // last row
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) maxVal[l]=1.0f/z[l+L*3];

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++)
    {
        float32_t t;
        z[l+L*2]*=maxVal[l];
        z[l+L*3] =maxVal[l];
        t=z[l+L*1];
        z[l+L*0] -= z[l+L*2] * t; 
        z[l+L*1] = -z[l+L*3] * t; 
    }
    
    /* final reverse permulation of raws/columns  */
    for (l=0; l<L; l++)
    {
        int prow=scol[0*L+l];
        int pcol=srow[0*L+l];
        t=z[l+L*(prow*2+0)];  z[l+L*(prow*2+0)]=z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(prow*2+1)];  z[l+L*(prow*2+1)]=z[l+L*1]; z[l+L*1]=t;
        t=z[l+L*(0+pcol)];    z[l+L*(0+pcol)]  =z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(2+pcol)];    z[l+L*(2+pcol)]  =z[l+L*2]; z[l+L*2]=t;
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pZrd;
          xb_vecN_2xf32 * restrict pZwr;
    vboolN_2 * browcol;  /* pairs of row/col boolean flags for permutations, [2*L]*/
    vboolN_2 * restrict bRowCol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 );

    if(L<=0 ) return;

    browcol=(vboolN_2*)pScr;

    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    pZrd   =(const xb_vecN_2xf32*)z;
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 b1,b2,b3,brow,bcol;
        xb_vecN_2xf32 z0,z1,z2,z3,maxval,t;
        z1=BBE_LVN_2XF32_X(pZrd,1*L*sizeof(float32_t));
        z2=BBE_LVN_2XF32_X(pZrd,2*L*sizeof(float32_t));
        z3=BBE_LVN_2XF32_X(pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_IP(z0,pZrd,2*BBE_SIMD_WIDTH);

        maxval=BBE_ABSN_2XF32(z0);
        t=BBE_ABSN_2XF32(z1); b1=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXN_2XF32(t,maxval); 
        t=BBE_ABSN_2XF32(z2); b2=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXN_2XF32(t,maxval); 
        t=BBE_ABSN_2XF32(z3); b3=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXN_2XF32(t,maxval); 

        brow=b2;
        bcol=b1&~b2;
        brow|=b3;
        bcol|=b3;

        BBE_SBN_2_IP(brow,bRowCol,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol,bRowCol,sizeof(vboolN_2));
    }
    __Pragma("no_reorder")

    /* swap rows/columns */
    pZrd   =(const xb_vecN_2xf32*)z;
    pZwr   =(      xb_vecN_2xf32*)z;
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        vboolN_2 brow,bcol;
        xb_vecN_2xf32 t,z0,z1,z2,z3;
        BBE_LBN_2_IP(brow,bRowCol,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol,bRowCol,sizeof(vboolN_2));
        z1=BBE_LVN_2XF32_X(pZrd,1*L*sizeof(float32_t));
        z2=BBE_LVN_2XF32_X(pZrd,2*L*sizeof(float32_t));
        z3=BBE_LVN_2XF32_X(pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_IP(z0,pZrd,2*BBE_SIMD_WIDTH);

        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,brow); z2=BBE_MOVN_2XF32T(t,z2,brow);
        t=z1;  z1=BBE_MOVN_2XF32T(z3,z1,brow); z3=BBE_MOVN_2XF32T(t,z3,brow);
        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol); z1=BBE_MOVN_2XF32T(t,z1,bcol);
        t=z2;  z2=BBE_MOVN_2XF32T(z3,z2,bcol); z3=BBE_MOVN_2XF32T(t,z3,bcol);
        z0=BBE_RECIPN_2XF32(z0);

        BBE_SVN_2XF32_X (z1,pZwr,1*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z2,pZwr,2*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z3,pZwr,3*L*sizeof(float32_t));
        BBE_SVN_2XF32_IP(z0,pZwr,2*BBE_SIMD_WIDTH);
    }

    /* process pivot row &  GS elimination : for all rows excluding the pivot one */
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_2xf32*)(z);
    pZwr   =(      xb_vecN_2xf32*)(z);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3;
        z1=BBE_LVN_2XF32_X(pZrd,1*L*sizeof(float32_t));
        z2=BBE_LVN_2XF32_X(pZrd,2*L*sizeof(float32_t));
        z3=BBE_LVN_2XF32_X(pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_IP(z0,pZrd,2*BBE_SIMD_WIDTH);

        z1=BBE_MULN_2XF32(z1,z0);
        BBE_MULSN_2XF32(z3,z1,z2);
        z2=BBE_NEGN_2XF32(BBE_MULN_2XF32(z0,z2));

        BBE_SVN_2XF32_X (z1,pZwr,1*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z2,pZwr,2*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z3,pZwr,3*L*sizeof(float32_t));
        BBE_SVN_2XF32_IP(z0,pZwr,2*BBE_SIMD_WIDTH);
    }
    // last row
    __Pragma("no_reorder")
    /* reciprocal of main diagonal */
    pZrd   =(const xb_vecN_2xf32*)(z+L*3);
    pZwr   =(      xb_vecN_2xf32*)(z+L*3);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 x;
        BBE_LVN_2XF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
        x=BBE_RECIPN_2XF32(x);
        BBE_SVN_2XF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
    }

    /* process pivot row &  GS elimination : for all rows excluding the pivot one */
    __Pragma("no_reorder")
    /* reciprocal of main diagonal */
    pZrd   =(const xb_vecN_2xf32*)(z);
    pZwr   =(      xb_vecN_2xf32*)(z);
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 t,z0,z1,z2,z3;
        vboolN_2 brow,bcol;
        z1=BBE_LVN_2XF32_X(pZrd,1*L*sizeof(float32_t));
        z2=BBE_LVN_2XF32_X(pZrd,2*L*sizeof(float32_t));
        z3=BBE_LVN_2XF32_X(pZrd,3*L*sizeof(float32_t));
        BBE_LVN_2XF32_IP(z0,pZrd,2*BBE_SIMD_WIDTH);
        z2=BBE_MULN_2XF32(z2,z3);
        BBE_MULSN_2XF32(z0,z2,z1);
        z1=BBE_NEGN_2XF32(BBE_MULN_2XF32(z3,z1));
        /* final reverse permulation of rows/columns  */
        BBE_LBN_2_IP(bcol,bRowCol,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow,bRowCol,sizeof(vboolN_2));

        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,brow); z2=BBE_MOVN_2XF32T(t,z2,brow);
        t=z1;  z1=BBE_MOVN_2XF32T(z3,z1,brow); z3=BBE_MOVN_2XF32T(t,z3,brow);
        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol); z1=BBE_MOVN_2XF32T(t,z1,bcol);
        t=z2;  z2=BBE_MOVN_2XF32T(z3,z2,bcol); z3=BBE_MOVN_2XF32T(t,z3,bcol);

        BBE_SVN_2XF32_X (z1,pZwr,1*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z2,pZwr,2*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (z3,pZwr,3*L*sizeof(float32_t));
        BBE_SVN_2XF32_IP(z0,pZwr,2*BBE_SIMD_WIDTH);
    }
}
#endif

size_t matinvgj2x2sf_getScratchSize (int N, int L  )
{
    NASSERT(N==2);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N!=2) return 0;
    (void)L; (void)N;
    return sizeof(vboolN_2)*2*(L>>(LOG2_BBE_SIMD_WIDTH-1));
}
#else
DISCARD_FUN(void, matinvgj2x2sf,( void* pScr, float32_t * restrict z, int L  ))

size_t matinvgj2x2sf_getScratchSize (int N, int L  )
{
    NASSERT(N==2);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

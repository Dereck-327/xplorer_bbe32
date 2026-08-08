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
#if 0
#include <math.h>
#include <float.h>
#include <complex.h>
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
void cmatinvgj2x2sf ( void* pScr, complex_float * restrict z, int L  )
#if 0
{
    int32_t *srow,*scol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 );

    if(L<=0 ) return;

    srow=(int32_t*)(pScr);
    scol=srow+2*L;

    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (l=0; l<L; l++)
    {
        float32_t t;
        float32_t maxval;       
        maxval=sqrc(z[l+L*0]);          srow[0*L+l]=0; scol[0*L+l]=0; 
        t=sqrc(z[l+L*1]);          if(t>maxval) { maxval=t; scol[0*L+l]=1; }
        t=sqrc(z[l+L*2]);          if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=0; }
        t=sqrc(z[l+L*3]);          if(t>maxval) { maxval=t; srow[0*L+l]=1; scol[0*L+l]=1; }
    }
    /* swap rows/columns */
    for (l=0; l<L; l++)
    {
        complex_float t;
        int prow=srow[0*L+l];
        int pcol=scol[0*L+l];
        t=z[l+L*(prow*2+0)];  z[l+L*(prow*2+0)]=z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(prow*2+1)];  z[l+L*(prow*2+1)]=z[l+L*1]; z[l+L*1]=t;
        t=z[l+L*(0+pcol)];    z[l+L*(0+pcol)]  =z[l+L*0]; z[l+L*0]=t;
        t=z[l+L*(2+pcol)];    z[l+L*(2+pcol)]  =z[l+L*2]; z[l+L*2]=t;
    }
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*0]=recipc(z[l+L*0]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++)
    {
        z[l+L*1] =                mulc(z[l+L*1] , z[l+L*0]);
        z[l+L*3] = subc(z[l+L*3], mulc(z[l+L*1] , z[l+L*2])); 
        z[l+L*2] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*0] , z[l+L*2])); 
    }
    // last row
    /* reciprocal of main diagonal */
    for (l=0; l<L; l++) z[l+L*3]=recipc(z[l+L*3]);

    /* process pivot row */
    /* GS elimination : for all rows excluding the pivot one */
    for(l=0; l<L; l++)
    {
        z[l+L*2] = mulc(z[l+L*2] , z[l+L*3]);
        z[l+L*0] = subc(z[l+L*0], mulc(z[l+L*2] , z[l+L*1]));
        z[l+L*1] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*3] , z[l+L*1])); 
    }
    
    /* final reverse permulation of raws/columns  */
    for (l=0; l<L; l++)
    {
        complex_float t;
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
    const xb_vecN_4xcf32 * restrict pZrd;
          xb_vecN_4xcf32 * restrict pZwr;
    vboolN_2 * browcol;  /* pairs of row/col boolean flags for permutations, [6*L>>((LOG2_BBE_SIMD_WIDTH-2))/]*/
    vboolN_2 * restrict bRowCol;
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 );

    if(L<=0 ) return;

    browcol=(vboolN_2 *)(pScr);

    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    pZrd   =(const xb_vecN_4xcf32*)(z);
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 a,z0,z1,z2,z3;
        xb_vecN_2xf32 t,maxval;
        vboolN_2 b1,b2,b3;
        vboolN_2 brow1,bcol1;
        BBE_LVN_4XCF32_XP(z0,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2,pZrd,  L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3,pZrd,2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));
        a=BBE_MULMN_4XCF32(z0 ,z0 ,0,0); BBE_MULMASN_4XCF32(a,z0 ,z0 ,0,15); maxval=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a));       
        a=BBE_MULMN_4XCF32(z1 ,z1 ,0,0); BBE_MULMASN_4XCF32(a,z1 ,z1 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b1=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        a=BBE_MULMN_4XCF32(z2 ,z2 ,0,0); BBE_MULMASN_4XCF32(a,z2 ,z2 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b2=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 
        a=BBE_MULMN_4XCF32(z3 ,z3 ,0,0); BBE_MULMASN_4XCF32(a,z3 ,z3 ,0,15); t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a)); b3=BBE_OGTN_2XF32(t,maxval); maxval=BBE_MAXNUMN_2XF32(t,maxval); 

        brow1=b2;
        bcol1=b1&~b2;
        brow1|=b3;
        bcol1|=b3;

        BBE_SBN_2_IP(brow1,bRowCol,sizeof(vboolN_2));
        BBE_SBN_2_IP(bcol1,bRowCol,sizeof(vboolN_2));
    }

    /* swap rows/columns */
    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z);
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,t;
        vboolN_2 brow2,bcol2;
        BBE_LVN_2XF32_XP(z0,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z1,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));

        BBE_LBN_2_IP(brow2,bRowCol,sizeof(vboolN_2));
        BBE_LBN_2_IP(bcol2,bRowCol,sizeof(vboolN_2));

        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,brow2); z2=BBE_MOVN_2XF32T(t,z2,brow2);
        t=z1;  z1=BBE_MOVN_2XF32T(z3,z1,brow2); z3=BBE_MOVN_2XF32T(t,z3,brow2);
        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol2); z1=BBE_MOVN_2XF32T(t,z1,bcol2);
        t=z2;  z2=BBE_MOVN_2XF32T(z3,z2,bcol2); z3=BBE_MOVN_2XF32T(t,z3,bcol2);

        BBE_SVN_2XF32_XP(z0,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z1,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));
    }

    /* reciprocal of main diagonal */
    cmatinvgjnxnsf_recip(z,L);
    /* process pivot row & GS elimination : for all rows excluding the pivot one */
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z+1*L);
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z0,z1,z2,z3,t;
        BBE_LVN_4XCF32_XP(z0,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z1,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z2,pZrd,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(z3,pZrd,2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));

        z1=BBE_MULN_4XCF32(z1,z0);
        BBE_MULSN_4XCF32(z3,z1,z2);
        t=z2 ; z2 =BBE_MULMN_4XCF32 (z0,t, 3, 4); BBE_MULMASN_4XCF32 (z2 , z0,t, 2, 11);
    z3=BBE_RECIPN_4XCF32(z3);

        BBE_SVN_4XCF32_XP(z1,pZwr,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z2,pZwr,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(z3,pZwr,2*BBE_SIMD_WIDTH-2*L*sizeof(complex_float));
    }
    // last row
    /* reciprocal of main diagonal */

    __Pragma("no_reorder")
    pZrd   =(const xb_vecN_4xcf32*)(z);
    pZwr   =(      xb_vecN_4xcf32*)(z);
    bRowCol=browcol;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 z0,z1,z2,z3,t;
        vboolN_2 brow2,bcol2;
        BBE_LVN_2XF32_XP(z0,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z1,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z2,castxcc(xb_vecN_2xf32,pZrd),L*sizeof(complex_float));
        BBE_LVN_2XF32_XP(z3,castxcc(xb_vecN_2xf32,pZrd),2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));

        /* GS elimination : for all rows excluding the pivot one */
        t=BBE_MULMN_2XF32 (z2,z3, 0, 4); BBE_MULMASN_2XF32(t, z2,z3, 1, 11); z2=t;
        BBE_MULMASN_2XF32(z0,z2,z1, 3,  4);
        BBE_MULMASN_2XF32(z0,z2,z1, 2, 11);
        t=z1 ; z1 =BBE_MULMN_2XF32 (z3,t, 3, 4); BBE_MULMASN_2XF32(z1 , z3,t, 2, 11);

        /* final reverse permulation of rows/columns  */
        BBE_LBN_2_IP(bcol2,bRowCol,sizeof(vboolN_2));
        BBE_LBN_2_IP(brow2,bRowCol,sizeof(vboolN_2));

        t=z0;  z0=BBE_MOVN_2XF32T(z1,z0,bcol2); z1=BBE_MOVN_2XF32T(t,z1,bcol2);
        t=z2;  z2=BBE_MOVN_2XF32T(z3,z2,bcol2); z3=BBE_MOVN_2XF32T(t,z3,bcol2);
        t=z0;  z0=BBE_MOVN_2XF32T(z2,z0,brow2); z2=BBE_MOVN_2XF32T(t,z2,brow2);
        t=z1;  z1=BBE_MOVN_2XF32T(z3,z1,brow2); z3=BBE_MOVN_2XF32T(t,z3,brow2);

        BBE_SVN_2XF32_XP(z0,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z1,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z2,castxcc(xb_vecN_2xf32,pZwr),L*sizeof(complex_float));
        BBE_SVN_2XF32_XP(z3,castxcc(xb_vecN_2xf32,pZwr),2*BBE_SIMD_WIDTH-3*L*sizeof(complex_float));
    }

}
#endif

size_t cmatinvgj2x2sf_getScratchSize (int N, int L  )
{
    NASSERT(N==2);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return //sizeof(int32_t)*2*2*L+
            sizeof(vboolN_2)*2*(L>>(LOG2_BBE_SIMD_WIDTH-2)); /* boolean masks */
}
#else
DISCARD_FUN(void, cmatinvgj2x2sf ,( void* pScr, complex_float * restrict z, int L  ))

size_t cmatinvgj2x2sf_getScratchSize (int N, int L  )
{
    NASSERT(N==2);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

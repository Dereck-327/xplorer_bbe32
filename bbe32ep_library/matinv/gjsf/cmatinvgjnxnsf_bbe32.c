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

static void swapRows(complex_float * z, const int32_t* srow, int N,int k,int L)
#if 0
{
    int l,n;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(srow,BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    for (l=0; l<L; l++)
    {
        int prow=srow[l];
        /* swap rows */
        for(n=0; n<N; n++) 
        {
            complex_float t;
            t=z[l+L*(prow*N+n)];
            z[l+L*(prow*N+n)]=z[l+L*(k*N+n)];
            z[l+L*(k*N+n   )]=t;
        }
    }
}
#else
{
    int n,l;
    xb_vecN_4xcf32* restrict pAw;
    xtcomplexfloat* restrict pA0;
    xtcomplexfloat* restrict pA1;
    xtcomplexfloat* restrict pA2;
    xtcomplexfloat* restrict pA3;
    valign aRow;
    const xb_vecNx16 * restrict pRow;
    pRow=(const xb_vecNx16 *)srow;
    aRow=BBE_LA_PP(pRow);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(srow,BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);

    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        xb_vecNx40 w;
        xb_vecNx16 xrow,addr,al,ah; 
        BBE_LAVNX16_XP(xrow,aRow,pRow,BBE_SIMD_WIDTH);
        xrow=BBE_SELNX16I(xrow,xrow,BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        w=BBE_SEQNX40();
        BBE_MULANX16(w,xrow,L*N);
        w=BBE_SLLINX40(w,3);
        w=BBE_ADDNX40(w,(int32_t)(z+l));
        al=BBE_PACKLNX40(w);
        ah=BBE_PACKVNX40(w,16);
        addr=BBE_SELNX16I(ah,al,BBE_SELI_INTERLEAVE_1_LO);
        pA0=(xtcomplexfloat *)BBE_EXTRNX16C(addr,0);
        pA1=(xtcomplexfloat *)BBE_EXTRNX16C(addr,1);
        pA2=(xtcomplexfloat *)BBE_EXTRNX16C(addr,2);
        pA3=(xtcomplexfloat *)BBE_EXTRNX16C(addr,3);
        pAw=(xb_vecN_4xcf32*)&z[l+L*k*N];
        for (n = 0; n<N; n++)
        {
            xb_vecN_4xcf32 t0,t1,t2,t3,T0,T1;
            T1=BBE_LVN_4XCF32_I(pAw,0);
            t0=BBE_LSN_4XCF32_I(pA0,0);
            t1=BBE_LSN_4XCF32_I(pA1,0);
            t2=BBE_LSN_4XCF32_I(pA2,0);
            t3=BBE_LSN_4XCF32_I(pA3,0);
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T0=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            BBE_SVN_4XCF32_XP(T0,pAw,L*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(T1,pA0,L*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_4 ),pA1,L*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_8 ),pA2,L*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_12),pA3,L*sizeof(complex_float));
        }
    }
}
#endif
/* swap k-th column with columns given in scol[L] */
static void swapCols(complex_float * z, const int32_t* scol, int N,int k,int L)
#if 0
{
    int l,n;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(scol,BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    for (l=0; l<L; l++)
    {
        int pcol;
        pcol=scol[l];
        for(n=0; n<N; n++) 
        {
            complex_float t;
            t=z[l+L*(n*N+pcol)];
            z[l+L*(n*N+pcol)]=z[l+L*(n*N+k)];
            z[l+L*(n*N+k   )]=t;
        }
    }
}
#else
{
    int n,l;
    xb_vecN_4xcf32* restrict pAw;
    xtcomplexfloat* restrict pA0;
    xtcomplexfloat* restrict pA1;
    xtcomplexfloat* restrict pA2;
    xtcomplexfloat* restrict pA3;
    valign aCol;
    const xb_vecNx16 * restrict pCol;
    pCol=(const xb_vecNx16 *)scol;
    aCol=BBE_LA_PP(pCol);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(scol,BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);

    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        xb_vecNx40 w;
        xb_vecNx16 xcol,addr,al,ah; 
        BBE_LAVNX16_XP(xcol,aCol,pCol,BBE_SIMD_WIDTH);
        xcol=BBE_SELNX16I(xcol,xcol,BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        w=BBE_SEQNX40();
        BBE_MULANX16(w,xcol,L);
        w=BBE_SLLINX40(w,3);
        w=BBE_ADDNX40(w,(int32_t)(z+l));
        al=BBE_PACKLNX40(w);
        ah=BBE_PACKVNX40(w,16);
        addr=BBE_SELNX16I(ah,al,BBE_SELI_INTERLEAVE_1_LO);
        pA0=(xtcomplexfloat *)BBE_EXTRNX16C(addr,0);
        pA1=(xtcomplexfloat *)BBE_EXTRNX16C(addr,1);
        pA2=(xtcomplexfloat *)BBE_EXTRNX16C(addr,2);
        pA3=(xtcomplexfloat *)BBE_EXTRNX16C(addr,3);
        pAw=(xb_vecN_4xcf32*)&z[l+L*k];
        for (n = 0; n<N; n++)
        {
            xb_vecN_4xcf32 t0,t1,t2,t3,T0,T1;
            T1=BBE_LVN_4XCF32_I(pAw,0);
            t0=BBE_LSN_4XCF32_I(pA0,0);
            t1=BBE_LSN_4XCF32_I(pA1,0);
            t2=BBE_LSN_4XCF32_I(pA2,0);
            t3=BBE_LSN_4XCF32_I(pA3,0);
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T0=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            BBE_SVN_4XCF32_XP(T0,pAw,L*N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(T1,pA0,L*N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_4 ),pA1,L*N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_8 ),pA2,L*N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(BBE_SELN_4XCF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_12),pA3,L*N*sizeof(complex_float));
        }
    }
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

void cmatinvgjnxnsf ( void* pScr, complex_float * restrict z, int N, int L  )
#if 0
{
    complex_float* maxVal;
    int32_t *srow,*scol;
    int l,k,m,n;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);

    if(L<=0 || N<=1) return;

    maxVal=(complex_float*)pScr;
    srow=(int32_t*)(maxVal+L);
    scol=srow+N*L;

    /* pivoting: search the absolute maximum and its position in the matrix 
        at all positions excluding previously used */
    for (k=0; k<N*L; k++) srow[k]=scol[k]=0;
    for (k=0; k<N; k++)
    {
        for (l=0; l<L; l++)
        {
            float32_t t,maxval=FLT_MIN;
            int prow,pcol;
            prow=k; pcol=k;
            for (m=k; m<N; m++)
            {
                for(n=k; n<N; n++) 
                {
                    t=sqrc(z[l+L*(m*N+n)]);  
                    if(t>maxval) { maxval=t; prow=m; pcol=n; }
                }
            }
            srow[k*L+l]=prow; scol[k*L+l]=pcol;
        }
        swapRows(z, srow+k*L, N,k,L);
        swapCols(z, scol+k*L, N,k,L);
        /* reciprocal of main diagonal */
        for (l=0; l<L; l++) maxVal[l]=recipc(z[l+L*(k*N+k)]);

        /* process pivot row */
        for(l=0; l<L; l++)
        {
            complex_float t;
            t=maxVal[l];
            for (n=0; n<N; n++) 
            {
                z[l+L*(k*N+n)]=mulc(z[l+L*(k*N+n)],t);
            }
        }
        for(l=0; l<L; l++) z[l+L*(k*N+k)]=maxVal[l];

            /* GS elimination : for all rows excluding the pivot one */
        for (m=0; m<N; m++)
        {
            if(m==k) continue; 
            for(l=0; l<L; l++)
            {
                complex_float t;
                t=z[l+L*(m*N+k)];
                maxVal[l]=t;
                for (n=0; n<N; n++)
                {
                     z[l+L*(m*N+n)] = subc(z[l+L*(m*N+n)],mulc(z[l+L*(k*N+n)] , t)); 
                }
            }
            for(l=0; l<L; l++) z[l+L*(m*N+k)] = subc(makecomplexf(0.f,0.f),mulc(z[l+L*(k*N+k)] , maxVal[l]));
        }
    }
    
    /* final reverse permulation of columns  */
    for (k=N-2; k>=0; k--)
    {
        swapCols(z, srow+k*L, N,k,L);
        swapRows(z, scol+k*L, N,k,L);
    }
}
#else
{
    const xb_vecN_4xcf32* restrict pZrd;
    const xb_vecN_4xcf32* restrict pZrw;
          xb_vecN_4xcf32* restrict pZwr;
    xb_vecN_2xc16 * restrict pRow;
    xb_vecN_2xc16 * restrict pCol;
          xb_vecN_4xcf32* restrict pMaxVal;
    complex_float* maxVal;
    int32_t *srow,*scol;
    int l,k,m,n;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);

    if(L<=0 || N<=1) return;

    maxVal=(complex_float*)pScr;
    srow=(int32_t*)(maxVal+L);
    scol=srow+N*L;

    for (k=0; k<N; k++)
    {
        if (k!=N-1)
        {
            /* pivoting: search the absolute maximum and its position in the matrix 
                at all positions excluding previously used */
            valign aRow,aCol;
            pRow=(xb_vecN_2xc16 *)(srow+k*L);
            pCol=(xb_vecN_2xc16 *)(scol+k*L);
            pZrd=(const xb_vecN_4xcf32*)(z+L*(k*N+k));
            aCol=aRow=BBE_ZALIGN();
            for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
            {
                xb_vecN_2xf32 maxval=FLT_MIN,t,a;
                xb_vecN_2xc16 prow,pcol,vm,vn;
                prow=k; pcol=k;
                vm=k;
                __Pragma("loop_count min=2")
                for (m=k; m<N; m++)
                {
                    vn=k;
                    __Pragma("loop_count min=2")
                    for(n=k; n<N; n++) 
                    {
                        vboolN_2 b;
                        xb_vecN_4xcf32 tt;
                        BBE_LVN_4XCF32_XP(tt,pZrd,L*sizeof(complex_float));
                        t=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tt));
                        a=BBE_MULMN_2XF32(t,t,0,0); 
                        BBE_MULMASN_2XF32(a,t,t,0,15); 
                        t=a;
                        b=BBE_OGTN_2XF32(t,maxval);
                        maxval=BBE_MAXNUMN_2XF32(t,maxval);
                        prow=BBE_MOVN_2XC16T(vm,prow,b);
                        pcol=BBE_MOVN_2XC16T(vn,pcol,b);
                        vn=BBE_ADDN_2XC16(vn,1);
                    }
                    vm=BBE_ADDN_2XC16(vm,1);
                    pZrd+= (k*L*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                }
                prow=BBE_SELN_2XC16I(prow,prow,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                pcol=BBE_SELN_2XC16I(pcol,pcol,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                BBE_SAVN_2XC16_XP(prow,aRow,pRow,BBE_SIMD_WIDTH);
                BBE_SAVN_2XC16_XP(pcol,aCol,pCol,BBE_SIMD_WIDTH);
                pZrd+= ((BBE_SIMD_WIDTH/4-N*(N-k)*L)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            BBE_SAN_2XC16POS_FP(aRow,pRow);
            BBE_SAN_2XC16POS_FP(aCol,pCol);
            swapRows(z, srow+k*L, N,k,L);
            swapCols(z, scol+k*L, N,k,L);
        }
        /* reciprocal of main diagonal */
        pMaxVal=(xb_vecN_4xcf32*)(maxVal);
        pZrd   =(const xb_vecN_4xcf32*)&z[L*(k*N+k)];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
        {
            xb_vecN_4xcf32 t;
            BBE_LVN_4XCF32_IP(t,pZrd,2*BBE_SIMD_WIDTH);
            t=BBE_RECIPN_4XCF32(t);
            BBE_SVN_4XCF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
        }

        /* process pivot row */
        __Pragma("no_reorder")
        pMaxVal=(xb_vecN_4xcf32*)(maxVal);
        pZrd   =(const xb_vecN_4xcf32*)&z[L*N*k];
        pZwr   =(      xb_vecN_4xcf32*)&z[L*N*k];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
        {
            xb_vecN_4xcf32 x,t;
            BBE_LVN_4XCF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
            __Pragma("loop_count min=2")
            for (n=0; n<N; n++) 
            {
                BBE_LVN_4XCF32_XP(x,pZrd,L*sizeof(complex_float));
                x=BBE_MULN_4XCF32(x,t);
                BBE_SVN_4XCF32_XP(x,pZwr,L*sizeof(complex_float));
            }
            pZrd+= ((BBE_SIMD_WIDTH/4-L*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            pZwr+= ((BBE_SIMD_WIDTH/4-L*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
        }
        __Pragma("no_reorder")
        pMaxVal=(xb_vecN_4xcf32*)(maxVal);
        pZwr   =(xb_vecN_4xcf32*)&z[L*(k*N+k)];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
        {
            xb_vecN_4xcf32 t;
            BBE_LVN_4XCF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(t,pZwr   ,2*BBE_SIMD_WIDTH);
        }

            /* GS elimination : for all rows excluding the pivot one */
        __Pragma("no_reorder")
        for (m=0; m<N; m++)
        {
            xb_vecN_4xcf32 x,y,t;
            if(m==k) continue; 
            pZrd=(const xb_vecN_4xcf32*)&z[L*N*k];
            pZwr=(      xb_vecN_4xcf32*)&z[L*N*m];
            pZrw=(const xb_vecN_4xcf32*)pZwr;
            __Pragma("loop_count min=1")
            for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
            {
                t = BBE_LVN_4XCF32_X  (pZrw,k*L*sizeof(complex_float));
                __Pragma("loop_count min=2")
                for (n=0; n<N; n++)
                {
                    BBE_LVN_4XCF32_XP (x,pZrd,L*sizeof(complex_float));
                    BBE_LVN_4XCF32_XP (y,pZrw,L*sizeof(complex_float));
                    BBE_MULSN_4XCF32(y,x,t);
                    BBE_SVN_4XCF32_XP (y,pZwr,L*sizeof(complex_float));
                }
                x = BBE_LVN_4XCF32_X  (pZrd,(k-N)*L*sizeof(complex_float));
                y=BBE_MULMN_4XCF32 (   x, t, 3, 4);
                BBE_MULMASN_4XCF32 (y, x, t, 2, 11);
                BBE_SVN_4XCF32_X  (y,pZwr,(k-N)*L*sizeof(complex_float));
                pZrd+= ((BBE_SIMD_WIDTH/4-N*L)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pZwr+= ((BBE_SIMD_WIDTH/4-N*L)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pZrw+= ((BBE_SIMD_WIDTH/4-N*L)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
           }
        }
    }
    
    /* final reverse permulation of columns  */
    for (k=N-2; k>=0; k--)
    {
        swapCols(z, srow+k*L, N,k,L);
        swapRows(z, scol+k*L, N,k,L);
    }
}
#endif
size_t cmatinvgjnxnsf_getScratchSize (int N, int L  )
{
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int32_t)*N*2*L+sizeof(complex_float)*L;
}

#else
DISCARD_FUN(void, cmatinvgjnxnsf ,( void* pScr, complex_float * restrict z, int N, int L  ))

size_t cmatinvgjnxnsf_getScratchSize (int N, int L  )
{
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

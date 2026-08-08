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
#include "cmatinvgjnxnnf_common.h"
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

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    m=XT_MIN(m,LOG2_BBE_SIMD_WIDTH-1);
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}
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

inline_ void cmatinvgjnxnnf_updPivot(complex_float *z, complex_float* recip, int N, int L, int SA)
#if 0
{
    int l,n;
    for (l=0; l<L; l++)
    {
        for (n=0; n<N; n++) z[SA*l+n]=mulc(z[SA*l+n],recip[l]); 
    }
}
#else
{
    const xtcomplexfloat* restrict pRecip;
    const xb_vecN_4xcf32* restrict pZrd;
          xb_vecN_4xcf32* restrict pZwr;
    int n,l,NN;

    NASSERT(N%4==0 && N>0);
    NASSERT_ALIGN32(z);
    NN = N>>(LOG2_BBE_SIMD_WIDTH-2);
    pRecip=(const xtcomplexfloat*)recip;
    pZrd=(const xb_vecN_4xcf32*)(z);
    pZwr=(      xb_vecN_4xcf32*)(z);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t,x;
        BBE_LSN_4XCF32_IP(t,pRecip,sizeof(complex_float));
        t=BBE_REPN_4XCF32(t,0);
        for (n=0; n<NN; n++) 
        {
            BBE_LVN_4XCF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
            x=BBE_MULN_4XCF32(x,t);
            BBE_SVN_4XCF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
        }
        pZrd+= ((SA-N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
        pZwr+= ((SA-N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
    }
}
#endif

/* GS elimination : for all rows excluding the pivot one */
inline_ void cmatinvgjnxnnf_elimGS(complex_float *z, complex_float *zk, const complex_float* recip, int N, int L, int SA)
#if 0
{
    int l,m,n;
    for (l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            complex_float t;
            t=recip[m*L+l];
            for (n=0; n<N; n++) z[SA*l+m*N+n]=subc(z[SA*l+m*N+n],mulc(zk[SA*l+n],t)); 
        }
    }
}
#else
{
    const xtcomplexfloat* restrict pRecip;
    const xb_vecN_4xcf32* restrict pZrd;
    const xb_vecN_4xcf32* restrict pZk;
          xb_vecN_4xcf32* restrict pZwr;
    int m,n,l,NN;

    NASSERT(N%4==0 && N>0);
    NASSERT_ALIGN32(z);
    NASSERT_ALIGN32(zk);
    NASSERT(N*N==SA);   // for SIMD_WIDTH<=16
    (void)SA;
    NN = N>>(LOG2_BBE_SIMD_WIDTH-2);
    pRecip=(const xtcomplexfloat*)recip;
    pZrd=(const xb_vecN_4xcf32*)(z);
    pZwr=(      xb_vecN_4xcf32*)(z);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t,x,y;
        pRecip=(const xtcomplexfloat*)(recip+l);
        for (m=0; m<N; m++)
        {
            pZk=(const xb_vecN_4xcf32*)(zk+SA*l);
            BBE_LSN_4XCF32_XP(t,pRecip,L*sizeof(complex_float));
            t=BBE_REPN_4XCF32(t,0);
            for (n=0; n<NN; n++) 
            {
                BBE_LVN_4XCF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_IP(y,pZk ,2*BBE_SIMD_WIDTH);
                BBE_MULSN_4XCF32(x,y,t);
                BBE_SVN_4XCF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
            }
        }
        pRecip-=N*L-1;
    }
}
#endif

// prepare elimination process: copy columns and clean them up
inline_ void cmatinvgjnxnnf_prepelimGS(complex_float * recip, complex_float * z, int k, int N, int L, int SA)
#if 0
{
    int m,l;

    for (l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            recip[m*L+l]= z[SA*l+m*N+k];
            if(m!=k) z[SA*l+m*N+k]=makecomplexf(0.f,0.f);
        }
    }
    for (l=0; l<L; l++)
    {
        recip[k*L+l]= makecomplexf(0.f,0.f);
    }
}
#else
{
    xb_vecN_2xc16 vm,vk;
    xb_vecN_4xcf32 x;
    const xtcomplexfloat* restrict  pZrd;
          xtcomplexfloat* restrict  pZwr;
          xtcomplexfloat* restrict  pRecip;
    xb_vecN_4xcf32 * restrict pWr;
    int m,l,nbytes,LL;
    valign awr;

    pZrd   =(const xtcomplexfloat*)(z+k);
    pZwr   =(      xtcomplexfloat*)(z+k);
    pRecip =(      xtcomplexfloat*)recip;

    vk=k; vm=0;
    for (m=0; m<N; m++)
    {
        vboolN_2 meqk;
        vboolN_4 meqk4;
        meqk=BBE_EQN_2XC16(vm,vk);
        {
            vboolN b,c0,c1;
            b=BBE_MOVN_FROMN_2(meqk);
            BBE_EXTRACTB(c1,c0,b);
            meqk4=BBE_MOVN_4_FROMN(c0);
        }
        vm=BBE_ADDN_2XC16(m,1);
        for (l=0; l<L; l++)
        {
            BBE_LSN_4XCF32_XP(x,pZrd,SA*sizeof(complex_float));
            BBE_SSN_4XCF32_IP(x,pRecip ,sizeof(complex_float));
            x=BBE_MOVN_4XCF32T(x,0.f,meqk4);
            BBE_SSN_4XCF32_XP(x,pZwr,SA*sizeof(complex_float));
        }
        pZrd+=-L*SA+N;
        pZwr+=-L*SA+N;
    }
    __Pragma("no_reorder")
    pWr=(xb_vecN_4xcf32 *)(recip+k*L);
    nbytes = (L * sizeof(complex_float)) & (2*BBE_SIMD_WIDTH-1);
    LL = L>>(LOG2_BBE_SIMD_WIDTH-2);
    awr=BBE_ZALIGN();
    for (l=0; l<LL; l++)
    {
        BBE_SAN_4XCF32_IP(0.f,awr,pWr);
    }
    BBE_SAVN_4XCF32_XP(0.f,awr,pWr,nbytes);
    BBE_SAN_4XCF32POS_FP(awr,pWr);
}
#endif

/* swap rows */
static void cmatinvgjnxnnf_swapRows(complex_float * z, const int16_t *srow, int k, int N, int L, int SA)
#if 0
{
    int l,n;
    for (l=0; l<L; l++)
    {
        int prow;
        prow=srow[k+N*l];
        for(n=0; n<N; n++) 
        {
            complex_float t;
            t=z[SA*l+prow*N+n];
            z[SA*l+prow*N+n]=z[SA*l+k*N+n];
            z[SA*l+k*N+n]=t;
        }
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pRd0;
    const xb_vecN_4xcf32 * restrict pRd1;
          xb_vecN_4xcf32 * restrict pWr0;
          xb_vecN_4xcf32 * restrict pWr1;
    int l,n,NN;
    NASSERT(N*N==SA);   // for SIMD_WIDTH<=16
    (void)SA;
    NN = N>>(LOG2_BBE_SIMD_WIDTH-2);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t0,t1;
        int prow;
        prow=srow[k+N*l];
        pRd0=(const xb_vecN_4xcf32 *)(z+SA*l+prow*N);
        pRd1=(const xb_vecN_4xcf32 *)(z+SA*l+k*N   );
        pWr0=(      xb_vecN_4xcf32 *)(z+SA*l+prow*N);
        pWr1=(      xb_vecN_4xcf32 *)(z+SA*l+k*N   );
        for(n=0; n<NN; n++) 
        {
            BBE_LVN_4XCF32_IP(t0,pRd0,2*BBE_SIMD_WIDTH);
            BBE_LVN_4XCF32_IP(t1,pRd1,2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(t1,pWr0,2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(t0,pWr1,2*BBE_SIMD_WIDTH);
        }
    }
}
#endif

static void cmatinvgjnxnnf_swapCols(complex_float * z, const int16_t *scol, int k, int N, int L, int SA)
#if 0
{
    complex_float t;
    int l,n;
    for (l=0; l<L; l++)
    {
        int pcol;
        pcol=scol[k+N*l];
        /* swap columns */
        for(n=0; n<N; n++) 
        {
            t=z[SA*l+n*N+pcol];
            z[SA*l+n*N+pcol]=z[SA*l+n*N+k];
            z[SA*l+n*N+k]=t;
        }
    }
}
#else
{
    const xtcomplexfloat * restrict pRd0;
    const xtcomplexfloat * restrict pRd1;
          xtcomplexfloat * restrict pWr0;
          xtcomplexfloat * restrict pWr1;
    int l,n;
    NASSERT(N*N==SA);   // for SIMD_WIDTH<=16
    (void)SA;
    pRd1=(const xtcomplexfloat *)(z+k   );
    pWr1=(      xtcomplexfloat *)(z+k   );
    scol+=k;
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t0,t1;
        int pcol;
        pcol=scol[0]; scol+=N;
        pRd0=(const xtcomplexfloat *)(z+SA*l+pcol);
        pWr0=(      xtcomplexfloat *)(z+SA*l+pcol);
        /* swap columns */
        __Pragma("loop_count min=4")
        for(n=0; n<N; n++) 
        {
            BBE_LSN_4XCF32_XP(t0,pRd0,N*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t1,pRd1,N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(t1,pWr0,N*sizeof(complex_float));
            BBE_SSN_4XCF32_XP(t0,pWr1,N*sizeof(complex_float));
        }
       // pRd1+=SA-N*N; /* not needed for bbe32 because N is a multiple of 4 and matrix size will always be a multiple of simd width */
       // pWr1+=SA-N*N;
    }
}
#endif

/* pivoting: search the absolute maximum and its position in the matrix 
    at all positions excluding previously used 
*/
inline_  void cmatinvgjnxnnf_searchPivot(int16_t *srow,int16_t *scol,const complex_float * z,int k,int N,int L,int SA)
#if 1
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    xb_vecN_2xc16 idx,maxidx,addidx,seqc;
    const xb_vecN_2xf32 * restrict pZ;
    vboolN_4 bmask;
    int l,m,n;
    seqc   = BBE_SELN_2XC16I(0,BBE_MOVN_2XC16_FROMNX16(BBE_SEQNX16()),BBE_SELI_INTERLEAVE_1_LO);
    addidx = BBE_SELN_2XC16I(0,BBE_SIMD_WIDTH/4,BBE_SELI_INTERLEAVE_2_LO);
    xtfloat maxT;
    vselN vmax;
    int dummy;
    xb_vecN_2xf32 MAXVAL,A,T;
    xb_vecNx16 pivot;
    vboolN_2 bmax;
    vboolN_2 cond;

    pZ=(const xb_vecN_2xf32*)(z);
    for (l=0; l<L; l++)
    {
        MAXVAL=FLT_MIN;
        pZ+=(k*N)>>(LOG2_BBE_SIMD_WIDTH-2);
        for (m=k; m<N; m++)
        {
            bmask=BBE_LTRN_4(k&3);
            idx=BBE_SELN_2XC16I(m,BBE_ADDN_2XC16(seqc,(k&~3)),BBE_SELI_INTERLEAVE_2_LO);
            pZ+=k>>(LOG2_BBE_SIMD_WIDTH-2);
            __Pragma("loop_count min=1")
            for(n=0; n<((N-k+3)>>2); n++) 
            {
                BBE_LVN_2XF32F_IP(A,pZ,2*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bmask)));
                T=BBE_MULMN_2XF32(A,A,0,0); 
                BBE_MULMASN_2XF32(T,A,A,0,15); 
                cond=BBE_OGTN_2XF32(T,MAXVAL);
                MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
                maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
                idx=BBE_ADDN_2XC16(idx,addidx);
                bmask=BBE_LTRN_4I(0);
            }
        }
        BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
        (void)maxT;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax)); (void)dummy;
        pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
        pivot=BBE_MINUNX16(pivot,N-1);
        BBE_SSNX16_XP(pivot,pCol,N*sizeof(int16_t));
        BBE_SSNX16_XP(BBE_REPNX16(pivot,2),pRow,N*sizeof(int16_t));
    }
}
#else
{
    int l,m,n;
    float32_t t,maxVal;
    for (l=0; l<L; l++)
    {
        int prow,pcol;
        prow=pcol=0; maxVal=FLT_MIN;
        for (m=k; m<N; m++)
        {
            for(n=k; n<N; n++) 
            {
                t=xqrc(z[SA*l+m*N+n]);  
                if(t>maxVal) { maxVal=t; prow=m; pcol=n; }
            }
        }
        srow[k+N*l]=prow;
        scol[k+N*l]=pcol;
    }
}
#endif
/*-------------------------------------------------------------------------
Inversion of Block Ordered Matrices By Gauss-Jordan Algortihm

Description: perform in-place inversion of real/complex matrices by Gauss-
Jordan elimination method, with full pivoting. The algorithm is applied to
a sequence of input matrices stored in block order.

Inversion result is not defined for a close to singular input matrix.

Storage size SA denote the number of data elements required to store a
matrix an NxN matrix A in block order. If matrix size is less than the
SIMD vector size, then the storage_size(matrix_size) equals the matrix_size
rounded up to the next power of two, otherwise it is matrix_size rounded up
to the next multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by 
            functions [c]matinvgj<size>nf_getScratchSize(N,L)
Input:
  N         Matrix size
  L         Number of matrices
Input/Output:
  A[L][SA]  Input matrices, inverted matrices on output
Restrictions:
  pScr,A    Must not overlap and must be aligned on 32-byte boundary 
  N         Must be a positive multiple of 4
---------------------------------------------------------------------------*/
void cmatinvgjnxnnf ( void* pScr, complex_float * restrict z, int N, int L  )
{
    int16_t *srow,*scol;
    complex_float * temp;   // [N*L]
    int k,SA;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(N%4==0 && N>0);

    if(L<=0 || N<=0) return;

    temp=(complex_float*)(pScr);
    srow=(int16_t*)(temp+N*L);
    scol=srow+N*L;
    NASSERT_ALIGN32(temp);

    SA=getSpace(N*N);
    for (k=0; k<N; k++)
    {
        cmatinvgjnxnnf_searchPivot(srow,scol,z,k,N,L,SA);
        cmatinvgjnxnnf_swapRows(z,srow,k,N,L,SA);
        cmatinvgjnxnnf_swapCols(z,scol,k,N,L, SA);
        cmatinvgjnxnnf_recipPivot(temp, (z+k*N)+k, L, SA);
        cmatinvgjnxnnf_updPivot( z+k*N, temp, N, L, SA);
        cmatinvgjnxnnf_prepelimGS(temp,z,k,N,L,SA);
        cmatinvgjnxnnf_elimGS(z,z+k*N, temp,N,L,SA);
    }
    /* final reverse permulation of columns/rows  */
    for (k=N-2; k>=0; k--)
    {
        cmatinvgjnxnnf_swapCols(z,srow,k,N,L, SA);
        cmatinvgjnxnnf_swapRows(z,scol,k,N,L,SA);
    }
}

size_t cmatinvgjnxnnf_getScratchSize (int N, int L  )
{
    NASSERT(N%4==0 && N>0);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int16_t)*N*L*2+sizeof(complex_float)*N*L;
}


#else
DISCARD_FUN(void, cmatinvgjnxnnf,  (void* pScr, complex_float * restrict z, int N, int L  ))

size_t cmatinvgjnxnnf_getScratchSize (int N, int L  )
{
    NASSERT(N%4==0 && N>0);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}

#endif

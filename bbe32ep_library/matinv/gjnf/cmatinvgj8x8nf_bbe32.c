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
    format, 8x8
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "cmatinvgjnxnnf_common.h"
#include "common.h"
#include <math.h>
#include <float.h>

#undef BBE_SELN_4XCF32
#define BBE_SELN_4XCF32(b,c,d) BBE_MOVN_4XCF32_FROMNX16((BBE_SELNX16((BBE_MOVNX16_FROMN_4XCF32(b)),(BBE_MOVNX16_FROMN_4XCF32(c)),d)))



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

#define METHOD_ELIM 1 // elimination method (0 - clean up columns except for pivot and make faster update, 1 - update via circular addressing)

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


inline_ void cmatinvgj8x8nf_updPivot(complex_float *z, complex_float* recip,  int L)
#if 0
{
    const int N=8,SA=64;
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
    int l;

    NASSERT_ALIGN32(z);
    pRecip=(const xtcomplexfloat*)recip;
    pZrd=(const xb_vecN_4xcf32*)(z);
    pZwr=(      xb_vecN_4xcf32*)(z);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t,x;
        BBE_LSN_4XCF32_IP(t,pRecip,sizeof(complex_float));
        t=BBE_REPN_4XCF32(t,0);
        BBE_LVN_4XCF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
        x=BBE_MULN_4XCF32(x,t);
        BBE_SVN_4XCF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(x,pZrd,15*2*BBE_SIMD_WIDTH);
        x=BBE_MULN_4XCF32(x,t);
        BBE_SVN_4XCF32_XP(x,pZwr,15*2*BBE_SIMD_WIDTH);
    }
}
#endif

#if METHOD_ELIM==0
/* GS elimination : for all rows excluding the pivot one */
inline_ void cmatinvgj8x8nf_elimGS(complex_float *z, complex_float *zk, const complex_float* recip, int L)
#if 0
{
    const int N=8, SA=64;
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
    int m,l;

    NASSERT_ALIGN32(z);
    NASSERT_ALIGN32(zk);
    pRecip=(const xtcomplexfloat*)recip;
    pZrd=(const xb_vecN_4xcf32*)(z);
    pZwr=(      xb_vecN_4xcf32*)(z);
    pZk=(const xb_vecN_4xcf32*)(zk);
    pRecip=(const xtcomplexfloat*)(recip);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t,x,y0,y1;
        BBE_LVN_4XCF32_IP(y0,pZk ,   2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(y1,pZk ,15*2*BBE_SIMD_WIDTH );
        __Pragma("loop_count min=8,max=8")
        for (m=0; m<8; m++)
        {
            BBE_LSN_4XCF32_XP(t,pRecip,L*sizeof(complex_float));
            t=BBE_REPN_4XCF32(t,0);
            BBE_LVN_4XCF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
            BBE_MULSN_4XCF32(x,y0,t);
            BBE_SVN_4XCF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
            BBE_LVN_4XCF32_IP(x,pZrd,2*BBE_SIMD_WIDTH);
            BBE_MULSN_4XCF32(x,y1,t);
            BBE_SVN_4XCF32_IP(x,pZwr,2*BBE_SIMD_WIDTH);
        }
        pRecip-=8*L-1;
    }
}
#endif

// prepare elimination process: copy columns and clean them up
inline_ void cmatinvgj8x8nf_prepelimGS(complex_float * recip, complex_float * z, int k, int L)
#if 0
{
    const int N=8, SA=64;
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
    const int N=8, SA=64;
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

#elif METHOD_ELIM==1
inline_ void cmatinvgj8x8nf_elimGS_alt(complex_float *z, int k, int L)
#if 0
{
    const int N=8,SA=64;
    int m,n,l;
    NASSERT_ALIGN32(z);
    for(l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            complex_float t;
            if(m==k) continue; 
            t=z[l*SA+m*N+k];
            z[l*SA+m*N+k]=makecomplexf(0.0f,0.0f);
            for (n=0; n<N; n++) 
            {
                z[l*SA+m*N+n]=subc(z[l*SA+m*N+n],mulc(z[l*SA+k*N+n],t)); 
            }
        }
    }
}
#else
{
    uintptr_t zbegin,zend;
    xb_vecNx16 seq0=BBE_SEQNX16();
    xb_vecNx16 seq1=BBE_ADDNX16(seq0,BBE_SIMD_WIDTH);
    vboolN_4 bk0,bk1;
    vselN selk;
    const xb_vecN_4xcf32 * restrict pZrd;
          xb_vecN_4xcf32 * restrict pZwr;
    int l;

    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    selk=BBE_MOVVSV(BBE_MOVNX16_FROMN_4X64(BBE_REPN_4X64(BBE_MOVN_4X64_FROMNX16(BBE_ADDNX16(BBE_SEQNX16(),(k<<2))),0)),0); /* select for replicaiton of k-th element */
    bk0 =BBE_MOVN_4_FROMN(BBE_EQNX16(seq0,BBE_SELNX16(seq1,seq0,selk)));  /* boolean 1 in k-th position */
    bk1 =BBE_MOVN_4_FROMN(BBE_EQNX16(seq1,BBE_SELNX16(seq1,seq0,selk)));  

    zbegin=(uintptr_t)z;
    zend  =zbegin+64*sizeof(complex_float);
    pZrd = (const xb_vecN_4xcf32 *)&z[k*8];
    WUR_CBEGIN(zbegin);
    WUR_CEND(zend  );
    for(l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t,t0,t1,zm0,zm1,zk0,zk1;
        BBE_LVN_4XCF32_IP(zk0,pZrd,2*BBE_SIMD_WIDTH); 
        BBE_LVN_4XCF32_IC(zk1,pZrd); 
        pZwr=(xb_vecN_4xcf32*)pZrd;
#if 0
        int n;
        __Pragma("loop_count min=7,max=7")
        for (n=0;n<7;n++)
        {
            BBE_LVN_4XCF32_IC(t0,pZrd); 
            BBE_LVN_4XCF32_IC(t1,pZrd); 
            zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
            zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
            t=BBE_SELN_4XCF32(t1,t0,selk);
            BBE_MULSN_4XCF32(zm0,zk0,t);
            BBE_MULSN_4XCF32(zm1,zk1,t);
            BBE_SVN_4XCF32_IC(zm0,pZwr); 
            BBE_SVN_4XCF32_IC(zm1,pZwr); 
        }
#else
        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 

        BBE_LVN_4XCF32_IC(t0,pZrd); 
        BBE_LVN_4XCF32_IC(t1,pZrd); 
        zm0=BBE_MOVN_4XCF32T(0.f,t0,bk0);
        zm1=BBE_MOVN_4XCF32T(0.f,t1,bk1);
        t=BBE_SELN_4XCF32(t1,t0,selk);
        BBE_MULSN_4XCF32(zm0,zk0,t);
        BBE_MULSN_4XCF32(zm1,zk1,t);
        BBE_SVN_4XCF32_IC(zm0,pZwr); 
        BBE_SVN_4XCF32_IC(zm1,pZwr); 
#endif
        zbegin+=64*sizeof(complex_float);
        zend  +=64*sizeof(complex_float);
        WUR_CBEGIN(zbegin);
        WUR_CEND(zend  );
        pZrd+=16;
    }
}
#endif
#endif

/* swap rows */
static void cmatinvgj8x8nf_swapRows(complex_float * z, const int16_t *srow, int k, int L)
#if 0
{
    const int N=8, SA=64;
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
    const xb_vecN_4xcf32 * restrict pRd;
          xb_vecN_4xcf32 * restrict pWr;
    int l;
    srow+=k;
    pRd=(const xb_vecN_4xcf32 *)(z+k*8);
    pWr=(      xb_vecN_4xcf32 *)(z+k*8);
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 t0,t1,t2,t3;
        int prow;
        prow=srow[0]-k; srow+=8;
        t1=BBE_LVN_4XCF32_X (pRd,prow*8*sizeof(complex_float));
        t3=BBE_LVN_4XCF32_X (pRd,prow*8*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        t2=BBE_LVN_4XCF32_I (pRd,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(t0,pRd,64*sizeof(complex_float));
        BBE_SVN_4XCF32_X (t0,pWr,prow*8*sizeof(complex_float));
        BBE_SVN_4XCF32_X (t2,pWr,prow*8*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_I (t3,pWr,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(t1,pWr,64*sizeof(complex_float));
    }
}
#endif

static void cmatinvgj8x8nf_swapCols(complex_float * z, const int16_t *scol, int k, int L)
#if 0
{
    const int N=8, SA=64;
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
    const xtcomplexfloat * restrict pRd;
          xtcomplexfloat * restrict pWr;
    int l;
    pRd=(const xtcomplexfloat *)(z+k   );
    pWr=(      xtcomplexfloat *)(z+k   );
    scol+=k;
    for (l=0; l<L; l++,z+=64)
    {
        xb_vecN_4xcf32 t0,t1;
        int pcol;
        pcol=scol[0]-k; scol+=8;
        /* swap columns */
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
        t1=BBE_LSN_4XCF32_X (pRd,pcol*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(t0,pRd,8*sizeof(complex_float));
        BBE_SSN_4XCF32_X (t0,pWr,pcol*sizeof(complex_float));
        BBE_SSN_4XCF32_XP(t1,pWr,8*sizeof(complex_float));
    }
}
#endif

/* pivoting: search the absolute maximum and its position in the matrix 
    at all positions excluding previously used 
*/
inline_  void cmatinvgj8x8nf_searchPivot(int16_t *srow,int16_t *scol,const complex_float * z,int k,int L)
#if 1
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    xb_vecN_2xc16 idx,maxidx,addidx,addidx1,seqc;
    const xb_vecN_2xf32 * restrict pZ;
    vboolN_4 bmask;
    int l,m;
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

    if (k<4)
    {
        /*--------------------------------------------------------*/
        /* k=0...3                                                */
        /*--------------------------------------------------------*/
        bmask=BBE_LTRN_4(k&3);
        addidx1=BBE_SELN_2XC16I(1,0,BBE_SELI_INTERLEAVE_2_LO);
        for (l=0; l<L; l++)
        {
            MAXVAL=FLT_MIN;
            pZ+=(k*8)>>(LOG2_BBE_SIMD_WIDTH-2);
            idx=BBE_SELN_2XC16I(k,BBE_ADDN_2XC16(seqc,0),BBE_SELI_INTERLEAVE_2_LO);
            __Pragma("loop_count min=5,max=8")
            for (m=k; m<8; m++)
            {
                BBE_LVN_2XF32F_IP(A,pZ,2*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bmask)));
                T=BBE_MULMN_2XF32(A,A,0,0); 
                BBE_MULMASN_2XF32(T,A,A,0,15); 
                cond=BBE_OGTN_2XF32(T,MAXVAL);
                MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
                maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);

                BBE_LVN_2XF32_IP(A,pZ,2*BBE_SIMD_WIDTH);
                T=BBE_MULMN_2XF32(A,A,0,0); 
                BBE_MULMASN_2XF32(T,A,A,0,15); 
                cond=BBE_OGTN_2XF32(T,MAXVAL);
                MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
                maxidx=BBE_MOVN_2XC16T(BBE_ADDN_2XC16(idx,addidx),maxidx,cond);
                idx=BBE_ADDN_2XC16(idx,addidx1);
            }
            BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
            (void)maxT;
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
            pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
            pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,8-1));
            BBE_SSNX16_IP(pivot,pCol,8*sizeof(int16_t));
            BBE_SSNX16_IP(BBE_REPNX16(pivot,2),pRow,8*sizeof(int16_t));
        }
    }
    else
    {
        /*--------------------------------------------------------*/
        /* k=4...7                                                */
        /*--------------------------------------------------------*/
        vboolN_4 maskk,bltk0,bltk1,bltk2,bltk3;
        maskk=~BBE_LTRN_4(k&3);
        bltk0 = ((k&3)<=0 ? ~BBE_LTRN_4(0): BBE_LTRN_4(0)) & maskk;
        bltk1 = ((k&3)<=1 ? ~BBE_LTRN_4(0): BBE_LTRN_4(0)) & maskk;
        bltk2 = ((k&3)<=2 ? ~BBE_LTRN_4(0): BBE_LTRN_4(0)) & maskk;
        bltk3 = ((k&3)<=3 ? ~BBE_LTRN_4(0): BBE_LTRN_4(0)) & maskk;

        seqc=BBE_ADDN_2XC16(seqc,4);
        addidx=BBE_SELN_2XC16I(1,0,BBE_SELI_INTERLEAVE_2_LO);
        pZ+=1;
        pZ+=8;
        for (l=0; l<L; l++)
        {
            MAXVAL=FLT_MIN;
            idx=BBE_SELN_2XC16I(4,seqc,BBE_SELI_INTERLEAVE_2_LO);
            BBE_LVN_2XF32T_IP(A,pZ,4*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bltk0)));
            T=BBE_MULMN_2XF32(A,A,0,0);  BBE_MULMASN_2XF32(T,A,A,0,15); 
            MAXVAL=T;
            maxidx=idx;

            idx=BBE_ADDN_2XC16(idx,addidx);
            BBE_LVN_2XF32T_IP(A,pZ,4*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bltk1)));
            T=BBE_MULMN_2XF32(A,A,0,0); BBE_MULMASN_2XF32(T,A,A,0,15); 
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx);
            BBE_LVN_2XF32T_IP(A,pZ,4*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bltk2)));
            T=BBE_MULMN_2XF32(A,A,0,0); BBE_MULMASN_2XF32(T,A,A,0,15); 
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx);
            BBE_LVN_2XF32T_IP(A,pZ,(8+2)*2*BBE_SIMD_WIDTH,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bltk3)));
            T=BBE_MULMN_2XF32(A,A,0,0); BBE_MULMASN_2XF32(T,A,A,0,15); 
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx);
            BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
            (void)maxT;
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
            pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
            pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,8-1));
            BBE_SSNX16_XP(pivot,pCol,8*sizeof(int16_t));
            BBE_SSNX16_XP(BBE_REPNX16(pivot,2),pRow,8*sizeof(int16_t));
        }
    }
}
#else
{
    const int N=8, SA=64;
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
void cmatinvgj8x8nf ( void* pScr, complex_float * restrict z, int L  )
{
    int16_t *srow,*scol;
    complex_float * temp;   // [N*L]
    int k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);

    if(L<=0 ) return;

    temp=(complex_float*)(pScr);
    srow=(int16_t*)(temp+8*L);
    scol=srow+8*L;
    NASSERT_ALIGN32(temp);

    for (k=0; k<8; k++)
    {
        cmatinvgj8x8nf_searchPivot(srow,scol,z,k,L);
        cmatinvgj8x8nf_swapRows(z,srow,k,L);
        cmatinvgj8x8nf_swapCols(z,scol,k,L);
        cmatinvgjnxnnf_recipPivot(temp, (z+k*8)+k, L, 8*8);
        cmatinvgj8x8nf_updPivot( z+k*8, temp, L);
#if METHOD_ELIM==0
        cmatinvgj8x8nf_prepelimGS(temp,z,k,L);
        cmatinvgj8x8nf_elimGS(z,z+k*8, temp,L);
#elif METHOD_ELIM==1
        cmatinvgj8x8nf_elimGS_alt(z,k,L);
#else
#error unsupported
#endif
    }
    /* final reverse permulation of columns/rows  */
    for (k=8-2; k>=0; k--)
    {
        cmatinvgj8x8nf_swapCols(z,srow,k,L);
        cmatinvgj8x8nf_swapRows(z,scol,k,L);
    }
}


size_t cmatinvgj8x8nf_getScratchSize (int N, int L  )
{
    NASSERT(N==8);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int16_t)*8*L*2+sizeof(complex_float)*8*L;
}


#else
DISCARD_FUN(void, cmatinvgj8x8nf,  (void* pScr, complex_float * restrict z, int L  ))

size_t cmatinvgj8x8nf_getScratchSize (int N, int L  )
{
    NASSERT(N==8);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}

#endif

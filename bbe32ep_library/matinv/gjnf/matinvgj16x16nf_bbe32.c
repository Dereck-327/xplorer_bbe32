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
    format, 16x16
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_Math.h"
#include "matinvgjnxnnf_common.h"
#include "common.h"
#include <math.h>
#include <float.h>

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

#undef BBE_SELN_2XF32
#define BBE_SELN_2XF32(b,c,d) BBE_MOVN_2XF32_FROMNX16((BBE_SELNX16((BBE_MOVNX16_FROMN_2XF32(b)),(BBE_MOVNX16_FROMN_2XF32(c)),d)))


inline_ void matinvgj16x16nf_updPivot(float32_t *z, float32_t* recip, int L)
#if 0
{
    const int SA=16*16,N=16;
    int n,l;
    NASSERT_ALIGN32(z);
    NASSERT(N==16 && SA==256 && L>0);
    for (l=0; l<L; l++)
    {
        for (n=0; n<N; n++) z[SA*l+n]*=recip[l]; 
    }
}
#else
{
    const xtfloat* restrict pRecip;
    const xb_vecN_2xf32* restrict pZrd;
          xb_vecN_2xf32* restrict pZwr;
    int l;

    NASSERT_ALIGN32(z);
    NASSERT( L>0);
    pRecip=(const xtfloat*)recip;
    pZrd=(const xb_vecN_2xf32*)(z);
    pZwr=(      xb_vecN_2xf32*)(z);
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,x0,x1;
        BBE_LSN_2XF32_IP(t,pRecip,sizeof(float32_t));
        t=BBE_REPN_2XF32(t,0);
        x1=BBE_LVN_2XF32_I(pZrd,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(x0,pZrd,256*sizeof(float32_t));
        x0=BBE_MULN_2XF32(x0,t);
        x1=BBE_MULN_2XF32(x1,t);
        BBE_SVN_2XF32_I (x1,pZwr,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_XP(x0,pZwr,256*sizeof(float32_t));
    }
}
#endif

/* GS elimination : for all rows excluding the pivot one */
inline_ void matinvgj16x16nf_elimGS(float32_t *z, int k, int L)
#if 0
{
    const int SA=16*16,N=16;
    int m,n,l;
    for(l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            float32_t t;
            if(m==k) continue; 
            t=z[l*SA+m*N+k];
            z[l*SA+m*N+k]=0.0f;
            for (n=0; n<N; n++) 
            {
                z[l*SA+m*N+n]-=z[l*SA+k*N+n]*t; 
            }
        }
    }
}
#elif 0
{
    const int N=16,SA=16*16;
    xb_vecNx16 seq0=BBE_SEQNX16();
    xb_vecNx16 seq1=BBE_ADDNX16(seq0,BBE_SIMD_WIDTH);
    vboolN_2 bk0,bk1;
    const xb_vecN_2xf32 *pZk ;
    const xb_vecN_2xf32 *pZrd;
          xb_vecN_2xf32 *pZwr;
    int m,n,l;
    vselN selk;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    (void)N,(void)SA;

    selk=BBE_MOVVSV(BBE_REPNX16C(BBE_ADDNX16(seq0,(k<<1)),0),0);          /* select for replicaiton of k-th element */
    bk0 =BBE_MOVN_2_FROMN(BBE_EQNX16(seq0,BBE_SELNX16(seq1,seq0,selk)));  /* boolean 1 in k-th position */
    bk1 =BBE_MOVN_2_FROMN(BBE_EQNX16(seq1,BBE_SELNX16(seq1,seq0,selk)));  
    for(l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,t0,t1,zm0,zm1,zk0,zk1;
        pZk = (const xb_vecN_2xf32 *)&z[l*SA+k*N];
        zk0=BBE_LVN_2XF32_I(pZk,0);
        zk1=BBE_LVN_2XF32_I(pZk,2*BBE_SIMD_WIDTH);
        for (m=0; m<16; m++)
        {
            if (m==k) continue;
            pZrd = (const xb_vecN_2xf32 *)&z[l*SA+m*N];
            pZwr = (      xb_vecN_2xf32 *)&z[l*SA+m*N];
            t0=BBE_LVN_2XF32_I(pZrd,0); 
            t1=BBE_LVN_2XF32_I(pZrd,2*BBE_SIMD_WIDTH);
            zm0=BBE_MOVN_2XF32T(0.f,t0,bk0);
            zm1=BBE_MOVN_2XF32T(0.f,t1,bk1);
            t=BBE_SELN_2XF32(t1,t0,selk);
            BBE_MULSN_2XF32(zm0,zk0,t);
            BBE_MULSN_2XF32(zm1,zk1,t);
            BBE_SVN_2XF32_I(zm0,pZwr,0);
            BBE_SVN_2XF32_I(zm1,pZwr,2*BBE_SIMD_WIDTH);
        }
    }
}
#else
{
    const int N=16,SA=16*16;
    xb_vecNx16 seq0=BBE_SEQNX16();
    xb_vecNx16 seq1=BBE_ADDNX16(seq0,BBE_SIMD_WIDTH);
    vboolN_2 bk0,bk1;
    const xb_vecN_2xf32 *restrict pZrd;
          xb_vecN_2xf32 *restrict pZwr;
    int m,l;
    vselN selk;
    uintptr_t zbegin,zend;

    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    (void)N,(void)SA;

    selk=BBE_MOVVSV(BBE_REPNX16C(BBE_ADDNX16(seq0,(k<<1)),0),0);          /* select for replicaiton of k-th element */
    bk0 =BBE_MOVN_2_FROMN(BBE_EQNX16(seq0,BBE_SELNX16(seq1,seq0,selk)));  /* boolean 1 in k-th position */
    bk1 =BBE_MOVN_2_FROMN(BBE_EQNX16(seq1,BBE_SELNX16(seq1,seq0,selk)));  
    zbegin=(uintptr_t)z;
    zend  =zbegin+256*sizeof(float32_t);
    pZrd = (const xb_vecN_2xf32 *)&z[k*16];
    WUR_CBEGIN(zbegin);
    WUR_CEND(zend  );
    for(l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,t0,t1,zm0,zm1,zk0,zk1;
        BBE_LVN_2XF32_IC(zk0,pZrd);
        BBE_LVN_2XF32_IC(zk1,pZrd);
        pZwr=(xb_vecN_2xf32*)pZrd;
        __Pragma("loop_count min=15,max=15")
        for (m=0; m<15; m++)
        {
            BBE_LVN_2XF32_IC(t0,pZrd);
            BBE_LVN_2XF32_IC(t1,pZrd);
            zm0=BBE_MOVN_2XF32T(0.f,t0,bk0);
            zm1=BBE_MOVN_2XF32T(0.f,t1,bk1);
            t=BBE_SELN_2XF32(t1,t0,selk);
            BBE_MULSN_2XF32(zm0,zk0,t);
            BBE_MULSN_2XF32(zm1,zk1,t);
            BBE_SVN_2XF32_IC(zm0,pZwr);
            BBE_SVN_2XF32_IC(zm1,pZwr);
        }
        zbegin+=256*sizeof(float32_t);
        zend  +=256*sizeof(float32_t);
        WUR_CBEGIN(zbegin);
        WUR_CEND(zend  );
        pZrd+=32;
    }
}
#endif

/* swap rows */
static void matinvgj16x16nf_swapRows(float32_t * z, const int16_t *srow, int k, int L)
#if 0
{
    const int SA=16*16,N=16;
    int l,n;
    NASSERT_ALIGN32(z);
    NASSERT(N==16 && SA==256 && L>0);
    for (l=0; l<L; l++)
    {
        int prow;
        prow=srow[k+N*l];
        for(n=0; n<N; n++) 
        {
            float32_t t;
            t=z[SA*l+prow*N+n];
            z[SA*l+prow*N+n]=z[SA*l+k*N+n];
            z[SA*l+k*N+n]=t;
        }
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pRd;
          xb_vecN_2xf32 * restrict pWr;
    xb_vecN_2xf32 t0,t1,t2,t3;
    int l;
    NASSERT_ALIGN32(z);
    NASSERT(L>0);
    pRd=(const xb_vecN_2xf32 *)(z+k*16);
    pWr=(      xb_vecN_2xf32 *)(z+k*16);
    z+=k*16;
    for (l=0; l<L; l++)
    {
        int prow;
        prow=srow[k+16*l]-k;
        t1=BBE_LVN_2XF32_X (pRd,prow*16*sizeof(float32_t));
        t3=BBE_LVN_2XF32_X (pRd,prow*16*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        t2=BBE_LVN_2XF32_I (pRd,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(t0,pRd,256*sizeof(float32_t));
        BBE_SVN_2XF32_X (t0,pWr,prow*16*sizeof(float32_t));
        BBE_SVN_2XF32_X (t2,pWr,prow*16*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_I (t3,pWr,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_XP(t1,pWr,256*sizeof(float32_t));
    }
}
#endif

static void matinvgj16x16nf_swapCols(float32_t * z, const int16_t *scol, int k, int L)
#if 0
{
    const int SA=16*16,N=16;
    int l,n;
    for (l=0; l<L; l++)
    {
        int pcol;
        pcol=scol[k+N*l];
        /* swap columns */
        for(n=0; n<N; n++) 
        {
            float32_t t;
            t=z[SA*l+n*N+pcol];
            z[SA*l+n*N+pcol]=z[SA*l+n*N+k];
            z[SA*l+n*N+k]=t;
        }
    }
}
#else
{
    const xtfloat * restrict pRd0;
    const xtfloat * restrict pRd1;
          xtfloat * restrict pWr0;
          xtfloat * restrict pWr1;
    int l,n;
    pRd1=(const xtfloat *)(z+k   );
    pWr1=(      xtfloat *)(z+k   );
    scol+=k;
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 t0,t1;
        int pcol;
        pcol=scol[0]; scol+=16;
        pRd0=(const xtfloat *)(z+256*l+pcol);
        pWr0=(      xtfloat *)(z+256*l+pcol);
        /* swap columns */
        __Pragma("loop_count min=16, max=16")
        for(n=0; n<16; n++) 
        {
            BBE_LSN_2XF32_XP(t0,pRd0,16*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t1,pRd1,16*sizeof(float32_t));
            BBE_SSN_2XF32_XP(t1,pWr0,16*sizeof(float32_t));
            BBE_SSN_2XF32_XP(t0,pWr1,16*sizeof(float32_t));
        }
    }
}
#endif


/* pivoting: search the absolute maximum and its position in the matrix 
    at all positions excluding previously used 
*/
inline_  void matinvgj16x16nf_searchPivot(int16_t *srow,int16_t *scol,const float32_t * z,int k,int L)
#if 0
{
    const int SA=16*16,N=16;
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
                t=fabsf(z[SA*l+m*N+n]);  
                if(t>maxVal) { maxVal=t; prow=m; pcol=n; }
            }
        }
        srow[k+N*l]=prow;
        scol[k+N*l]=pcol;
    }
}
#else
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    xb_vecN_2xc16 idx,maxidx,addidx0,addidx1;
    const xb_vecN_2xf32 * restrict pZ;
    int l,m;
    xtfloat maxT;
    vselN vmax;
    int dummy;
    xb_vecN_2xf32 MAXVAL,T;
    xb_vecNx16 pivot;
    vboolN_2 bmax;
    vboolN_2 cond,maskk;
    vboolN_2 bltk[8];
    vboolN_2 * restrict pbltk;
    maskk = ~BBE_LTRN_2(k);

    /*--------------------------------------------------------*/
    /* big k(>=8)                                             */
    /*--------------------------------------------------------*/
    if (k>=8)
    {
        /*
            for this case, inner loop is unrolled. All unused data are masked with precomputed boolean 
            registers
        */
        addidx0 = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,BBE_SIMD_WIDTH/2,BBE_SELI_INTERLEAVE_1_LO));
        maskk = ~BBE_LTRN_2(k);
        pbltk=bltk;
        for (m=0; m<(k-8); m++) BBE_SBN_2_IP(maskk & BBE_LTRN_2(0),pbltk,sizeof(vboolN_2));
        for (; m<8; m++)        BBE_SBN_2_IP(maskk &~BBE_LTRN_2(0),pbltk,sizeof(vboolN_2));
        addidx1=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(1,0,BBE_SELI_INTERLEAVE_1_LO));

        __Pragma("no_reorder")
        pZ=(const xb_vecN_2xf32*)(&z[8*16+8]);
        for (l=0; l<L; l++)
        {
            vboolN_2 bk;
            idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(8,BBE_ADDNX16(BBE_SEQNX16(),8),BBE_SELI_INTERLEAVE_1_LO));
            pbltk=bltk;
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            MAXVAL=BBE_ABSN_2XF32(T);
            maxidx=idx;
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_IP(T,pZ,2*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
            BBE_LBN_2_IP(bk,pbltk,sizeof(vboolN_2));
            BBE_LVN_2XF32T_XP(T,pZ,2*2*BBE_SIMD_WIDTH+16*2*BBE_SIMD_WIDTH,bk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
            (void)maxT;
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
            (void)dummy;
            pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
            pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,16-1));
            BBE_SSNX16_XP(pivot,pCol,16*sizeof(int16_t));
            BBE_SSNX16_XP(BBE_REPNX16(pivot,1),pRow,16*sizeof(int16_t));
        }
        return;
    }

    /*--------------------------------------------------------*/
    /* small k (<8)                                           */
    /*--------------------------------------------------------*/
    addidx0 = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,BBE_SIMD_WIDTH/2,BBE_SELI_INTERLEAVE_1_LO));
    addidx1 = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(1,0,BBE_SELI_INTERLEAVE_1_LO));
    for (l=0; l<L; l++,z+=256)
    {
        MAXVAL=FLT_MIN;
        pZ=(const xb_vecN_2xf32*)(&z[k*16]);
        idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(k,BBE_SEQNX16(),BBE_SELI_INTERLEAVE_1_LO));
        for (m=k; m<16; m++)
        {
            BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,maskk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);

            BBE_LVN_2XF32_IP(T,pZ,2*BBE_SIMD_WIDTH);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(BBE_ADDN_2XC16(idx,addidx0),maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx1);
        }
        BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
        (void)maxT;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
        (void)dummy;
        pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
        pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,16-1));
        BBE_SSNX16_XP(pivot,pCol,16*sizeof(int16_t));
        BBE_SSNX16_XP(BBE_REPNX16(pivot,1),pRow,16*sizeof(int16_t));
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
void matinvgj16x16nf ( void* pScr, float32_t * restrict z, int L  )
{
    int16_t *srow,*scol;
    float32_t * temp;   // [N*L]
    int k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);

    if(L<=0 ) return;
    temp=(float32_t*)(pScr);
    srow=(int16_t*)(temp+16*L);
    scol=srow+16*L;
    NASSERT_ALIGN32(temp);

    for (k=0; k<16; k++)
    {
        matinvgj16x16nf_searchPivot(srow,scol,z,k,L);
        matinvgj16x16nf_swapRows(z,srow,k,L);
        matinvgj16x16nf_swapCols(z,scol,k,L);
        matinvgjnxnnf_recipPivot(temp, (z+k*16)+k, L,16*16);
        matinvgj16x16nf_updPivot( z+k*16, temp, L);
        matinvgj16x16nf_elimGS(z,k,L);
    }
    /* final reverse permulation of columns/rows  */
    for (k=16-2; k>=0; k--)
    {
        matinvgj16x16nf_swapCols(z,srow,k,L);
        matinvgj16x16nf_swapRows(z,scol,k,L);
    }
}

size_t matinvgj16x16nf_getScratchSize (int N, int L  )
{
    NASSERT(N==16);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int16_t)*N*L*2+sizeof(float32_t)*N*L;
}
#else
DISCARD_FUN(void, matinvgj16x16nf,(void* pScr, float32_t * restrict z, int L  ))

size_t matinvgj16x16nf_getScratchSize (int N, int L  )
{
    NASSERT(N==16);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}

#endif

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
    format, 8x8
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

#undef BBE_SELN_2XF32
#define BBE_SELN_2XF32(b,c,d) BBE_MOVN_2XF32_FROMNX16((BBE_SELNX16((BBE_MOVNX16_FROMN_2XF32(b)),(BBE_MOVNX16_FROMN_2XF32(c)),d)))


inline_ void matinvgj8x8nf_updPivot(float32_t *z, float32_t* recip, int L)
#if 0
{
    int n,l;
    NASSERT_ALIGN32(z);
    NASSERT(N==8 && SA==64 && L>0);
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
        xb_vecN_2xf32 t,x;
        BBE_LSN_2XF32_IP(t,pRecip,sizeof(float32_t));
        t=BBE_REPN_2XF32(t,0);
        BBE_LVN_2XF32_XP(x,pZrd,64*sizeof(float32_t));
        x=BBE_MULN_2XF32(x,t);
        BBE_SVN_2XF32_XP(x,pZwr,64*sizeof(float32_t));
    }
}
#endif

/* GS elimination : for all rows excluding the pivot one */
inline_ void matinvgj8x8nf_elimGS(float32_t *z, int k, int L)
#if 0
{
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
#elif 1
{
    uintptr_t zbegin,zend;
    vboolN_2 bk;
    vselN selk;
    const xb_vecN_2xf32 * restrict pZrd;
          xb_vecN_2xf32 * restrict pZwr;
    int l;

    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    bk = (k==BBE_SIMD_WIDTH/2-1) ?  ~BBE_LTRN_2(BBE_SIMD_WIDTH/2-1) : BBE_LTRN_2(k+1) &~ BBE_LTRN_2(k);
    selk=BBE_MOVVSV(BBE_REPNX16C(BBE_ADDNX16(BBE_SEQNX16(),(k<<1)),0),0); /* select for replicaiton of k-th element */

    zbegin=(uintptr_t)z;
    zend  =zbegin+64*sizeof(float32_t);
    pZrd = (const xb_vecN_2xf32 *)&z[k*8];
    WUR_CBEGIN(zbegin);
    WUR_CEND(zend  );
    for(l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,zm,zk;
        BBE_LVN_2XF32_IC(zk,pZrd); 
        pZwr=(xb_vecN_2xf32*)pZrd;
        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        BBE_LVN_2XF32_IC(t,pZrd); 
        zm=BBE_MOVN_2XF32T(0.f,t,bk);
        t=BBE_SELN_2XF32(t,t,selk);
        BBE_MULSN_2XF32(zm,zk,t);
        BBE_SVN_2XF32_IC(zm,pZwr); 

        zbegin+=64*sizeof(float32_t);
        zend  +=64*sizeof(float32_t);
        WUR_CBEGIN(zbegin);
        WUR_CEND(zend  );
        pZrd+=8;
    }
}
#else
{
    vboolN_2 bk;
    const xb_vecN_2xf32 *pZk ;
    const xb_vecN_2xf32 *pZrd;
          xb_vecN_2xf32 *pZwr;
    int m,n,l;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(N==8 && SA==64 && L>0);
    (void)N,(void)SA;

    bk = (k==BBE_SIMD_WIDTH/2-1) ?  ~BBE_LTRN_2(BBE_SIMD_WIDTH/2-1) : BBE_LTRN_2(k+1) &~ BBE_LTRN_2(k);
    for(l=0; l<L; l++)
    {
        pZk = (const xb_vecN_2xf32 *)&z[l*SA+k*N];
        for (m=0; m<k; m++)
        {
            xb_vecN_2xf32 t,zm,zk;
            pZrd = (const xb_vecN_2xf32 *)&z[l*SA+m*N];
            pZwr = (      xb_vecN_2xf32 *)&z[l*SA+m*N];
            t=BBE_LVN_2XF32_I(pZrd,0); 
            zk=BBE_LVN_2XF32_I(pZk,0);
            zm=BBE_MOVN_2XF32T(0.f,t,bk);
            t=BBE_REPN_2XF32(t,k);
            BBE_MULSN_2XF32(zm,zk,t);
            BBE_SVN_2XF32_I(zm,pZwr,0);
        }

        for (m=k+1; m<N; m++)
        {
            xb_vecN_2xf32 t,zm,zk;
            pZrd = (const xb_vecN_2xf32 *)&z[l*SA+m*N];
            pZwr = (      xb_vecN_2xf32 *)&z[l*SA+m*N];
            t=BBE_LVN_2XF32_I(pZrd,0); 
            zk=BBE_LVN_2XF32_I(pZk,0);
            zm=BBE_MOVN_2XF32T(0.f,t,bk);
            t=BBE_REPN_2XF32(t,k);
            BBE_MULSN_2XF32(zm,zk,t);
            BBE_SVN_2XF32_I(zm,pZwr,0);
        }
    }
}
#endif

/* swap rows */
static void matinvgj8x8nf_swapRows(float32_t * z, const int16_t *srow, int k, int L)
#if 0
{
    int l,n;
    NASSERT_ALIGN32(z);
    NASSERT(N==8 && SA==64 && L>0);
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
    xb_vecN_2xf32 t0,t1;
    int l;
    NASSERT_ALIGN32(z);
    NASSERT(L>0);
    pRd=(const xb_vecN_2xf32 *)(z+k*8);
    pWr=(      xb_vecN_2xf32 *)(z+k*8);
    z+=k*8;
    for (l=0; l<L; l++)
    {
        int prow;
        prow=srow[k+8*l]-k;
        t1=BBE_LVN_2XF32_X (pRd,prow*8*sizeof(float32_t));
        BBE_LVN_2XF32_XP(t0,pRd,64*sizeof(float32_t));
        BBE_SVN_2XF32_X (t0,pWr,prow*8*sizeof(float32_t));
        BBE_SVN_2XF32_XP(t1,pWr,64*sizeof(float32_t));
    }
}
#endif

static void matinvgj8x8nf_swapCols(float32_t * z, const int16_t *scol, int k, int L)
#if 0
{
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
    int j,l;
    const xb_vecNx16 * restrict pRd;
          xb_vecNx16 * restrict pWr;
    const short * restrict pCol;
    xb_vecNx16 perm,x;
    vselN vperm;
    vboolN bk,bj;
    xb_vecNx16 vj,vk;

    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    pCol=(const short *)(scol+k);
    pRd=(const xb_vecNx16*)z;
    pWr=(      xb_vecNx16*)z;
    bk=BBE_LTRN(k+1) & ~BBE_LTRN(k);
    vk=k;
    for (l=0; l<L; l++)
    {
        BBE_LSNX16_XP(vj,pCol,8*sizeof(int16_t));
        j=BBE_EXTRNX16C(vj,0);
        vj=BBE_REPNX16(vj,0);
        bj=BBE_LTRN(j+1) & ~BBE_LTRN(j);
        perm=BBE_MOVNX16T(k,BBE_SEQNX16(),bj);
        perm=BBE_MOVNX16T(vj,perm        ,bk);
        perm=BBE_SLLINX16(perm,1);
        perm=BBE_SELNX16I(BBE_ADDNX16(perm,1),perm,BBE_SELI_INTERLEAVE_1_LO);
        vperm=BBE_MOVVSV(perm,0);

        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x,pRd,2*BBE_SIMD_WIDTH); x=BBE_SHFLNX16(x,vperm); BBE_SVNX16_IP(x,pWr,2*BBE_SIMD_WIDTH);
    }
}
#endif

/* pivoting: search the absolute maximum and its position in the matrix 
    at all positions excluding previously used 
*/
inline_  void matinvgj8x8nf_searchPivot(int16_t *srow,int16_t *scol,const float32_t * z,int k,int L)
#if 0
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
                t=fabsf(z[SA*l+m*N+n]);  
                if(t>maxVal) { maxVal=t; prow=m; pcol=n; }
            }
        }
        srow[k+N*l]=prow;
        scol[k+N*l]=pcol;
    }
}
#elif 0
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    const xb_vecN_2xf32* restrict pZ;
    vboolN_2 maskk; // mask for masking first k elements
    xb_vecN_2xc16 idx,maxidx,addidx;
    xtfloat maxT;
    vselN vmax;
    int dummy;
    xb_vecN_2xf32 MAXVAL,T;
    xb_vecNx16 pivot;
    vboolN_2 bmax;
    vboolN_2 cond;
    int l,m;

    maskk = BBE_LTRN_2(k);
    addidx = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(1,0,BBE_SELI_INTERLEAVE_1_LO));

    pZ=(const xb_vecN_2xf32*)(z);
    for (l=0; l<L; l++)
    {
        pZ+=k;
        MAXVAL=FLT_MIN;
        maxidx=0;
        idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(k, BBE_SEQNX16(),BBE_SELI_INTERLEAVE_1_LO));
        for (m=k; m<8; m++)
        {
            BBE_LVN_2XF32F_IP(T,pZ,2*BBE_SIMD_WIDTH,maskk);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
            idx=BBE_ADDN_2XC16(idx,addidx);
        }
        BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
        (void)maxT;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
        pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
        BBE_SSNX16_XP(pivot,pCol,N*sizeof(int16_t));
        BBE_SSNX16_XP(BBE_REPNX16(pivot,1),pRow,N*sizeof(int16_t));
    }
}
#else
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    const xb_vecN_2xf32* restrict pZ;
    vboolN_2 maskk; // mask for masking first k elements
    xb_vecN_2xc16 idx,maxidx,addidx;
    xtfloat maxT;
    vselN vmax;
    int dummy;
    xb_vecN_2xf32 MAXVAL,T;
    xb_vecNx16 pivot;
    vboolN_2 bm,bmax;
    vboolN_2 cond;
    int l,m;
    vboolN_2 bltk[8];
    vboolN_2 * restrict pbltk;

    maskk = ~BBE_LTRN_2(k);
    addidx = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(1,0,BBE_SELI_INTERLEAVE_1_LO));

    pbltk=bltk;
    for (m=0; m<k; m++) BBE_SBN_2_IP(maskk & BBE_LTRN_2(0),pbltk,sizeof(vboolN_2));
    for (m=k; m<8; m++) BBE_SBN_2_IP(maskk &~BBE_LTRN_2(0),pbltk,sizeof(vboolN_2));

    __Pragma("no_reorder")
    pZ=(const xb_vecN_2xf32*)(z);
    pbltk=bltk;
    for (l=0; l<L; l++)
    {
        idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0, BBE_SEQNX16(),BBE_SELI_INTERLEAVE_1_LO));
        pbltk=bltk;

        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        MAXVAL=T;
        maxidx=idx;
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);

        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);
        BBE_LBN_2_IP(bm,pbltk,sizeof(vboolN_2));
        BBE_LVN_2XF32T_IP(T,pZ,2*BBE_SIMD_WIDTH,bm);
        T=BBE_ABSN_2XF32(T);
        cond=BBE_OGTN_2XF32(T,MAXVAL);
        MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
        maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        idx=BBE_ADDN_2XC16(idx,addidx);

        BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
        (void)maxT;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
        (void)dummy;
        pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
        pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,8-1));
        BBE_SSNX16_XP(pivot,pCol,8*sizeof(int16_t));
        BBE_SSNX16_XP(BBE_REPNX16(pivot,1),pRow,8*sizeof(int16_t));
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
void matinvgj8x8nf ( void* pScr, float32_t * restrict z, int L  )
{
    int16_t *srow,*scol;
    float32_t * temp;   // [N*L]
    int k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);

    if(L<=0 ) return;

    temp=(float32_t*)(pScr);
    srow=(int16_t*)(temp+8*L);
    scol=srow+8*L;
    NASSERT_ALIGN32(temp);

    for (k=0; k<8; k++)
    {
        matinvgj8x8nf_searchPivot(srow,scol,z,k,L);
        matinvgj8x8nf_swapRows(z,srow,k,L);
        matinvgj8x8nf_swapCols(z,scol,k,L);
        matinvgjnxnnf_recipPivot(temp, (z+k*8)+k, L,8*8);
        matinvgj8x8nf_updPivot( z+k*8, temp, L);
        matinvgj8x8nf_elimGS(z,k,L);
    }
    /* final reverse permulation of columns/rows  */
    for (k=8-2; k>=0; k--)
    {
        matinvgj8x8nf_swapCols(z,srow,k,L);
        matinvgj8x8nf_swapRows(z,scol,k,L);
    }
}

size_t matinvgj8x8nf_getScratchSize (int N, int L  )
{
    NASSERT(N==8);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int16_t)*N*L*2+sizeof(float32_t)*N*L;
}
#else
DISCARD_FUN(void, matinvgj8x8nf,  (void* pScr, float32_t * restrict z, int L  ))

size_t matinvgj8x8nf_getScratchSize (int N, int L  )
{
    NASSERT(N==8);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}

#endif

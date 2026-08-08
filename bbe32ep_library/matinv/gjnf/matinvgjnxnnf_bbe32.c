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
#include "matinvgjnxnnf_common.h"
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
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

inline_ void matinvgjnxnnf_updPivot(float32_t *z, float32_t* recip, int N, int L, int SA)
#if 0
{
    int l,n;
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
    valign ard,awr;
    int n,l,nbytes,NN;

    NASSERT(N%4==0 && N>0);
    awr=BBE_ZALIGN();
    nbytes = (N * sizeof(float32_t)) & (2*BBE_SIMD_WIDTH-1);
    NN = N>>(LOG2_BBE_SIMD_WIDTH-1);
    pRecip=(const xtfloat*)recip;
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,x;
        BBE_LSN_2XF32_IP(t,pRecip,sizeof(float32_t));
        t=BBE_REPN_2XF32(t,0);
        pZrd=(const xb_vecN_2xf32*)(z+SA*l);
        pZwr=(      xb_vecN_2xf32*)(z+SA*l);
        ard=BBE_LAN_2XF32_PP(pZrd);
        for (n=0; n<NN; n++) 
        {
            BBE_LAN_2XF32_IP(x,ard,pZrd);
            x=BBE_MULN_2XF32(x,t);
            BBE_SAN_2XF32_IP(x,awr,pZwr);
        }
        BBE_LAVN_2XF32_XP(x,ard,pZrd,nbytes);
        x=BBE_MULN_2XF32(x,t);
        BBE_SAVN_2XF32_XP(x,awr,pZwr,nbytes);
        BBE_SAN_2XF32POS_FP(awr,pZwr);
    }
}
#endif

/* GS elimination : for all rows excluding the pivot one */
inline_ void matinvgjnxnnf_elimGS(float32_t *z, float32_t *zk, float32_t* recip, int N, int L, int SA)
#if 0
{
    int l,m,n;
    for (l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            float32_t t;
            t=recip[m*L+l];
            for (n=0; n<N; n++) z[SA*l+m*N+n]-=zk[SA*l+n]*t; 
        }
    }
}
#else
{
    const xtfloat* restrict pRecip;
    const xb_vecN_2xf32* restrict pZrd;
    const xb_vecN_2xf32* restrict pZk;
          xb_vecN_2xf32* restrict pZwr;
    valign ak,ard,awr;
    int m,n,l,nbytes,NN;

    NASSERT(N%4==0 && N>0);
    awr=BBE_ZALIGN();
    nbytes = (N * sizeof(float32_t)) & (2*BBE_SIMD_WIDTH-1);
    NN = N>>(LOG2_BBE_SIMD_WIDTH-1);
    pRecip=(const xtfloat*)recip;
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 t,x,y;
        pRecip=(const xtfloat*)(recip+l);
        pZrd=(const xb_vecN_2xf32*)(z+SA*l);
        pZwr=(      xb_vecN_2xf32*)(z+SA*l);
        ard=BBE_LAN_2XF32_PP(pZrd);
        for (m=0; m<N; m++)
        {
            pZk=(const xb_vecN_2xf32*)(zk+SA*l);
            ak=BBE_LAN_2XF32_PP(pZk);
            BBE_LSN_2XF32_XP(t,pRecip,L*sizeof(float32_t));
            t=BBE_REPN_2XF32(t,0);
            for (n=0; n<NN; n++) 
            {
                BBE_LAN_2XF32_IP(x,ard,pZrd);
                BBE_LAN_2XF32_IP(y,ak ,pZk );
                BBE_MULSN_2XF32(x,y,t);
                BBE_SAN_2XF32_IP(x,awr,pZwr);
            }
            BBE_LAVN_2XF32_XP(x,ard,pZrd,nbytes);
            BBE_LAVN_2XF32_XP(y,ak ,pZk ,nbytes);
            BBE_MULSN_2XF32(x,y,t);
            BBE_SAVN_2XF32_XP(x,awr,pZwr,nbytes);
        }
        BBE_SAN_2XF32POS_FP(awr,pZwr);
        pRecip-=N*L-1;
    }
}
#endif

// prepare elimination process: copy columns and clean them up
inline_ void matinvgjnxnnf_prepelimGS(float32_t * recip, float32_t * z, int k, int N, int L, int SA)
#if 0
{
    int m,l;

    for (l=0; l<L; l++)
    {
        for (m=0; m<N; m++)
        {
            recip[m*L+l]= z[SA*l+m*N+k];
            if(m!=k) z[SA*l+m*N+k]=0.0f;
        }
    }
    for (l=0; l<L; l++)
    {
        recip[k*L+l]= 0.f;
    }
}
#else
{
    xb_vecN_2xc16 vm,vk;
    xb_vecN_2xf32 x;
    const xtfloat* restrict  pZrd;
          xtfloat* restrict  pZwr;
          xtfloat* restrict  pRecip;
    xb_vecN_2xf32 * restrict pWr;
    int m,l,nbytes,LL;
    valign awr;

    pZrd   =(const xtfloat*)(z+k);
    pZwr   =(      xtfloat*)(z+k);
    pRecip =(      xtfloat*)recip;

    vk=k; vm=0;
    for (m=0; m<N; m++)
    {
        vboolN_2 meqk;
        meqk=BBE_EQN_2XC16(vm,vk);
        vm=BBE_ADDN_2XC16(m,1);
        for (l=0; l<L; l++)
        {
            BBE_LSN_2XF32_XP(x,pZrd,SA*sizeof(float32_t));
            BBE_SSN_2XF32_IP(x,pRecip,sizeof(float32_t));
            x=BBE_MOVN_2XF32T(x,0.f,meqk);
            BBE_SSN_2XF32_XP(x,pZwr,SA*sizeof(float32_t));
        }
        pZrd+=-L*SA+N;
        pZwr+=-L*SA+N;
    }
    __Pragma("no_reorder")
    pWr=(xb_vecN_2xf32 *)(recip+k*L);
    nbytes = (L * sizeof(float32_t)) & (2*BBE_SIMD_WIDTH-1);
    LL = L>>(LOG2_BBE_SIMD_WIDTH-1);
    awr=BBE_ZALIGN();
    for (l=0; l<LL; l++)
    {
        BBE_SAN_2XF32_IP(0.f,awr,pWr);
    }
    BBE_SAVN_2XF32_XP(0.f,awr,pWr,nbytes);
    BBE_SAN_2XF32POS_FP(awr,pWr);
}
#endif

/* swap rows */
static void matinvgjnxnnf_swapRows(float32_t * z, const int16_t *srow, int k, int N, int L, int SA)
#if 0
{
    int l,n;
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
    const xb_vecN_2xf32 * restrict pRd0;
    const xb_vecN_2xf32 * restrict pRd1;
          xb_vecN_2xf32 * restrict pWr0;
          xb_vecN_2xf32 * restrict pWr1;
    valign aRd0,aRd1,aWr0,aWr1;
    int l,n,nbytes,NN;
    nbytes = (N * sizeof(float32_t)) & (2*BBE_SIMD_WIDTH-1);
    NN = N>>(LOG2_BBE_SIMD_WIDTH-1);
    aWr0=BBE_ZALIGN();
    aWr1=BBE_ZALIGN();
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 t0,t1;
        int prow;
        prow=srow[k+N*l];
        pRd0=(const xb_vecN_2xf32 *)(z+SA*l+prow*N);
        pRd1=(const xb_vecN_2xf32 *)(z+SA*l+k*N   );
        pWr0=(      xb_vecN_2xf32 *)(z+SA*l+prow*N);
        pWr1=(      xb_vecN_2xf32 *)(z+SA*l+k*N   );
        aRd0=BBE_LAN_2XF32_PP(pRd0);
        aRd1=BBE_LAN_2XF32_PP(pRd1);
        for(n=0; n<NN; n++) 
        {
            BBE_LAN_2XF32_IP(t0,aRd0,pRd0);
            BBE_LAN_2XF32_IP(t1,aRd1,pRd1);
            BBE_SAN_2XF32_IP(t1,aWr0,pWr0);
            BBE_SAN_2XF32_IP(t0,aWr1,pWr1);
        }
        BBE_LAVN_2XF32_XP(t0,aRd0,pRd0,nbytes);
        BBE_LAVN_2XF32_XP(t1,aRd1,pRd1,nbytes);
        BBE_SAVN_2XF32_XP(t1,aWr0,pWr0,nbytes);
        BBE_SAVN_2XF32_XP(t0,aWr1,pWr1,nbytes);
        BBE_SAN_2XF32POS_FP(aWr0,pWr0);
        BBE_SAN_2XF32POS_FP(aWr1,pWr1);
    }
}
#endif

static void matinvgjnxnnf_swapCols(float32_t * z, const int16_t *scol, int k, int N, int L, int SA)
#if 0
{
    float32_t t;
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
        pcol=scol[0]; scol+=N;
        pRd0=(const xtfloat *)(z+SA*l+pcol);
        pWr0=(      xtfloat *)(z+SA*l+pcol);
        /* swap columns */
        __Pragma("loop_count min=4")
        for(n=0; n<N; n++) 
        {
            BBE_LSN_2XF32_XP(t0,pRd0,N*sizeof(float32_t));
            BBE_LSN_2XF32_XP(t1,pRd1,N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(t1,pWr0,N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(t0,pWr1,N*sizeof(float32_t));
        }
        pRd1+=SA-N*N;
        pWr1+=SA-N*N;
    }
}
#endif

/* pivoting: search the absolute maximum and its position in the matrix 
    at all positions excluding previously used 
*/
inline_  void matinvgjnxnnf_searchPivot(int16_t *srow,int16_t *scol,const float32_t * z,int k,int N,int L,int SA)
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
#else
{
    short * restrict pCol=(short *)(scol+k);
    short * restrict pRow=(short *)(srow+k);
    xb_vecN_2xc16 idx,maxidx,addidx;
    const xb_vecN_2xf32 * restrict pZ;
    valign aZ;
    int l,m,n,nbytes,NN;
    nbytes = ((N-k) * sizeof(float32_t)) & (2*BBE_SIMD_WIDTH-1);
    NN = (N-k)>>(LOG2_BBE_SIMD_WIDTH-1);
    addidx = BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,BBE_SIMD_WIDTH/2,BBE_SELI_INTERLEAVE_1_LO));
    for (l=0; l<L; l++,z+=SA)
    {
        xtfloat maxT;
        vselN vmax;
        int dummy;
        xb_vecN_2xf32 MAXVAL,T;
        xb_vecNx16 pivot;
        vboolN_2 bmax;
        vboolN_2 cond;
        MAXVAL=FLT_MIN;
        for (m=k; m<N; m++)
        {
            idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(m,BBE_ADDNX16(BBE_SEQNX16(),k),BBE_SELI_INTERLEAVE_1_LO));
            pZ=(const xb_vecN_2xf32*)(&z[m*N+k]);
            aZ=BBE_LAN_2XF32_PP(pZ);
            for(n=0; n<NN; n++) 
            {
                BBE_LAN_2XF32_IP(T,aZ,pZ);
                T=BBE_ABSN_2XF32(T);
                cond=BBE_OGTN_2XF32(T,MAXVAL);
                MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
                maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
                idx=BBE_ADDN_2XC16(idx,addidx);
            }
            BBE_LAVN_2XF32_XP(T,aZ,pZ,nbytes);
            T=BBE_ABSN_2XF32(T);
            cond=BBE_OGTN_2XF32(T,MAXVAL);
            MAXVAL=BBE_MAXNUMN_2XF32(T,MAXVAL);
            maxidx=BBE_MOVN_2XC16T(idx,maxidx,cond);
        }
        BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXVAL);
        (void)maxT;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax)); (void)dummy;
        pivot=BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax);
        pivot=BBE_MAXNX16(0,BBE_MINUNX16(pivot,N-1));
        BBE_SSNX16_XP(pivot,pCol,N*sizeof(int16_t));
        BBE_SSNX16_XP(BBE_REPNX16(pivot,1),pRow,N*sizeof(int16_t));
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
void matinvgjnxnnf ( void* pScr, float32_t * restrict z, int N, int L  )
{
    int16_t *srow,*scol;
    float32_t * temp;   // [N*L]
    int k,SA;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(N%4==0 && N>0);

    if(L<=0 || N<=0) return;

    temp=(float32_t*)(pScr);
    srow=(int16_t*)(temp+N*L);
    scol=srow+N*L;
    NASSERT_ALIGN32(temp);

    SA=getSpace(N*N);
    for (k=0; k<N; k++)
    {
        matinvgjnxnnf_searchPivot(srow,scol,z,k,N,L,SA);
        matinvgjnxnnf_swapRows(z,srow,k,N,L,SA);
        matinvgjnxnnf_swapCols(z,scol,k,N,L, SA);
        matinvgjnxnnf_recipPivot(temp, (z+k*N)+k, L, SA);
        matinvgjnxnnf_updPivot( z+k*N, temp, N, L, SA);
        matinvgjnxnnf_prepelimGS(temp,z,k,N,L,SA);
        matinvgjnxnnf_elimGS(z,z+k*N, temp,N,L,SA);
    }
    /* final reverse permulation of columns/rows  */
    for (k=N-2; k>=0; k--)
    {
        matinvgjnxnnf_swapCols(z,srow,k,N,L, SA);
        matinvgjnxnnf_swapRows(z,scol,k,N,L,SA);
    }

}

size_t matinvgjnxnnf_getScratchSize (int N, int L  )
{
    NASSERT(N%4==0 && N>0);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int16_t)*N*L*2+sizeof(float32_t)*N*L;
}

#else
DISCARD_FUN(void, matinvgjnxnnf,  (void* pScr, float32_t * restrict z, int N, int L  ))

size_t matinvgjnxnnf_getScratchSize (int N, int L  )
{
    NASSERT(N%4==0 && N>0);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif

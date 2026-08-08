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
 * NatureDSP_Baseband Library API
 * Matrix Decomposition and Inversion Functions

    LU decomposition for complex matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include <math.h>
#include <float.h>

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
#endif

/*-------------------------------------------------------------------------
LU Decomposition For Block Ordered Matrices

Description: compute LU decomposition of a square matrix using partial pivoting
with row interchanges: P*A = L*U, where P is the permutation matrix, L is the 
lower triangular factor (with diagnoal elements equal to 1), U is the upper 
triangular factor.

Algorithm is applied in-place to a sequence of real/complex matrices stored
in block order. For each input matrix A[k], the resutling factor U[k] replaces
the upper triangle of input matrix A[k], and subdiagonal elements of factor L[k]
replace the lower triangle of matrix A[k]. Ones on the main diagonal of resulting
factor L[k] are discarded.

Row ordering used by LU decomposition algorithm for input matrix A[k] is stored
to the vector of permutation indices P[k][0..N-1]: i-th row of the matrix was 
interchanged with row P[k][i], i=0..N-1.

Decomposition result is not defined for a close to singular input matrix.

Storage sizes SA,SP denote the number of data elements required to store a
matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)
SP = storage_size(N)

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by 
            functions [c]lu<size>nf_getScratchSize(N,L)
Input:
  N         Matrix size
  L         Number of matrices
Input/Output:
  A[L][SA]  Input matrices, packed L and U factors on output
Output:
  P[L][SP]  Permutation index vectors
Restrictions:
  pScr,A,P  Must not overlap and must be aligned on 32-byte boundary 
  N         Must be a positive multiple of 4
---------------------------------------------------------------------------*/
void clunxnnf ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int N, int L )
#if 0
{
    int SA, SP;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0) return;

    SA=getSpace(N*N<<1)>>1;
    SP=getSpace(N>>1)<<1;   // additional shift right by 1 takes into account that P sizeof(P[0])=sizeof(A[0])/2


    /*
      reference Matlab code
      function [P,A]=lu_fullpivot_2(A)
      sz=size(A);
      N=sz(1);
      assert(sz(1)==sz(2),'for square matrices only');
      P=feval(class(A),eye(N,N));
      for i=1:N
          % find those raw permutation which maximizes A(i,i) 
          % after the row update
          s=A(i:N,i);
          [amax,imax]=max(abs(s));
          imax=imax+i-1;
          % permutation
          t=A(imax,:);    A(imax,:)=A(i,:);    A(i,:)=t;
          t=P(imax,:);    P(imax,:)=P(i,:);    P(i,:)=t;
          for j=(i+1):N
              A(j,i)=A(j,i)/A(i,i);
              k=(i+1:N);
              A(j,k) = A(j,k) - A(j,i)*A(i,k);
          end
      end
      % L=tril(A,-1)+eye(N,N);
      % U=triu(A);
    */
    for (l=0; l<L; l++)
    {
        for (i=0; i<N; i++) P[i]=i;
        for (i=0; i<N-1; i++)
        {
            float32_t maxVal, s;
            complex_float t;
            int pivot;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            maxVal = FLT_MIN; pivot = 0;
            for (j = i; j<N; j++)
            {
                t = A[j*N + i];
                s = fabs(crealf(t)*crealf(t)+cimagf(t)*cimagf(t));
                if (s>maxVal) { maxVal = s; pivot = j; }
            }
            /* permute i-th and pivot-th rows */
            if (i != pivot)
            {
                for (k = 0; k<N; k++)
                {
                    t = A[pivot*N + k]; A[pivot*N + k] = A[i*N + k]; A[i*N + k] = t;
                }
                { int16_t t; t=P[pivot]; P[pivot] = P[i]; P[i] = t; }
            }

            for (j = i + 1; j < N; j++)
            {
                complex_float norm;
                norm=recipc(A[i*N + i]);
                A[j*N + i] = mulc(A[j*N + i],norm);
                for (k = i + 1; k < N; k++)
                {
                    A[j*N + k] = subc(A[j*N + k],mulc(A[j*N + i],A[i*N + k]));
                }
            }
        }
        A+=SA;
        P+=SP;
    }
}
#else
{
    complex_float * Al;
    complex_float *norms; // [L]
    complex_float *temp0,*temp1;   // [N]
    const xtcomplexfloat * restrict pAii;
          xb_vecN_4xcf32 * restrict pNormwr;
    int SA, SP, K, nbytes;
    int l,i,j,k;
    xb_vecN_4xcf32* restrict pAj;
    xb_vecN_4xcf32* restrict pAi;
    xb_vecN_4xcf32* restrict pT0;
    xb_vecN_4xcf32* restrict pT1;
    xtcomplexfloat* restrict pAji;

    xtcomplexfloat* restrict pNorm;
    const xb_vecN_4xcf32 * restrict pAjk;
    const xb_vecN_4xcf32 * restrict pAik;
          xb_vecN_4xcf32 * restrict pAw;
    vboolN_4 bk,bi;
    vboolN_4 abk[4],abi[4];
    vselN seli;

    xb_vecNx16 * pP;
    valign aP;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0) return;

    SA=getSpace(N*N<<1)>>1;
    SP=getSpace(N>>1)<<1;   // additional shift right by 1 takes into account that P sizeof(P[0])=sizeof(A[0])/2
    // allocate on stack
    {
        size_t szN;
        szN=sizeof(complex_float)*N;
        temp0=(complex_float*)pScr;
        temp1=(complex_float*)(((uintptr_t)temp0)+szN);
        norms=(complex_float*)(((uintptr_t)temp1)+szN);
    }
    NASSERT_ALIGN(temp0,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(temp1,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(norms,2*BBE_SIMD_WIDTH);

    // form bit masks
    abi[0]=BBE_LTRN_4I(1) &~ BBE_LTRN_4I(0);    /* 1 in i-th position */
    abi[1]=BBE_LTRN_4I(2) &~ BBE_LTRN_4I(1);
    abi[2]=BBE_LTRN_4I(3) &~ BBE_LTRN_4I(2);
    abi[3]=~ BBE_LTRN_4I(3);
    abk[0]=~BBE_LTRN_4I(1);                    /* 1 for k>i          */
    abk[1]=~BBE_LTRN_4I(2);                 
    abk[2]=~BBE_LTRN_4I(3);                 
    abk[3]= BBE_LTRN_4I(0); 

    // initialize permutation tables (write SP elements instead of N to fill out all holes between blocks)
    K=SP>>(LOG2_BBE_SIMD_WIDTH);
    nbytes=(SP&(BBE_SIMD_WIDTH-1))*sizeof(int16_t);
    for (l=0; l<L; l++)
    {
        xb_vecNx16 seq,inc;
        pP=(xb_vecNx16*)&P[l*SP];
        aP=BBE_ZALIGN();
        seq=BBE_SEQNX16(); inc=BBE_SIMD_WIDTH;
        for (k=0; k<K; k++)
        {
            BBE_SANX16_IP(seq,aP,pP);
            seq=BBE_ADDNX16(seq,inc);
        }
        BBE_SAVNX16_XP(seq,aP,pP,nbytes);
        BBE_SANX16POS_FP(aP,pP);
    }
    __Pragma("no_reorder")

    for (i=0; i<N-1; i++)
    {
        K=N>>(LOG2_BBE_SIMD_WIDTH-2);
        Al=A;
        for (l=0; l<L; l++,Al+=SA)
        {
            int pivot,dummy;
            xtfloat maxT;
            vselN vmax;
            xb_vecN_2xc16 maxidx,idx,inc;
            vboolN_2 bmask,bmax; /* for i&3==0: 0000, i&3==1: 1000,  i&3==2: 1100. i&3==3: 1110 */
            xb_vecN_2xf32 MAXT;
            (void)dummy;

            pAji=(xtcomplexfloat*)&Al[(i&~3)*N + i];
            bmask=BBE_LTRN_2(i&3);
            MAXT=FLT_MIN;
            idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,BBE_ADDNX16(BBE_SEQNX16(),(i&~3)),BBE_SELI_INTERLEAVE_1_LO));
            inc=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,4,BBE_SELI_INTERLEAVE_1_LO));
            maxidx=idx;
            for (j = (i&~3); j<N; j+=4)
            {
                xb_vecN_4xcf32 t0,t1,t2,t3,a0;
                xb_vecN_2xf32  T0;
                BBE_LSN_4XCF32_XP(t0,pAji,N*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAji,N*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAji,N*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAji,N*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                t0=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);

                a0=BBE_MULMN_4XCF32(t0,t0,0,0); 
                BBE_MULMASN_4XCF32(a0,t0,t0,0,15); 
                t0=BBE_SELN_4XCF32I(BBE_CONSTN_4XCF32(0),a0,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                T0=BBE_MOVN_2XF32_FROMN_4XCF32(t0);

                BBE_CONSTN_2XF32T(T0,0,bmask);
                MAXT=BBE_MAXN_2XF32(MAXT,T0);
                maxidx=BBE_MOVN_2XC16T(idx,maxidx,BBE_OEQN_2XF32(T0,MAXT));
                idx=BBE_ADDN_2XC16(idx,inc);
                bmask=BBE_LTRN_2(0);
            }
            BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXT);
            (void)maxT;
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
            pivot=BBE_EXTRNX16C(BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax),0);

            /* permute i-th and pivot-th rows */
            pAi=(xb_vecN_4xcf32*)&Al[i*N];
            pAj=(xb_vecN_4xcf32*)&Al[pivot*N];
            pT0=(xb_vecN_4xcf32*)temp0;
            pT1=(xb_vecN_4xcf32*)temp1;
            NASSERT_ALIGN(pAi,2*BBE_SIMD_WIDTH);
            NASSERT_ALIGN(pAj,2*BBE_SIMD_WIDTH);
            for (k=0; k<K; k++)
            {
                xb_vecN_4xcf32 t0,t1;
                BBE_LVN_4XCF32_IP(t0,pAj,2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_IP(t1,pAi,2*BBE_SIMD_WIDTH);
                BBE_SVN_4XCF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
                BBE_SVN_4XCF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
            }
            pAi=(xb_vecN_4xcf32*)&Al[i*N];
            pAj=(xb_vecN_4xcf32*)&Al[pivot*N];
            pT0=(xb_vecN_4xcf32*)temp0;
            pT1=(xb_vecN_4xcf32*)temp1;
            __Pragma("no_reorder")
            for (k=0; k<K; k++)
            {
                xb_vecN_4xcf32 t0,t1;
                BBE_LVN_4XCF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
                BBE_SVN_4XCF32_IP(t1,pAj,2*BBE_SIMD_WIDTH);
                BBE_SVN_4XCF32_IP(t0,pAi,2*BBE_SIMD_WIDTH);
            }
            { int16_t t; t=P[l*SP+pivot]; P[l*SP+pivot] = P[l*SP+i]; P[l*SP+i] = t; }
            __Pragma("no_reorder")
        }
        // compute normalization factor
        K=L>>(LOG2_BBE_SIMD_WIDTH-2);
        nbytes=(L&(BBE_SIMD_WIDTH/4-1))*sizeof(complex_float);
        pAii=(const xtcomplexfloat*)&A[i*N + i];
        pNormwr=(xb_vecN_4xcf32*)norms;
        for (l=0; l<K; l++)
        {
            xb_vecN_4xcf32 t0,t1,t2,t3,T;
            BBE_LSN_4XCF32_XP(t0,pAii,SA*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t1,pAii,SA*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t2,pAii,SA*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t3,pAii,SA*sizeof(complex_float));
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T =BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            T =BBE_RECIPN_4XCF32(T);
            BBE_SVN_4XCF32_IP(T,pNormwr,2*BBE_SIMD_WIDTH);
        }
        if (nbytes)
        {
            int Alast=(int)(uintptr_t)(A+L*SA-1);
            xb_vecN_4xcf32 t0,t1,t2,t3,T;
            valign aN=BBE_ZALIGN();
            BBE_LSN_4XCF32_XP(t0,pAii,SA*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t1,pAii,SA*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t2,pAii,SA*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t3,pAii,SA*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T =BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            T =BBE_RECIPN_4XCF32(T);
            BBE_SAVN_4XCF32_XP(T,aN,pNormwr,nbytes);
            BBE_SAN_4XCF32POS_FP(aN,pNormwr);
        }
        __Pragma("no_reorder")
        bi=abi[i&3];
        bk=abk[i&3];
        /* form select register for duplication of i-th element */
        {
            xb_vecNx16 t;
            int dummy;
            (void)dummy;
            BBE_SQZN(seli,dummy,BBE_MOVN_FROMN_4(bi));
            t=BBE_SELNX16(BBE_SEQNX16(),BBE_SEQNX16(),seli);
            seli=BBE_MOVVSV(BBE_MOVNX16_FROMN_4X64(BBE_REPN_4X64(BBE_MOVN_4X64_FROMNX16(t),0)),0);
        }
        /* update rows below current */
        K=((N-(i&~3))>>(LOG2_BBE_SIMD_WIDTH-2));
        pNorm=(xtcomplexfloat*)norms;
        pAik=(const xb_vecN_4xcf32*)&A[i*N + (i&~3)];
        pAjk=(const xb_vecN_4xcf32*)&A[((i+1))*N + (i&~3)];
        pAw =(xb_vecN_4xcf32*)pAjk;
        switch(K)
        {
        case 1:
            for (l=0; l<L; l++)
            {
                xb_vecN_4xcf32 Aji,NORM,Ajk,Aik;
                BBE_LVN_4XCF32_XP(Aik,pAik,SA*sizeof(complex_float));
                BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
                NORM=BBE_REPN_4XCF32(NORM,0);
                for (j = 1; j < N-i; j++)
                {
                    BBE_LVN_4XCF32_XP(Ajk,pAjk,N*sizeof(complex_float));
                    BBE_MULN_4XCF32T(Ajk,Ajk,NORM,bi);
                    Aji=BBE_SELN_4XCF32(Ajk,Ajk,seli);
                    BBE_MULSN_4XCF32T(Ajk,Aji,Aik,bk);
                    BBE_SVN_4XCF32_XP(Ajk,pAw,N*sizeof(complex_float));
                }
                pAjk+=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAw +=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            break;
        case 2:
            for (l=0; l<L; l++)
            {
                xb_vecN_4xcf32 Aji,NORM,Ajk0,Ajk1,Aik0,Aik1;
                BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
                Aik1=BBE_LVN_4XCF32_I (pAik,1*2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_XP(Aik0,pAik,SA*sizeof(complex_float));
                NORM=BBE_REPN_4XCF32(NORM,0);
                for (j = 1; j < N-i; j++)
                {
                    Ajk1=BBE_LVN_4XCF32_I (pAjk,1*2*BBE_SIMD_WIDTH);
                    BBE_LVN_4XCF32_XP(Ajk0,pAjk,N*sizeof(complex_float));
                    BBE_MULN_4XCF32T(Ajk0,Ajk0,NORM,bi);
                    Aji=BBE_SELN_4XCF32(Ajk0,Ajk0,seli);
                    BBE_MULSN_4XCF32T(Ajk0,Aji,Aik0,bk);
                    BBE_MULSN_4XCF32(Ajk1,Aji,Aik1);
                    BBE_SVN_4XCF32_I (Ajk1,pAw,1*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_XP(Ajk0,pAw,N*sizeof(complex_float));
                }
                pAjk+=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAw +=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            break;
        case 3:
            for (l=0; l<L; l++)
            {
                xb_vecN_4xcf32 Aji,NORM,Ajk0,Ajk1,Ajk2,Aik0,Aik1,Aik2;
                BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
                Aik1=BBE_LVN_4XCF32_I (pAik,1*2*BBE_SIMD_WIDTH);
                Aik2=BBE_LVN_4XCF32_I (pAik,2*2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_XP(Aik0,pAik,SA*sizeof(complex_float));
                NORM=BBE_REPN_4XCF32(NORM,0);
                for (j = 1; j < N-i; j++)
                {
                    Ajk1=BBE_LVN_4XCF32_I (pAjk,1*2*BBE_SIMD_WIDTH);
                    Ajk2=BBE_LVN_4XCF32_I (pAjk,2*2*BBE_SIMD_WIDTH);
                    BBE_LVN_4XCF32_XP(Ajk0,pAjk,N*sizeof(complex_float));
                    BBE_MULN_4XCF32T(Ajk0,Ajk0,NORM,bi);
                    Aji=BBE_SELN_4XCF32(Ajk0,Ajk0,seli);
                    BBE_MULSN_4XCF32T(Ajk0,Aji,Aik0,bk);
                    BBE_MULSN_4XCF32(Ajk1,Aji,Aik1);
                    BBE_MULSN_4XCF32(Ajk2,Aji,Aik2);
                    BBE_SVN_4XCF32_I (Ajk1,pAw,1*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_I (Ajk2,pAw,2*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_XP(Ajk0,pAw,N*sizeof(complex_float));
                }
                pAjk+=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAw +=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            break;
        case 4:
            for (l=0; l<L; l++)
            {
                xb_vecN_4xcf32 Aji,NORM,Ajk0,Ajk1,Ajk2,Ajk3,Aik0,Aik1,Aik2,Aik3;
                BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
                Aik1=BBE_LVN_4XCF32_I (pAik,1*2*BBE_SIMD_WIDTH);
                Aik2=BBE_LVN_4XCF32_I (pAik,2*2*BBE_SIMD_WIDTH);
                Aik3=BBE_LVN_4XCF32_I (pAik,3*2*BBE_SIMD_WIDTH);
                BBE_LVN_4XCF32_XP(Aik0,pAik,SA*sizeof(complex_float));
                NORM=BBE_REPN_4XCF32(NORM,0);
                for (j = 1; j < N-i; j++)
                {
                    Ajk1=BBE_LVN_4XCF32_I (pAjk,1*2*BBE_SIMD_WIDTH);
                    Ajk2=BBE_LVN_4XCF32_I (pAjk,2*2*BBE_SIMD_WIDTH);
                    Ajk3=BBE_LVN_4XCF32_I (pAjk,3*2*BBE_SIMD_WIDTH);
                    BBE_LVN_4XCF32_XP(Ajk0,pAjk,N*sizeof(complex_float));
                    BBE_MULN_4XCF32T(Ajk0,Ajk0,NORM,bi);
                    Aji=BBE_SELN_4XCF32(Ajk0,Ajk0,seli);
                    BBE_MULSN_4XCF32T(Ajk0,Aji,Aik0,bk);
                    BBE_MULSN_4XCF32(Ajk1,Aji,Aik1);
                    BBE_MULSN_4XCF32(Ajk2,Aji,Aik2);
                    BBE_MULSN_4XCF32(Ajk3,Aji,Aik3);
                    BBE_SVN_4XCF32_I (Ajk1,pAw,1*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_I (Ajk2,pAw,2*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_I (Ajk3,pAw,3*2*BBE_SIMD_WIDTH);
                    BBE_SVN_4XCF32_XP(Ajk0,pAw,N*sizeof(complex_float));
                }
                pAjk+=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAw +=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            break;
        default:
            for (l=0; l<L; l++)
            {
                xb_vecN_4xcf32 Aji,NORM,Ajk,Aik;
                BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
                NORM=BBE_REPN_4XCF32(NORM,0);
                for (j = 1; j < N-i; j++)
                {
                    BBE_LVN_4XCF32_IP(Ajk,pAjk,2*BBE_SIMD_WIDTH);
                    BBE_LVN_4XCF32_IP(Aik,pAik,2*BBE_SIMD_WIDTH);
                    BBE_MULN_4XCF32T(Ajk,Ajk,NORM,bi);
                    Aji=BBE_SELN_4XCF32(Ajk,Ajk,seli);
                    BBE_MULSN_4XCF32T(Ajk,Aji,Aik,bk);
                    BBE_SVN_4XCF32_IP(Ajk,pAw,2*BBE_SIMD_WIDTH);
                    for (k = 1; k < K; k++)
                    {
                        BBE_LVN_4XCF32_IP(Ajk,pAjk,2*BBE_SIMD_WIDTH);
                        BBE_LVN_4XCF32_IP(Aik,pAik,2*BBE_SIMD_WIDTH);
                        BBE_MULSN_4XCF32(Ajk,Aji,Aik);
                        BBE_SVN_4XCF32_IP(Ajk,pAw,2*BBE_SIMD_WIDTH);
                    }
                    pAjk+=-K+(N*sizeof(complex_float)/sizeof(xb_vecN_4xcf32));
                    pAw +=-K+(N*sizeof(complex_float)/sizeof(xb_vecN_4xcf32));
                    pAik+=-K;
                }
                pAjk+=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAw +=((SA-(N-i-1)*N)*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
                pAik+=(SA*sizeof(complex_float))/sizeof(xb_vecN_4xcf32);
            }
            break;
        }
    }
}
#endif

/* Return the scratch area size, in bytes. */
size_t clunxnnf_getScratchSize   ( int N, int L )
{
    size_t szL,szN;
    (void)N, (void)L;
    NASSERT(N>0 && N%4==0);
    L=XT_MAX(L,0);
    N=XT_MAX(N,0);
    szL=sizeof(complex_float)*L;
    szN=sizeof(complex_float)*N;
    szL=(szL+2*BBE_SIMD_WIDTH  )&~(2*BBE_SIMD_WIDTH-1);
    return szL+2*szN;
}
#else
DISCARD_FUN(void, clunxnnf, ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int N, int L ))

size_t clunxnnf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    if (L<=0) return 0;
    return 0;
}
#endif

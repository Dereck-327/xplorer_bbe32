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

    LU decomposition for real matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include <math.h>
#include <float.h>

#if HAVE_VFPU

#undef BBE_SELN_2XF32
#define BBE_SELN_2XF32(b,c,d) BBE_MOVN_2XF32_FROMNX16((BBE_SELNX16((BBE_MOVNX16_FROMN_2XF32(b)),(BBE_MOVNX16_FROMN_2XF32(c)),d)))


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
void lunxnnf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int N, int L )
#if 0
{
    float32_t *norms; // [L]
    float32_t *temp0,*temp1;   // [N]
    int SA,SP;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0 || N<=1) return;

    SA=getSpace(N*N);
    SP=getSpace(N>>1)<<1;   // additional shift right by 1 takes into account that P sizeof(P[0])=sizeof(A[0])/2

    // allocate on stack
    {
        size_t szL,szN;
        szL=sizeof(float32_t)*L;
        szN=sizeof(float32_t)*N;
        szN=(szN+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
        temp0=(float32_t*)pScr;
        temp1=(float32_t*)(((uintptr_t)temp0)+szN);
        norms=(float32_t*)(((uintptr_t)temp1)+szN);
    }
    NASSERT_ALIGN(temp0,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(temp1,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(norms,2*BBE_SIMD_WIDTH);


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
        for (i=0; i<N; i++) P[l*SP+i]=i;
    }
    for (i=0; i<N-1; i++)
    {
        for (l=0; l<L; l++)
        {
            float32_t maxVal, s, t;
            int pivot;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            maxVal = FLT_MIN; pivot = 0;
            for (j = i; j<N; j++)
            {
                s = A[l*SA+j*N + i];
                s = fabs(s);
                if (s>maxVal) { maxVal = s; pivot = j; }
            }
            /* permute i-th and pivot-th rows */
            for (k = 0; k<N; k++)
            {
                temp0[k]=A[l*SA+pivot*N + k];
                temp1[k]=A[l*SA+i*N + k];
            }
            for (k = 0; k<N; k++)
            {
                A[l*SA+pivot*N + k]=temp1[k];
                A[l*SA+i*N + k]    =temp0[k];
            }

            { int16_t t; t=P[l*SP+pivot]; P[l*SP+pivot] = P[l*SP+i]; P[l*SP+i] = t; }
        }
        for (l=0; l<L; l++)
        {
            norms[l]=1.0f/A[l*SA+i*N + i];
        }
        for (l=0; l<L; l++)
        {
            float32_t norm=norms[l];
            for (j = i + 1; j < N; j++)
            {
                A[l*SA+j*N + i] *= norm;
                for (k = i + 1; k < N; k++)
                {
                    A[l*SA+j*N + k] -= A[l*SA+j*N + i] * A[l*SA+i*N + k];
                }
            }
        }
    }
}
#else

{
    float32_t *norms; // [L]
    float32_t *temp0,*temp1;   // [N]
    float32_t *Al;
    xtfloat* restrict pNorm;
    xtfloat* restrict pAii;
    const xb_vecN_2xf32* restrict pNormr;
          xb_vecN_2xf32* restrict pNormw;
            xtfloat      * restrict pAji;
    const xb_vecN_2xf32* restrict pAjk;
    const xb_vecN_2xf32* restrict pAik;
            xb_vecN_2xf32* restrict pAw;
    int K,nbytes;
    xb_vecN_2xf32* restrict pAj;
    xb_vecN_2xf32* restrict pAi;
    xb_vecN_2xf32* restrict pT0;
    xb_vecN_2xf32* restrict pT1;
    xb_vecNx16 * pP;
    valign aj,ai,aP;
    vboolN_2 bi,bk,abi[8],abk[8];
    vselN seli;

    int SA,SP;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0 || N<=1) return;

    SA=getSpace(N*N);
    SP=getSpace(N>>1)<<1;   // additional shift right by 1 takes into account that P sizeof(P[0])=sizeof(A[0])/2

    // allocate on stack
    {
        size_t szN;
        szN=sizeof(float32_t)*N;
        szN=(szN+2*BBE_SIMD_WIDTH)&~(2*BBE_SIMD_WIDTH-1);
        temp0=(float32_t*)pScr;
        temp1=(float32_t*)(((uintptr_t)temp0)+szN);
        norms=(float32_t*)(((uintptr_t)temp1)+szN);
    }
    NASSERT_ALIGN(temp0,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(temp1,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(norms,2*BBE_SIMD_WIDTH);
    // form bit masks
    abi[0]= BBE_LTRN_2I(1) &~ BBE_LTRN_2I(0);    /* 1 in i-th position */
    abi[1]= BBE_LTRN_2I(2) &~ BBE_LTRN_2I(1);
    abi[2]= BBE_LTRN_2I(3) &~ BBE_LTRN_2I(2);
    abi[3]= BBE_LTRN_2I(4) &~ BBE_LTRN_2I(3);
    abi[4]= BBE_LTRN_2I(5) &~ BBE_LTRN_2I(4);
    abi[5]= BBE_LTRN_2I(6) &~ BBE_LTRN_2I(5);
    abi[6]= BBE_LTRN_2I(7) &~ BBE_LTRN_2I(6);
    abi[7]=~BBE_LTRN_2I(7);
    abk[0]=~BBE_LTRN_2I(1);                    /* 1 for k>i          */
    abk[1]=~BBE_LTRN_2I(2);                 
    abk[2]=~BBE_LTRN_2I(3);                 
    abk[3]=~BBE_LTRN_2I(4);                 
    abk[4]=~BBE_LTRN_2I(5);                 
    abk[5]=~BBE_LTRN_2I(6);                 
    abk[6]=~BBE_LTRN_2I(7);                 
    abk[7]= BBE_LTRN_2I(0); 

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

    for (i=0; i<N-1; i++)
    {
        K=N>>(LOG2_BBE_SIMD_WIDTH-1);
        nbytes=(N&(BBE_SIMD_WIDTH/2-1))*sizeof(float32_t);
        Al=A;
        for (l=0; l<L; l++,Al+=SA)
        {
            xtfloat maxT;
            int pivot;
            vselN vmax;
            xb_vecN_2xc16 maxidx,idx,inc;
            vboolN_2 bmask,bmax; /* for i&3==0: 0000, i&3==1: 1000,  i&3==2: 1100. i&3==3: 1110 */
            xb_vecN_2xf32 t0,t1,t2,t3,MAXT;
            int dummy;
            (void)dummy;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            pAji=(xtfloat*)&Al[(i&~3)*N + i];
            bmask=BBE_LTRN_2(i&3);
            MAXT=FLT_MIN;
            idx=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,BBE_ADDNX16(BBE_SEQNX16(),(i&~3)),BBE_SELI_INTERLEAVE_1_LO));
            inc=BBE_MOVN_2XC16_FROMNX16(BBE_SELNX16I(0,4,BBE_SELI_INTERLEAVE_1_LO));
            maxidx=idx;
            for (j = (i&~3); j<N; j+=4)
            {
                BBE_LSN_2XF32_XP(t0,pAji,N*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t1,pAji,N*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t2,pAji,N*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t3,pAji,N*sizeof(float32_t));
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t0=BBE_ABSN_2XF32(t0);
                BBE_CONSTN_2XF32T(t0,0,bmask);
                MAXT=BBE_MAXN_2XF32(MAXT,t0);
                maxidx=BBE_MOVN_2XC16T(idx,maxidx,BBE_OEQN_2XF32(t0,MAXT));
                idx=BBE_ADDN_2XC16(idx,inc);
                bmask=BBE_LTRN_2(0);
            }
            BBE_RBMAXNUMN_2XF32(bmax,maxT,MAXT);
            (void)maxT;
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
            pivot=BBE_EXTRNX16C(BBE_SELNX16(BBE_MOVNX16_FROMN_2XC16(maxidx),BBE_MOVNX16_FROMN_2XC16(maxidx),vmax),0);

            /* permute i-th and pivot-th rows */
            pAi=(xb_vecN_2xf32*)&Al[i*N];
            pAj=(xb_vecN_2xf32*)&Al[pivot*N];
            pT0=(xb_vecN_2xf32*)temp0;
            pT1=(xb_vecN_2xf32*)temp1;
            ai=BBE_LAN_2XF32_PP(pAi);
            aj=BBE_LAN_2XF32_PP(pAj);
            for (k=0; k<K; k++)
            {
                BBE_LAN_2XF32_IP(t0,aj,pAj);
                BBE_LAN_2XF32_IP(t1,ai,pAi);
                BBE_SVN_2XF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
            }
            BBE_LAVN_2XF32_XP(t0,aj,pAj,nbytes);
            BBE_LAVN_2XF32_XP(t1,ai,pAi,nbytes);
            BBE_SVN_2XF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
            ai=aj=BBE_ZALIGN();
            pAi=(xb_vecN_2xf32*)&Al[i*N];
            pAj=(xb_vecN_2xf32*)&Al[pivot*N];
            pT0=(xb_vecN_2xf32*)temp0;
            pT1=(xb_vecN_2xf32*)temp1;
            __Pragma("no_reorder")
            for (k=0; k<K; k++)
            {
                BBE_LVN_2XF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
                BBE_LVN_2XF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
                BBE_SAN_2XF32_IP(t1,aj,pAj);
                BBE_SAN_2XF32_IP(t0,ai,pAi);
            }
            BBE_LVN_2XF32_IP(t0,pT0,2*BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(t1,pT1,2*BBE_SIMD_WIDTH);
            BBE_SAVN_2XF32_XP(t1,aj,pAj,nbytes);
            BBE_SAVN_2XF32_XP(t0,ai,pAi,nbytes);
            BBE_SAN_2XF32POS_FP(aj,pAj);
            BBE_SAN_2XF32POS_FP(ai,pAi);
            { int16_t t; t=P[l*SP+pivot]; P[l*SP+pivot] = P[l*SP+i]; P[l*SP+i] = t; }
        }
        // compute normalization factors
        __Pragma("no_reorder")
        pNorm=(xtfloat*)norms;
        pAii=(xtfloat*)&A[i*N+i];
        for (l=0; l<L; l++)
        {
            xb_vecN_2xf32 Ai;
            BBE_LSN_2XF32_XP(Ai,pAii,SA*sizeof(float32_t));
            BBE_SSN_2XF32_IP(Ai,pNorm,sizeof(float32_t));
        }
        // update rows below current
        __Pragma("no_reorder")
        pNormr=(const xb_vecN_2xf32*)norms;
        pNormw=(      xb_vecN_2xf32*)norms;
        for (l=0; l<((L+BBE_SIMD_WIDTH/2-1)>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
        {
            xb_vecN_2xf32 norm;
            BBE_LVN_2XF32_IP(norm,pNormr,2*BBE_SIMD_WIDTH);
            norm=BBE_RECIPN_2XF32(norm);
            BBE_SVN_2XF32_IP(norm,pNormw,2*BBE_SIMD_WIDTH);
        }
        __Pragma("no_reorder")

        /*
        optimization note: comparing with complex routines, we can not avoid unaligned access because 
        N is a multiple of 4, so rows in the matrices not neccesarily to be aligned by 
        2*BBE_SIMD_WIDTH==32 bytes
        So, here we hav 2 algorithms: one with unaligned access for N not a multiple of 8 and another 
        one - for aligning access 
        */
        if (N&7)
        {
            pAji=(xtfloat *)&A[(i + 1)*N+i] ;
            pNorm=(xtfloat*)norms;
            K=((N-i -1)>>(LOG2_BBE_SIMD_WIDTH-1));
            nbytes=((N-i-1)&(BBE_SIMD_WIDTH/2-1))*sizeof(float32_t);
            for (Al=A,l=0; l<L; l++,Al+=SA)
            {
                xb_vecN_2xf32 Aji,NORM,Ajk,Aik;
                valign ajk,aik,aw;
                BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
                switch (K)  // specialized variants for small K
                {
                case 0:
                    for (j = i + 1; j < N; j++)
                    {
                        Aji=BBE_LSN_2XF32_I(pAji,0);
                        Aji=BBE_MULN_2XF32(Aji,NORM);
                        BBE_SSN_2XF32_XP(Aji,pAji,N*sizeof(float32_t));
                        Aji=BBE_REPN_2XF32(Aji,0);
                        aw=BBE_ZALIGN();
                        pAjk=(const xb_vecN_2xf32*)&Al[j*N + i + 1];
                        pAw =(      xb_vecN_2xf32*)pAjk;
                        pAik=(const xb_vecN_2xf32*)&Al[i*N + i + 1];
                        ajk=BBE_LAN_2XF32_PP(pAjk);
                        aik=BBE_LAN_2XF32_PP(pAik);
                        BBE_LAVN_2XF32_XP(Ajk,ajk,pAjk,nbytes);
                        BBE_LAVN_2XF32_XP(Aik,aik,pAik,nbytes);
                        BBE_MULSN_2XF32(Ajk,Aji,Aik);
                        BBE_SAVN_2XF32_XP(Ajk,aw ,pAw,nbytes);
                        BBE_SAN_2XF32POS_FP(aw ,pAw);
                    }
                    break;
                case 1:
                    for (j = i + 1; j < N; j++)
                    {
                        Aji=BBE_LSN_2XF32_I(pAji,0);
                        Aji=BBE_MULN_2XF32(Aji,NORM);
                        BBE_SSN_2XF32_XP(Aji,pAji,N*sizeof(float32_t));
                        Aji=BBE_REPN_2XF32(Aji,0);
                        aw=BBE_ZALIGN();
                        pAjk=(const xb_vecN_2xf32*)&Al[j*N + i + 1];
                        pAw =(      xb_vecN_2xf32*)pAjk;
                        pAik=(const xb_vecN_2xf32*)&Al[i*N + i + 1];
                        ajk=BBE_LAN_2XF32_PP(pAjk);
                        aik=BBE_LAN_2XF32_PP(pAik);
                        BBE_LAN_2XF32_IP(Ajk,ajk,pAjk);
                        BBE_LAN_2XF32_IP(Aik,aik,pAik);
                        BBE_MULSN_2XF32(Ajk,Aji,Aik);
                        BBE_SAN_2XF32_IP(Ajk,aw ,pAw);
                        BBE_LAVN_2XF32_XP(Ajk,ajk,pAjk,nbytes);
                        BBE_LAVN_2XF32_XP(Aik,aik,pAik,nbytes);
                        BBE_MULSN_2XF32(Ajk,Aji,Aik);
                        BBE_SAVN_2XF32_XP(Ajk,aw ,pAw,nbytes);
                        BBE_SAN_2XF32POS_FP(aw ,pAw);
                    }
                    break;
                default:
                    for (j = i + 1; j < N; j++)
                    {
                        Aji=BBE_LSN_2XF32_I(pAji,0);
                        Aji=BBE_MULN_2XF32(Aji,NORM);
                        BBE_SSN_2XF32_XP(Aji,pAji,N*sizeof(float32_t));
                        Aji=BBE_REPN_2XF32(Aji,0);
                        aw=BBE_ZALIGN();
                        pAjk=(const xb_vecN_2xf32*)&Al[j*N + i + 1];
                        pAw =(      xb_vecN_2xf32*)pAjk;
                        pAik=(const xb_vecN_2xf32*)&Al[i*N + i + 1];
                        ajk=BBE_LAN_2XF32_PP(pAjk);
                        aik=BBE_LAN_2XF32_PP(pAik);
                        for (k = 0; k < K; k++)
                        {
                            BBE_LAN_2XF32_IP(Ajk,ajk,pAjk);
                            BBE_LAN_2XF32_IP(Aik,aik,pAik);
                            BBE_MULSN_2XF32(Ajk,Aji,Aik);
                            BBE_SAN_2XF32_IP(Ajk,aw ,pAw);
                        }
                        BBE_LAVN_2XF32_XP(Ajk,ajk,pAjk,nbytes);
                        BBE_LAVN_2XF32_XP(Aik,aik,pAik,nbytes);
                        BBE_MULSN_2XF32(Ajk,Aji,Aik);
                        BBE_SAVN_2XF32_XP(Ajk,aw ,pAw,nbytes);
                        BBE_SAN_2XF32POS_FP(aw ,pAw);
                    }
                    break;
                }
                pAji += (SA-(N-i-1)*N);
            }
        }
        else // N is multiple of 8: all rows will be aligned - another algorithm
        {
            bi=abi[i&7];
            bk=abk[i&7];
            /* form select register for duplication of i-th element */
            {
                xb_vecNx16 t;
                int dummy;
                (void)dummy;
                BBE_SQZN(seli,dummy,BBE_MOVN_FROMN_2(bi));
                t=BBE_SELNX16(BBE_SEQNX16(),BBE_SEQNX16(),seli);
                seli=BBE_MOVVSV(BBE_MOVNX16_FROMN_2X32(BBE_REPN_2X32(BBE_MOVN_2X32_FROMNX16(t),0)),0);
            }
            /* update rows below current */
            K=((N-(i&~7))>>(LOG2_BBE_SIMD_WIDTH-1));
            pNorm=(xtfloat*)norms;
            pAik=(const xb_vecN_2xf32*)&A[i*N + (i&~7)];
            pAjk=(const xb_vecN_2xf32*)&A[((i+1))*N + (i&~7)];
            pAw =(xb_vecN_2xf32*)pAjk;
            switch(K)
            {
            case 1:
                for (l=0; l<L; l++)
                {
                    xb_vecN_2xf32 Aji,NORM,Ajk,Aik;
                    BBE_LVN_2XF32_XP(Aik,pAik,SA*sizeof(float32_t));
                    BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
                    NORM=BBE_REPN_2XF32(NORM,0);
                    for (j = 1; j < N-i; j++)
                    {
                        BBE_LVN_2XF32_XP(Ajk,pAjk,N*sizeof(float32_t));
                        BBE_MULN_2XF32T(Ajk,Ajk,NORM,bi);
                        Aji=BBE_SELN_2XF32(Ajk,Ajk,seli);
                        BBE_MULSN_2XF32T(Ajk,Aji,Aik,bk);
                        BBE_SVN_2XF32_XP(Ajk,pAw,N*sizeof(float32_t));
                    }
                    pAjk+=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                    pAw +=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                }
                break;
            case 2:
                for (l=0; l<L; l++)
                {
                    xb_vecN_2xf32 Aji,NORM,Ajk0,Ajk1,Aik0,Aik1;
                    BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
                    Aik1=BBE_LVN_2XF32_I (pAik,1*2*BBE_SIMD_WIDTH);
                    BBE_LVN_2XF32_XP(Aik0,pAik,SA*sizeof(float32_t));
                    NORM=BBE_REPN_2XF32(NORM,0);
                    for (j = 1; j < N-i; j++)
                    {
                        Ajk1=BBE_LVN_2XF32_I (pAjk,1*2*BBE_SIMD_WIDTH);
                        BBE_LVN_2XF32_XP(Ajk0,pAjk,N*sizeof(float32_t));
                        BBE_MULN_2XF32T(Ajk0,Ajk0,NORM,bi);
                        Aji=BBE_SELN_2XF32(Ajk0,Ajk0,seli);
                        BBE_MULSN_2XF32T(Ajk0,Aji,Aik0,bk);
                        BBE_MULSN_2XF32(Ajk1,Aji,Aik1);
                        BBE_SVN_2XF32_I (Ajk1,pAw,1*2*BBE_SIMD_WIDTH);
                        BBE_SVN_2XF32_XP(Ajk0,pAw,N*sizeof(float32_t));
                    }
                    pAjk+=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                    pAw +=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                }
                break;
            default:
                for (l=0; l<L; l++)
                {
                    xb_vecN_2xf32 Aji,NORM,Ajk,Aik;
                    BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
                    NORM=BBE_REPN_2XF32(NORM,0);
                    for (j = 1; j < N-i; j++)
                    {
                        BBE_LVN_2XF32_IP(Ajk,pAjk,2*BBE_SIMD_WIDTH);
                        BBE_LVN_2XF32_IP(Aik,pAik,2*BBE_SIMD_WIDTH);
                        BBE_MULN_2XF32T(Ajk,Ajk,NORM,bi);
                        Aji=BBE_SELN_2XF32(Ajk,Ajk,seli);
                        BBE_MULSN_2XF32T(Ajk,Aji,Aik,bk);
                        BBE_SVN_2XF32_IP(Ajk,pAw,2*BBE_SIMD_WIDTH);
                        for (k = 1; k < K; k++)
                        {
                            BBE_LVN_2XF32_IP(Ajk,pAjk,2*BBE_SIMD_WIDTH);
                            BBE_LVN_2XF32_IP(Aik,pAik,2*BBE_SIMD_WIDTH);
                            BBE_MULSN_2XF32(Ajk,Aji,Aik);
                            BBE_SVN_2XF32_IP(Ajk,pAw,2*BBE_SIMD_WIDTH);
                        }
                        pAjk+=-K+(N*sizeof(float32_t)/sizeof(xb_vecN_2xf32));
                        pAw +=-K+(N*sizeof(float32_t)/sizeof(xb_vecN_2xf32));
                        pAik+=-K;
                    }
                    pAjk+=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                    pAw +=((SA-(N-i-1)*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                    pAik+=(SA*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                }
                break;
            }
        }
    }
}
#endif

/* Return the scratch area size, in bytes. */
size_t lunxnnf_getScratchSize   ( int N, int L )
{
    size_t szL,szN;
    (void)N, (void)L;
    NASSERT(N>0 && N%4==0);
    L=XT_MAX(L,0);
    N=XT_MAX(N,0);
    szL=sizeof(float32_t)*L;
    szN=sizeof(float32_t)*N;
    szN=(szN+2*BBE_SIMD_WIDTH)&~(2*BBE_SIMD_WIDTH-1);
    szL=(szL+2*BBE_SIMD_WIDTH)&~(2*BBE_SIMD_WIDTH-1);
    return szL+2*szN;
}
#else
DISCARD_FUN(void, lunxnnf, ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int N, int L ))

size_t lunxnnf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    if (L<=0) return 0;
    return 0;
}
#endif

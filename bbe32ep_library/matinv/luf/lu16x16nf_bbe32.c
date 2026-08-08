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
#include "lunf_common.h"
#include <math.h>
#include <float.h>

#if HAVE_VFPU

#undef BBE_SELN_2XF32
#define BBE_SELN_2XF32(b,c,d) BBE_MOVN_2XF32_FROMNX16((BBE_SELNX16((BBE_MOVNX16_FROMN_2XF32(b)),(BBE_MOVNX16_FROMN_2XF32(c)),d)))


typedef void (*fnupdate)(float32_t* A, const float32_t* restrict pNorm, int i, int N, int L);

#if 0
/*--------------------------------------------------
update rows below i-th in matrix NxN (8x8, 16x16 only)

Input:
norm[L]   - normalization value for diagonal Aji, j=i+1...N-1
N           matrix size (8 or 16)
L           number of matrices
Input/output:
A[(N-i)*N]  matrix - pointer to the i-th row
Returns 
none:
--------------------------------------------------*/
void lunf_update(float32_t* A, const float32_t* restrict norm, int i, int N, int L)
{
    int j,k,l,SA;
    SA=N*N;
    NASSERT(N==8 || N==16);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    for (l=0; l<L; l++)
    {
        for (j = i + 1; j < N; j++)
        {
            A[l*SA+(j-i)*N + i] *= norm[l];
            for (k = i + 1; k < N; k++)
            {
                A[l*SA+(j-i)*N + k] -= A[l*SA+(j-i)*N + i] * A[l*SA+ k];
            }
        }
    }
}
#endif
/*--------------------------------------------------
update rows below i-th in matrix NxN (8x8, 16x16 only)

Input:
norm[L]   - normalization value for diagonal Aji, j=i+1...N-1
N           matrix size (8 or 16)
L           number of matrices
Input/output:
A[(N-i)*N]  matrix - pointer to the i-th row
Returns 
none:
--------------------------------------------------*/
// specialized version N-i==12 10...16
void lunf_update10_16(float32_t* A, const float32_t* restrict norm, int i, int N, int L)
#if 1
{
    int l;
    int SA=N*N;
    vboolN_2 bk,bi;
    xtfloat* restrict pNorm;
    xb_vecN_2xf32 * restrict pAw;
    xb_vecN_2xf32 * restrict pAj;
    xb_vecN_2xf32 * restrict pAi;
    xb_vecN_2xf32 Aj0,Aj1,Ai0,Ai1,Aji,NORM;
    vselN seli;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16);
    NASSERT(N-i>=10 && N-i<=16);
    pNorm=(xtfloat*)norm;
    bi=BBE_LTRN_2(i+1) &~ BBE_LTRN_2(i);    /* 1 in i-th position */
    bk=~BBE_LTRN_2(i+1);                    /* 1 for k>i          */
    seli=BBE_MOVVSV(BBE_REPNX16C(BBE_ADDNX16(BBE_SEQNX16(),(i<<1)),0),0); /* select for replicaiton of i-th element */
    pAw=pAj=(xb_vecN_2xf32*)&A[N ];
    pAi=(xb_vecN_2xf32 *)&A[0];
    lastinc=(SA*sizeof(float32_t) - (N*sizeof(float32_t) * (N-i-1)));
    for (l=0; l<L; l++)
    {
        int j;
        Ai1=BBE_LVN_2XF32_I (pAi,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Ai0,pAi,SA*sizeof(float32_t));
        BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
        NORM=BBE_REPN_2XF32(NORM,0);
        __Pragma("loop_count min=9")
        for (j=0; j<(N-i-1); j++)
        {
            Aj1=BBE_LVN_2XF32_I (pAj,2*BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_XP(Aj0,pAj,N*sizeof(float32_t));
            BBE_MULN_2XF32T(Aj0,Aj0,NORM,bi);
            Aji=BBE_SELN_2XF32(Aj0,Aj0,seli);
            BBE_MULSN_2XF32T(Aj0,Aji,Ai0,bk);
            BBE_MULSN_2XF32 (Aj1,Aji,Ai1);
            BBE_SVN_2XF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_XP(Aj0,pAw,N*sizeof(float32_t));
        }
        pAw+=lastinc/sizeof(xb_vecN_2xf32);
        pAj+=lastinc/sizeof(xb_vecN_2xf32);
    }
}
#else
{
    int l;
    int SA=N*N;
    vboolN_2 bk,bi;
    xtfloat* restrict pNorm;
    xb_vecN_2xf32 * restrict pAw;
    xb_vecN_2xf32 * restrict pAj;
    xb_vecN_2xf32 * restrict pAi;
    xb_vecN_2xf32 Aj00,Aj10,Ai00,Ai10,Aji0,NORM0;
    xb_vecN_2xf32 Aj01,Aj11,Ai01,Ai11,Aji1,NORM1;
    vselN seli;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16);
    NASSERT(N-i>=10 && N-i<=16);
    pNorm=(xtfloat*)norm;
    bi=BBE_LTRN_2(i+1) &~ BBE_LTRN_2(i);    /* 1 in i-th position */
    bk=~BBE_LTRN_2(i+1);                    /* 1 for k>i          */
    seli=BBE_MOVVSV(BBE_REPNX16C(BBE_ADDNX16(BBE_SEQNX16(),(i<<1)),0),0); /* select for replicaiton of i-th element */
    pAw=pAj=(xb_vecN_2xf32*)&A[N ];
    pAi=(xb_vecN_2xf32 *)&A[0];
    lastinc=(SA*sizeof(float32_t) - (N*sizeof(float32_t) * (N-i-1)));
    for (l=0; l<(L&~1); l+=2)
    {
        int j;
        Ai11=BBE_LVN_2XF32_X (pAi,SA*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        Ai01=BBE_LVN_2XF32_X (pAi,SA*sizeof(float32_t)+0*BBE_SIMD_WIDTH);

        Ai10=BBE_LVN_2XF32_I (pAi,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Ai00,pAi,2*SA*sizeof(float32_t));
        BBE_LSN_2XF32_IP(NORM0,pNorm,sizeof(float32_t));
        BBE_LSN_2XF32_IP(NORM1,pNorm,sizeof(float32_t));
        NORM0=BBE_REPN_2XF32(NORM0,0);
        NORM1=BBE_REPN_2XF32(NORM1,0);
        __Pragma("loop_count min=9")
        for (j=0; j<(N-i-1); j++)
        {
            Aj11=BBE_LVN_2XF32_X (pAj,SA*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
            Aj01=BBE_LVN_2XF32_X (pAj,SA*sizeof(float32_t)+0*BBE_SIMD_WIDTH);

            Aj10=BBE_LVN_2XF32_I (pAj,2*BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_XP(Aj00,pAj,N*sizeof(float32_t));
            BBE_MULN_2XF32T(Aj00,Aj00,NORM0,bi);
            BBE_MULN_2XF32T(Aj01,Aj01,NORM1,bi);
            Aji0=BBE_SELN_2XF32(Aj00,Aj00,seli);
            Aji1=BBE_SELN_2XF32(Aj01,Aj01,seli);
            BBE_MULSN_2XF32T(Aj00,Aji0,Ai00,bk);
            BBE_MULSN_2XF32T(Aj01,Aji1,Ai01,bk);
            BBE_MULSN_2XF32 (Aj10,Aji0,Ai10);
            BBE_MULSN_2XF32 (Aj11,Aji1,Ai11);
            BBE_SVN_2XF32_X (Aj11,pAw,SA*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_X (Aj01,pAw,SA*sizeof(float32_t)+0*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_I (Aj10,pAw,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_XP(Aj00,pAw,N*sizeof(float32_t));
        }
        pAw+=(lastinc+SA*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
        pAj+=(lastinc+SA*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
    }
    if (L&1)
    {
        int j;
        Ai10=BBE_LVN_2XF32_I (pAi,2*BBE_SIMD_WIDTH);
        Ai00=BBE_LVN_2XF32_I (pAi,0);
        BBE_LSN_2XF32_IP(NORM0,pNorm,sizeof(float32_t));
        NORM0=BBE_REPN_2XF32(NORM0,0);
        __Pragma("loop_count min=9")
        for (j=0; j<(N-i-1); j++)
        {
            Aj10=BBE_LVN_2XF32_I (pAj,2*BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_XP(Aj00,pAj,N*sizeof(float32_t));
            BBE_MULN_2XF32T(Aj00,Aj00,NORM0,bi);
            Aji0=BBE_SELN_2XF32(Aj00,Aj00,seli);
            BBE_MULSN_2XF32T(Aj00,Aji0,Ai00,bk);
            BBE_MULSN_2XF32 (Aj10,Aji0,Ai10);
            BBE_SVN_2XF32_I (Aj10,pAw,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_XP(Aj00,pAw,N*sizeof(float32_t));
        }
    }

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
void lu16x16nf ( 
            void * pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int L )
#if 0
{
    float32_t *norm=(float32_t *)pScr;
    const int N=16;
    int SA=16*16,SP=16;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    if (L<=0) return;

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
            maxVal = FLT_MIN; pivot = i;
            for (j = i; j<N; j++)
            {
                s = A[l*SA+j*N + i];
                s = fabs(s);
                if (s>maxVal) { maxVal = s; pivot = j; }
            }
            /* permute i-th and pivot-th rows */
            for (k = 0; k<N; k++)
            {
                t = A[l*SA+pivot*N + k]; A[l*SA+pivot*N + k] = A[l*SA+i*N + k]; A[l*SA+i*N + k] = t;
            }
            { int16_t t; t=P[l*SP+pivot]; P[l*SP+pivot] = P[l*SP+i]; P[l*SP+i] = t; }
        }
        for (l=0; l<L; l++)
        {
            norm[l]=1.0f/A[l*SA+i*N + i];
        }
        for (l=0; l<L; l++)
        {
            for (j = i + 1; j < N; j++)
            {
                A[l*SA+j*N + i] *= norm[l];
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
    static const fnupdate upd[]=
    {
        lunf_update16,
        lunf_update15,
        lunf_update14,
        lunf_update13,
        lunf_update12,
        lunf_update11,
        lunf_update10,
        lunf_update9,
        lunf_update8,
        lunf_update7,
        lunf_update6,
        lunf_update5,
        lunf_update4,
        lunf_update3,
        lunf_update2
    };

    float32_t *norm=(float32_t *)pScr;
    int l,i;
    xtfloat* restrict pAii;
    xtfloat* restrict pNorm;
    xb_vecNx16* restrict pP;
    xb_vecN_2xf32 Ai,NORM;
    xb_vecN_2xf32 * pAw;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    if (L<=0) return;

    pP=(xb_vecNx16*)P;
    for (l=0; l<L; l++)
    {
        BBE_SVNX16_IP(BBE_SEQNX16(),pP,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    for (i=0; i<16-1; i++)
    {
        pAii=(xtfloat*)&A[i];
        pAw=(xb_vecN_2xf32*)&A[i*16];
        if (i<8)
        {
            for (l=0; l<L; l++)
            {
                int pivot;
                xb_vecN_2xf32 Ai0,Ai1,Aj0,Aj1;
                /* find those raw permutation which maximizes A(i,i) after the row update */
                vboolN_2 b0,b1;
                vboolN b;
                vselN vb;
                int dummy;
                (void)dummy;

                xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T0,T1,maxT;
                BBE_LSN_2XF32_XP(t0,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t1,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t2,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t3,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t4,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t5,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t6,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t7,pAii,16*sizeof(float32_t));
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
                t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
                T0=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
                BBE_LSN_2XF32_XP(t0,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t1,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t2,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t3,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t4,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t5,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t6,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t7,pAii,16*sizeof(float32_t));
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
                t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
                T1=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
                T0=BBE_ABSN_2XF32(T0);
                T1=BBE_ABSN_2XF32(T1);
                BBE_CONSTN_2XF32T(T0,0,BBE_LTRN_2(i));
                maxT=BBE_RMAXNUMN_2XF32(BBE_MAXNUMN_2XF32(T0,T1));
                b0=BBE_OEQN_2XF32(maxT,T0);
                b1=BBE_OEQN_2XF32(maxT,T1);
                b=BBE_JOINBN_2(b1,b0);
                BBE_SQZN(vb,dummy,b);
                pivot=BBE_EXTRNX16(BBE_MOVVVS(vb),0);
                /* permute i-th and pivot-th rows */
                { int16_t t; t=P[l*16+pivot]; P[l*16+pivot] = P[l*16+i]; P[l*16+i] = t; }
                pivot-=i;
                Ai0=BBE_LVN_2XF32_I(pAw,0);
                Ai1=BBE_LVN_2XF32_I(pAw,2*BBE_SIMD_WIDTH);
                Aj0=BBE_LVN_2XF32_I((pAw+2*pivot),0);
                Aj1=BBE_LVN_2XF32_I((pAw+2*pivot),2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_I(Ai0,(pAw+2*pivot),0);
                BBE_SVN_2XF32_I(Ai1,(pAw+2*pivot),2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_I(Aj1,pAw,2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_XP(Aj0,pAw,256*sizeof(float32_t));
            }
        }
        else// i>=8
        {
            pAii+=16*8;
            for (l=0; l<L; l++)
            {
                int pivot;
                xb_vecN_2xf32 Ai0,Ai1,Aj0,Aj1;
                /* find those raw permutation which maximizes A(i,i) after the row update */
                vboolN_2 b0,b1;
                vboolN b;
                vselN vb;
                int dummy;
                (void)dummy;

                xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T1;
                xtfloat maxT;
                (void)maxT;
                BBE_LSN_2XF32_XP(t0,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t1,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t2,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t3,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t4,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t5,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t6,pAii,16*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t7,pAii,(16+16*8)*sizeof(float32_t));
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
                t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
                T1=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
                T1=BBE_ABSN_2XF32(T1);
                BBE_CONSTN_2XF32T(T1,0,BBE_LTRN_2(i));
                BBE_RBMAXNUMN_2XF32(b1,maxT,T1);
                b0=BBE_LTRN_2I(0);
                b=BBE_JOINBN_2(b1,b0);
                BBE_SQZN(vb,dummy,b);
                pivot=BBE_EXTRNX16(BBE_MOVVVS(vb),0);
                /* permute i-th and pivot-th rows */
                { int16_t t; t=P[l*16+pivot]; P[l*16+pivot] = P[l*16+i]; P[l*16+i] = t; }
                pivot-=i;
                Ai0=BBE_LVN_2XF32_I(pAw,0);
                Ai1=BBE_LVN_2XF32_I(pAw,2*BBE_SIMD_WIDTH);
                Aj0=BBE_LVN_2XF32_I((pAw+2*pivot),0);
                Aj1=BBE_LVN_2XF32_I((pAw+2*pivot),2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_I(Ai0,(pAw+2*pivot),0);
                BBE_SVN_2XF32_I(Ai1,(pAw+2*pivot),2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_I(Aj1,pAw,2*BBE_SIMD_WIDTH);
                BBE_SVN_2XF32_XP(Aj0,pAw,256*sizeof(float32_t));
            }
        }
        __Pragma("no_reorder")
        pNorm=(xtfloat*)pScr;
        pAii=(xtfloat*)&A[i*16+i];
        for (l=0; l<L; l++)
        {
            BBE_LSN_2XF32_XP(Ai,pAii,256*sizeof(float32_t));
            NORM=BBE_RECIPN_2XF32(Ai);
            BBE_SSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
        }
#if 0
        for (l=0; l<L; l++)
        {
            int j,k;
            for (j = i + 1; j < 16; j++)
            {
                A[l*256+j*16 + i] *= norm[l];
                for (k = i + 1; k < 16; k++)
                {
                    A[l*256+j*16 + k] -= A[l*256+j*16+ i] * A[l*256+i*16 + k];
                }
            }
        }
#endif
        __Pragma("no_reorder")
        upd[i](A+i*16,norm,i,16,L);
    }
}
#endif


size_t lu16x16nf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==16);
    return sizeof(float32_t)*L;
}
#else
DISCARD_FUN(void, lu16x16nf, ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L ))
size_t lu16x16nf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==16);
    return 0;
}

#endif

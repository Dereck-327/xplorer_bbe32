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

typedef void (*fnupdate)(float32_t* A, const float32_t* restrict norm, int i, int N, int L);

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
void lu8x8nf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L )
#if 0
{
    const int N=8;
    int SA=8*8,SP=8;
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
        for (i=0; i<N; i++) P[i]=i;
        for (i=0; i<N-1; i++)
        {
            float32_t maxVal, s, t;
            int pivot;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            maxVal = FLT_MIN; pivot = 0;
            for (j = i; j<N; j++)
            {
                s = A[j*N + i];
                s = fabs(s);
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
                float32_t norm;
                norm=1.0f/A[i*N + i];
                A[j*N + i] *= norm;
                for (k = i + 1; k < N; k++)
                {
                    A[j*N + k] -= A[j*N + i] * A[i*N + k];
                }
            }
        }
        A+=SA;
        P+=SP;
    }
}
#else
{
    static const fnupdate upd[]=
    {
        lunf_update8,
        lunf_update7,
        lunf_update6,
        lunf_update5,
        lunf_update4,
        lunf_update3,
        lunf_update2
    };

    xb_vecN_2xf32 * restrict pAi;
    xtfloat* restrict pAii;
    xtfloat* restrict pNorm;
    int l,i;
    xb_vecN_2xf32 Ai,Aj,T;
    xb_vecN_2xf32 NORM;
    valign aP;
    xb_vecNx16 * restrict pP;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    if (L<=0) return;

    /* initialize permutation vector first */
    aP=BBE_ZALIGN();
    pP=(xb_vecNx16 *)P;
    for (l=0; l<L; l++)
    {
        BBE_SAVNX16_XP(BBE_SEQNX16(),aP,pP,8*sizeof(int16_t));
    }
    BBE_SAPOS_FP(aP,pP);

    for (i=0; i<8-1; i++)
    {
        pAii=(xtfloat*)&A[i];
        pAi=(xb_vecN_2xf32 *)(A);
        __Pragma("concurrent")
        for (l=0; l<L; l++)
        {
            xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7;
            xtfloat maxT;
            vboolN_2 bmax;
            vselN vmax;
            int dummy,pivot;
            (void)dummy;
            /* find those raw permutation which maximizes A(i,i) after the row update */
            BBE_LSN_2XF32_IP(t0,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t1,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t2,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t3,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t4,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t5,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t6,pAii,8*sizeof(float32_t));
            BBE_LSN_2XF32_IP(t7,pAii,8*sizeof(float32_t));
            t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
            t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
            t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
            t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
            t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
            t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
            T=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
            /* take index of element with maximum absolute value */
            T=BBE_ABSN_2XF32(T);
            BBE_CONSTN_2XF32T(T,0,BBE_LTRN_2(i));
            BBE_RBMAXNUMN_2XF32(bmax,maxT,T);
            BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));
            pivot=BBE_EXTRNX16(BBE_MOVVVS(vmax),0);
            pivot>>=1;
            /* permute i-th and pivot-th rows */
            Ai=BBE_LVN_2XF32_X(pAi,i*8*sizeof(float32_t));
            Aj=BBE_LVN_2XF32_X(pAi,pivot*8*sizeof(float32_t));
            BBE_SVN_2XF32_X(Aj,pAi,i*8*sizeof(float32_t));
            BBE_SVN_2XF32_X(Ai,pAi,pivot*8*sizeof(float32_t));
            { int16_t t; t=P[l*8+pivot]; P[l*8+pivot] = P[l*8+i]; P[l*8+i] = t; }
            pAi+=8;
        }
        __Pragma("no_reorder")
        pNorm=(xtfloat*)pScr;
        pAii=(xtfloat*)&A[i*8+i];
        for (l=0; l<L; l++)
        {
            BBE_LSN_2XF32_XP(Ai,pAii,64*sizeof(float32_t));
            NORM=BBE_RECIPN_2XF32(Ai);
            BBE_SSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
        }
        upd[i](A+i*8,(float32_t*)pScr,i,8,L);
    }
}
#endif

//-----------------------------------------------------------------------------------
size_t lu8x8nf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==8);
    return L*sizeof(float32_t); // for 1/Aii
}
#else
DISCARD_FUN(void, lu8x8nf, ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L ))

size_t lu8x8nf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==8);
    return 0;
}
#endif

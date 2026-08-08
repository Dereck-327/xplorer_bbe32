/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Apply the QR decomposition to the matrix of normal equations system
    Make QR decomposition for block ordered matrices.
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqrn_common.h"

#if HAVE_CQRN

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-1)) m=LOG2_BBE_SIMD_WIDTH-1;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*-------------------------------------------------------------------------
Make QR decomposition for block ordered matrices.
Matrix sizes SA,SV are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SV=size(((2*M-N+1)*N/2+N)*L)
SD=size(N)
Scratch size in bytes is defined by functions cqr_build_rmxnn_getScratchSize()

Input:
 M, N         Dimensional parameters
 L            Number of matrices
Input/output:
 A[L][SA]     On input it is the sequence of L complex matrices A. 
              At the end of the process, matrices R replace input
              matrices A. In a case of non-square matrices (N!=M),
              only N*N elements of each output matrix are valid.
Output:
 V[SV]        Sequence of L Housholder rotation vectors 
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in a special format

Restrictions:
1. A, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. Number of columns for input matrices A must not exceed the number
   of rows: N <= M.
---------------------------------------------------------------------------*/
void cqr_build_rmxnn (void *pScr,
                    complex_fract16* _A,
                    complex_fract16* _V,
                    complex_fract16* _D,
                    int M, int N,
                    int L)
{
    int16_t* A=(int16_t*)_A;
    int16_t* V=(int16_t*)_V;
    int16_t* D=(int16_t*)_D;
    static const tCqrn_buildr_Iteration it_gen= {cqrnTakeColumn,cqrnHousholder,cqrnUpdateR};    /* generic code for iteration */
    static const tCqrn_buildr_Iteration it[]=
    {
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn32,cqrnHousholder32,cqrnUpdateR},
        {cqrnTakeColumn16,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn16,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn16,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn16,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn12,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn12,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn12,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn12,cqrnHousholder16,cqrnUpdateR},
        {cqrnTakeColumn8,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn8,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn8,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn8,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn4,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn4,cqrnHousholder8  ,cqrnUpdateR},
        {cqrnTakeColumn4,cqrnHousholder8  ,cqrnUpdateR},
    };
    int idx;
    int m,Ncolumns=N==M?N-1:N;
    int SA=2*getSpace(M*N);
    int SD=2*getSpace(N);
    int16_t* Z=(int16_t*)pScr;
    int16_t* Fi;
    int16_t* pV;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D   ,2*BBE_SIMD_WIDTH);
    NASSERT(M%4==0 && M>0);
    NASSERT(N%4==0 && N>0);
    NASSERT(L>0);
    {
        size_t Zsize;
        Zsize=2*M*L*sizeof(int16_t);
        Zsize=(Zsize+2*BBE_SIMD_WIDTH-1) & ~(2*BBE_SIMD_WIDTH-1);
        pScr=(void*)(((uintptr_t)pScr)+Zsize);
    }
#if EXTRA_SHIFT
    {
        int32_t rnd=1<<(EXTRA_SHIFT-1);
        for (m=0;m<SA*L;m++) 
        {
            int32_t t;
            t=A[m]+rnd;
            t>>=EXTRA_SHIFT;
            A[m]=(int16_t)t;
        }
    }
#endif
    Fi=V+(2*M-N+1)*N*L;
    pV=V;
    idx = 32-M;
    for(m=0; m<Ncolumns; m++,idx++)
    {
        const tCqrn_buildr_Iteration *pItr;
        pItr = idx<0 ? &it_gen : it+idx;
        pItr->takeColumn (Z,A+2*m*(N+1),M-m,N,SA,L);
        pItr->housholder (pScr,pV,Fi+2*m*L,D+2*m,Z,M-m,SD,L);
        pItr->updateR (Z,A+2*m*(N+1),pV,SA,M-m,N-m,N,L);
        pV+=2*(M-m)*L;
    }
    if (M==N)
    {
        cqrnTakeColumn1(Z,A+2*m*(N+1),SA,L);
        cqrnHousholder1(pV,Fi+2*m*L,D+2*m,Z,SD,L);
        cqrnUpdateR1   (A+2*m*(N+1),pV,SA,L);
    }

    // final rotation Fi'*R
    cqrnRotateR(A,Fi,N,SA,L,0);
#if EXTRA_SHIFT
   for (m=0;m<(L*SD)/2;m++) D[2*m+1]+=EXTRA_SHIFT;
#endif
} /* cqr_build_rmxnn() */

size_t cqr_build_rmxnn_getScratchSize (int M, int N,int L)
{
    size_t Zsize,Hsize;
    (void)M;
    (void)N;
    (void)L;
    Zsize=2*M*L*sizeof(int16_t);
    Zsize=(Zsize+2*BBE_SIMD_WIDTH-1) & ~(2*BBE_SIMD_WIDTH-1);
    Hsize=cqrnHousholder_getScratchSz(M,L);
    return Zsize+Hsize;
} /* cqr_build_rmxnn_getScratchSize() */
#else
DISCARD_FUN(void,cqr_build_rmxnn,(void *pScr,
                    complex_fract16* A,
                    complex_fract16* V,
                    complex_fract16* D,
                    int M, int N,
                    int L))

size_t cqr_build_rmxnn_getScratchSize(int M, int N,int L) { (void)M,(void)N,(void)L; return 0; }
#endif

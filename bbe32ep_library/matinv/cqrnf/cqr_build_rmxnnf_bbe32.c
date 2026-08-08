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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, block format
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "cqrnf_common.h"
#include "common.h"

#if HAVE_VFPU

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    m=XT_MIN(m ,LOG2_BBE_SIMD_WIDTH-1);
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*-------------------------------------------------------------------------
Make QR decomposition for block ordered matrices.
Matrix sizes SA,SV are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SV=size(((2*M-N+1)*N/2+N)*L)
SD=size(N)
Scratch size in bytes is defined by functions xxxx_getScratchSize()

Input:
 M, N         Dimensional parameters
 L            Number of matrices
Input/output:
 A[L][SA]     On input it is the sequence of L matrices A. 
              At the end of the process, matrices R replace input
              matrices A. In a case of non-square matrices (N!=M),
              only N*N elements of each output matrix are valid.
Output:
 V[SV]        Sequence of L Housholder rotation vectors 
 D[L*SD]      Reciprocals of main diagonal in a special format

Restrictions:
1. A, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. Number of columns for input matrices A must not exceed the number
   of rows: N <= M.
---------------------------------------------------------------------------*/
void  cqr_build_rmxnnf  (void *pScr,complex_float* _A,complex_float* _V,complex_float* _D,int M, int N,int L)
{
    static const tCqrnIteration it_gen= {cqrnfTakeColumn,cqrnfHousholder,cqrnfUpdateR};    /* generic code for iteration */
    static const tCqrnIteration it[]=
    {
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR}
    };
    typedef void (*fnupdateR)(float32_t*,float32_t*,const float32_t*,int,int,int,int,int);
    static const fnupdateR updTbl[]=
    {
        cqrnfUpdateR4 ,cqrnfUpdateR4 ,cqrnfUpdateR4 ,cqrnfUpdateR4 , cqrnfUpdateR4 ,
        cqrnfUpdateR8 ,cqrnfUpdateR8 ,cqrnfUpdateR8 ,cqrnfUpdateR8 ,
        cqrnfUpdateR12,cqrnfUpdateR12,cqrnfUpdateR12,cqrnfUpdateR12,
        cqrnfUpdateR16,cqrnfUpdateR16,cqrnfUpdateR16,cqrnfUpdateR16
    };
    size_t Hsize;
    float32_t* A=(float32_t*)_A;
    float32_t* V=(float32_t*)_V;
    float32_t* D=(float32_t*)_D;
    int m,Ncolumns=N==M?N-1:N;
    int SA=getSpace(M*N<<1);
    int SD=getSpace(N  <<1);
    int idx;
    float32_t* Z;
    float32_t* Zt;
    float32_t* Fi;
    float32_t* pV;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT(N<=M);
    NASSERT(N%4==0 && M%4==0);
    NASSERT(N>0 && M>0);
    NASSERT(L>0);
    if (L<=0 || M<=0 || N<=0 || M<N) return;
    Hsize=cqrnfHousholder_getScratchSize(M,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Z=(float32_t*)(((uintptr_t)pScr)+Hsize);
    Zt=Z+2*M*L;
    Zt=(float32_t*)((((uintptr_t)Zt)+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1));
    Fi=V+(2*M-N+1)*N*L;
    pV=V;

    idx = 16-M;
    for(m=0; m<Ncolumns; m++,idx++)
    {
        const tCqrnIteration *pItr;
        fnupdateR upd;
        pItr = idx<0 ? &it_gen : it+idx;
        upd = N-m >= (int)(sizeof(updTbl)/sizeof(updTbl[0])) ? cqrnfUpdateR : updTbl[N-m];
        pItr->takeColumn (Z,Zt,A+2*m*(N+1),M-m,N,SA,L);
        pItr->housholder (pScr,pV,Fi+2*m*L,D+2*m,Z,Zt,M-m,SD,L);
        upd (Z,A+2*m*(N+1),pV,SA,M-m,N-m,N,L);
        pV+=2*(M-m)*L;
    }
    if (M==N)
    {
        cqrnfTakeColumn1(Z,A+2*m*(N+1),SA,L);
        cqrnfHousholder1(pV,Fi+2*m*L,D+2*m,Z,SD,L);
        cqrnfUpdateR1   (A+2*m*(N+1),pV,SA,L);
    }
    // final rotation Fi'*R
    cqrnfrotateR(A,Fi,N,SA,L);
}

/* scratch allocation functions */
size_t  cqr_build_rmxnnf_getScratchSize  (int M, int N,int L)
{
    size_t Zsize,Hsize;
    NASSERT(N<=M);
    NASSERT(N%4==0 && M%4==0);
    NASSERT(N>0 && M>0);
    NASSERT(L>0);
    L=XT_MAX(L,0);
    M=XT_MAX(M,0);
    N=XT_MAX(N,0);
    Hsize=cqrnfHousholder_getScratchSize(M,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Zsize=M*L*sizeof(complex_float);
    Zsize=(Zsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    return 2*Zsize+Hsize;
}
#else
DISCARD_FUN(void, cqr_build_rmxnnf, (void *pScr,complex_float* _A,complex_float* _V,complex_float* _D,int M, int N,int L))
size_t  cqr_build_rmxnnf_getScratchSize(int M, int N, int L)
{
    NASSERT(N<=M);
    NASSERT(N%4==0 && M%4==0);
    NASSERT(N>0 && M>0);
    NASSERT(L>0);
    L=XT_MAX(L,0);
    M = XT_MAX(M, 0);
    N = XT_MAX(N, 0);
    return 0;
}
#endif

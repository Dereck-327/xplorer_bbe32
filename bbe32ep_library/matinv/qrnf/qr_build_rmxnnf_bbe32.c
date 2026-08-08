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
  QR decomposition, floating point, real data, block format
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "qrnf_common.h"

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

/*
Reference Matlab code:
% build R matrix, Housholder vectors v and diagonal rotation matrix FI
% output:
% V  - sequence of Housholder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% D  - reciprocals of main diagonal elements
% R  - upper triangle decomposition
function [V,Fi,D,Q,R] = cqr_buildr(A)
A0=A;
[M, N] = size(A); 
if (M==N) Ncolumns=N;
else      Ncolumns=N;
end

V =[];
Q = eye(M); 
Fi = zeros(N,1);
D  = zeros(N,1);

for m=1: Ncolumns
    tmpA = A(m:end, m:end);
    e1 = zeros(M-m+1,1);
    e1(1) = 1;
    x = tmpA(:,1);

    fi = x(1)/(1E-10+abs(x(1)));
    D(m)=  1/sqrt(x'*x);
    alpha = - 1/D(m) * fi;
    v = x - alpha*e1;
    Fi(m) = fi;

    v = v/sqrt(v'*v);
    V=[V;v];
    P = eye(M-m+1) - 2*v*v';
    A(m:end, m:end) = P*tmpA;
    A(m,m)=alpha;
end
Fi=-Fi; % make diagonals positive 
F = diag([Fi;ones(M-N,1)]);
Q = Q*F;
R = F'*A;
*/

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
void  qr_build_rmxnnf  (void *pScr,float32_t* A,float32_t* V,float32_t* D,int M, int N,int L)
{
    static const tQrnIteration it_gen= {qrnfTakeColumn,qrnfHousholder,qrnfUpdateR};    /* generic code for iteration */
    static const tQrnIteration it[]=
    {
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR16},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR16},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR16},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR16},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR12},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR12},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR12},
        {qrnfTakeColumn16, qrnfHousholder16, qrnfUpdateR12},
        {qrnfTakeColumn8 , qrnfHousholder8 , qrnfUpdateR8 },
        {qrnfTakeColumn8 , qrnfHousholder8 , qrnfUpdateR8 },
        {qrnfTakeColumn8 , qrnfHousholder8 , qrnfUpdateR8 },
        {qrnfTakeColumn8 , qrnfHousholder8 , qrnfUpdateR8 },
        {qrnfTakeColumn8 , qrnfHousholder4 , qrnfUpdateR4 },
        {qrnfTakeColumn8 , qrnfHousholder4 , qrnfUpdateR4 },
        {qrnfTakeColumn8 , qrnfHousholder4 , qrnfUpdateR4 }
    };
    typedef void (*fnupdateR)(float32_t*,float32_t*,const float32_t*,int,int,int,int,int);
    static const fnupdateR updTbl[]=
    {
        qrnfUpdateR4 ,qrnfUpdateR4 ,qrnfUpdateR4 ,qrnfUpdateR4 , qrnfUpdateR4 ,
        qrnfUpdateR8 ,qrnfUpdateR8 ,qrnfUpdateR8 ,qrnfUpdateR8 ,
        qrnfUpdateR12,qrnfUpdateR12,qrnfUpdateR12,qrnfUpdateR12,
        qrnfUpdateR16,qrnfUpdateR16,qrnfUpdateR16,qrnfUpdateR16
    };

    int m,idx,Ncolumns=N==M?N-1:N;
    int SA=getSpace(M*N);
    int SD=getSpace(N);
    size_t Hsize;
    float32_t* Z,*Zt;
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
    Hsize=qrnfHousholder_getScratchSize(M,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Z=(float32_t*)(((uintptr_t)pScr)+Hsize);
    Zt=Z+M*L;
    Zt=(float32_t*)((((uintptr_t)Zt)+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1));

    Fi=V+(((2*M-N+1)*N)>>1)*L;
    pV=V;
    idx = 16-M;
    for(m=0; m<Ncolumns; m++,idx++)
    {
        const tQrnIteration *pItr;
        fnupdateR upd;
        pItr = idx<0 ? &it_gen : it+idx;
        upd = N-m >= (int)(sizeof(updTbl)/sizeof(updTbl[0])) ? qrnfUpdateR : updTbl[N-m];
        pItr->takeColumn (Z,Zt,A+m*(N+1),M-m,N,SA,L);
        pItr->housholder (pScr,pV,Fi+m*L,D+m,Z,Zt,M-m,SD,L);
        upd (Z,A+m*(N+1),pV,SA,M-m,N-m,N,L);
        pV+=(M-m)*L;
    }
    if (M==N)
    {
        qrnfTakeColumn1(Z,A+m*(N+1),SA,L);
        qrnfHousholder1(pV,Fi+m*L,D+m,Z,SD,L);
        qrnfUpdateR1   (A+m*(N+1),pV,SA,L);
    }

    // final rotation Fi'*R
    qrnfRotateR(A,Fi,N,SA,L);
}

/* scratch allocation functions */
size_t  qr_build_rmxnnf_getScratchSize  (int M, int N,int L)
{
    size_t Hsize,Zsize;
    NASSERT(N<=M);
    NASSERT(N%4==0 && M%4==0);
    NASSERT(N>0 && M>0);
    NASSERT(L>0);
    L=XT_MAX(L,0);
    M=XT_MAX(M,0);
    N=XT_MAX(N,0);
    Hsize=qrnfHousholder_getScratchSize(M,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Zsize=M*L*sizeof(float32_t);
    Zsize=(Zsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    return 2*Zsize+Hsize;
}

#else
DISCARD_FUN(void, qr_build_rmxnnf, (void *pScr,float32_t* A,float32_t* V,float32_t* D,int M, int N,int L))
size_t  qr_build_rmxnnf_getScratchSize  (int M, int N,int L)
{
    NASSERT(N<=M);
    NASSERT(N%4==0 && M%4==0);
    NASSERT(N>0 && M>0);
    NASSERT(L>0);
    return 0;
}
#endif

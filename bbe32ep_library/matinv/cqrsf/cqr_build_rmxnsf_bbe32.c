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
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

/*
Reference Matlab code:
% build R matrix, Householder vectors v and diagonal rotation matrix FI
% output:
% V  - sequence of Householder vectors  [(2*M-N+1)*N/2,1]
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
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cqrsf_common.h"
#if HAVE_VFPU

/*-----------------------------------------------------------------------
[c]qr_build_rMxNsf

QR decomposition of MxN real/complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

Data format: IEEE-754 Std single precision floating-point

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]               matrices A (L matrices of size MxN)
Output:
V[((2*M-N+1)*N/2+N)*L]  L sets of Householder vectors
R[M*N][L]               upper triangular matrices (L matrices
                        of size MxN)
D[N*L]                  reciprocals of main diagonal

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
void cqr_build_rmxnsf   ( void* pScr, complex_float* restrict _V, complex_float * restrict _R, complex_float * restrict _D, int M, int N,int L)
{
     typedef void (*fnhousholder)(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
     typedef void (*fnupdateR)(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
     static const fnhousholder htbl[]={cqrfsHousholder,cqrfsHousholder2,cqrfsHousholder3,cqrfsHousholder4};
     static const fnupdateR    utbl[]={cqrfsUpdateR1,cqrfsUpdateR2,cqrfsUpdateR3,cqrfsUpdateR4};

    float32_t *V=(float32_t *)_V;
    float32_t *R=(float32_t *)_R;
    float32_t *D=(float32_t *)_D;

    int m,Ncolumns=N==M?N-1:N;
    int SV=(2*M-N+1)*N; // number of elements in Householder without rotator
    float32_t* Fi;
    float32_t* pV;
    float32_t* pA;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(V   ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D   ,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(M>=N);

    Fi=V+SV*L;
    pV=V;
    pA=R; // pointer to the diagonal element
    for(m=0; m<Ncolumns; m++)
    {
        fnhousholder housholder=M-m > 4 ? cqrfsHousholder:htbl[M-m-1];
        fnupdateR    updateR   =M-m > 4 ? cqrfsUpdateR:utbl[M-m-1];
        housholder (pScr,pV,Fi,D,pA,SV,M-m,N,L);
        updateR    (pA,pV,m,M-m,N,L);
        pA+=2*(N+1)*L; 
        pV+=2*(M-m)*L;
        D+=2*L;
        Fi+=2*L;
    }
    if (M==N)
    {
        cqrfsHousholder1(pV,Fi,D,pA,SV,M-m,N,L);
        cqrfsUpdateR1(pA,pV,m,1,N,L);
    }

    Fi=V+SV*L;
    // final rotation Fi'*R
    cqrfsRotateR(R,Fi,M,N,L);
}

// scratch memory needed for cqrxxxs functions
size_t cqr_build_rmxnsf_getScratchSize  (int M, int N,int L) 
{ 
    (void)M,(void)N,(void)L;
    NASSERT(M>=N);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return cqrfsHousholder_getScratchSize(M,N,L);
}
#else
DISCARD_FUN(void, cqr_build_rmxnsf, ( void* pScr, complex_float* restrict _V, complex_float * restrict _R, complex_float * restrict _D, int M, int N,int L))
// scratch memory needed for cqrxxxs functions
size_t cqr_build_rmxnsf_getScratchSize  (int M, int N,int L) 
{ 
    (void)M,(void)N,(void)L;
    NASSERT(M>=N);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return 0;
}
#endif

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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*===========================================================================
  Cholesky MMSE solution API (real floating point data)
  C code optimized for BBE32

  Integrit 2006-2017
===========================================================================*/
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

/*-------------------------------------------------------------------------
Compute the MMSE solution for a system of linear equations A*x=b, where A
is an MxN real (complex) matrix with M>=N and rank(A)==N, x is an Nx1 vector
of unknowns, and b is an Mx1 right hand side vector. This task is accomplished
in 3 steps:
  1. Cholesky decomposition is applied to the matrix of normal equations
     system, which results in an upper triangular NxN matrix R with real
     and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
     where adj(...) denotes the (conjugate) transpose of a matrix, and
     sigma2*I is an NxN identity matrix multiplied by the regularization
     term.
  2. Forward recursion step: the system adj(R)*y=adj(A)*b is solved for an 
     Nx1 vector y.
  3. Backward recursion step: the system R*x=y is solved for the Nx1 vector
     of unknows x.

For a single MxN matrix A, these 3 steps may be done simultaneously for P
variants of Mx1 right hand side column vectors b gathered into an MxP input
matrix B. MMSE solution is computed independently for each of P columns,
with resulting column vectors forming the solution matrix X of size NxP.

[r]cholmmse<size>sf() functions process L pairs of MxN matrices A and MxP 
matrices B in a single call, which results in L solution matrices X of
size NxP. L matrices of each kind are stored as stream ordered sequences.

Data format: IEEE-754 Std. single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [r]cholmmse<size>sf_getScratchSize(M,N,P,L)
Input:
  M,N,P      Dimensional parameters
  L          Number of matrices
  sigma2[L]  Regularization term
  A[M*N][L]  Sequence of L matrices A
  B[M*P][L]  Sequence of L right hand side matrices B
Output:
  x[N*P][L]  Sequence of L solution matrices X
Restrictions:
  1. pScr,x,A,B,sigma2 must not overlap
  2. pScr,x,A,B,sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 4 for complex-valued 
     functions, or a multiple of 8 for real-valued functions.
  4. Matrix sizes M,N must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
#if !(HAVE_VFPU)
DISCARD_FUN(void, rcholmmsemxnxpsf,(
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int M, int N, int P, int L ))
size_t rcholmmsemxnxpsf_getScratchSize ( int M, int N, int P, int L )
{
    (void)M, (void)N, (void)P, (void)L;
    NASSERT(M>1);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    return 0;
}
#else

void rcholmmsemxnxpsf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int M, int N, int P, int L )
{
    float32_t *y;
    float32_t *R;
    float32_t *D;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2,(2*BBE_SIMD_WIDTH));
    NASSERT(M>1);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    /* allocate temporary data on scratch */
    R=(float32_t*)pScr;
    y=(float32_t*)(R+N*N*L);
    D=(float32_t*)(y+M*P*L);
    rcholmxnsf(R,D,A,sigma2,M,N,L);
    rcholfwdmxnxpsf(y,R,D,A,B,M,N,P,L);
    rcholbkwnxpsf(x,R,D,y,N,P,L);
}

size_t rcholmmsemxnxpsf_getScratchSize ( int M, int N, int P, int L )
{
    size_t sz;
    NASSERT(M>1);
    NASSERT(N>=1);
    NASSERT(P>0);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    sz=sizeof(float32_t)*(N*N*L+M*P*L+N*L);
    return sz;
}
#endif

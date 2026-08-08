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
  Cholesky MMSE solution API (complex floating point data)
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
DISCARD_FUN(void, cholmmse2x2x1sf,(
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L ))
size_t cholmmse2x2x1sf_getScratchSize ( int M, int N, int P, int L )
{
    (void)M, (void)N, (void)P, (void)L;
    NASSERT(M==2 && N==2 && P==1);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/4-1))==0);
    return 0;
}
#else

void cholmmse2x2x1sf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L )
{
    complex_float *y;
    complex_float *R;
    complex_float *D;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2,(2*BBE_SIMD_WIDTH));
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/4-1))==0);
    /* allocate temporary data on scratch */
    R=(complex_float*)pScr;
    y=(complex_float*)(R+2*2*L);
    D=(complex_float*)(y+2*1*L);
    chol2x2sf(R,D,A,sigma2,L);
    cholfwd2x2x1sf(y,R,D,A,B,L);
    cholbkw2x1sf(x,R,D,y,L);
}

/* Return the scratch area size, in bytes. */
size_t cholmmse2x2x1sf_getScratchSize ( int M, int N, int P, int L )
{
    size_t sz;
    (void)M, (void)N, (void)P, (void)L;
    NASSERT(M==2 && N==2 && P==1);
    NASSERT(L>1);
    NASSERT((L&(BBE_SIMD_WIDTH/4-1))==0);
    sz=sizeof(complex_float)*(2*2*L+2*1*L+2*L);
    return sz;
}
#endif

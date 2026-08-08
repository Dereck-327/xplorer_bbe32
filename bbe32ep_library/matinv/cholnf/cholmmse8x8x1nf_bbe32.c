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
/*
    Cholesky decomposition, floating point complex data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#if (HAVE_VFPU)

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

Storage sizes SA,SB,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SB = storage_size(M*P)
SX = storage_size(N*P)

[r]cholmmse<size>nf() functions process L pairs of MxN matrices A and MxP 
matrices B in a single call, which results in L solution matrices X of
size NxP. L matrices of each kind are stored as block ordered sequences.

Data format: IEEE-754 Std. single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [r]cholmmse<size>sf_getScratchSize(M,N,P,L)
Input:
  M,N,P      Dimensional parameters
  L          Number of matrices
  sigma2[L]  Regularization term
  A[L][SA]  Sequence of L matrices A
  B[L][SB]  Sequence of L right hand side matrices B
Output:
  x[L][SZ]  Sequence of L solution matrices X
Restrictions:
  1. pScr,x,A,B,sigma2 must not overlap
  2. pScr,x,A,B,sigma2 must be aligned on 32-byte boundary
  3. M and N should be a multiple of 4
  4. Matrix sizes M,N must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
void cholmmse8x8x1nf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L )
{
    xb_vecN_2xf32 * pZ;
    int n, N_;
    complex_float * restrict R;
    complex_float * restrict Y;
    complex_float * restrict D;
    int SR,SD,SY;
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2,2*BBE_SIMD_WIDTH);
    if (L<=0) return;
    SR=L*(((8*(8+1)))>>1)*sizeof(complex_float);
    SD=L*((8  *2)>>1    )*sizeof(complex_float);
    SY=L*((8*1*2)>>1    )*sizeof(complex_float);
    SR=(SR+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    SD=(SD+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    SY=(SY+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    pZ = (xb_vecN_2xf32 *)pScr;
    N_ = cholmmse8x8x1nf_getScratchSize(8, 8, 1, L) >> (LOG2_BBE_SIMD_WIDTH + 1);
    for (n = 0; n < N_; n++)
    {
        BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(), pZ, 2 * BBE_SIMD_WIDTH);
    }
    R   =(complex_float *)pScr;
    D   =(complex_float *)(((uintptr_t)R)+SR);
    Y   =(complex_float *)(((uintptr_t)D)+SD);
    pScr=(complex_float *)(((uintptr_t)Y)+SY);
    chol8x8nf(pScr,R,D,A,sigma2,L);
    cholfwd8x8x1nf(pScr,Y,R,D,A,B,L);
    cholbkw8x1nf(pScr,x,R,D,Y,L);
}


/* Return the scratch area size, in bytes. */
size_t cholmmse8x8x1nf_getScratchSize   ( int M, int N, int P, int L )
{
    NASSERT(M==8 && N==8 && P==1);
    size_t szLoc,szScr;
    int SR,SD,SY;
    NASSERT(M>1 && N>1 && N<=N && M%4==0 && N%4==0 && P>0);
    if (L<=0) return 0;
    SR=L*(((N*(N+1)))>>1)*sizeof(complex_float);
    SD=L*((N  *2)>>1    )*sizeof(complex_float);
    SY=L*((N*P*2)>>1    )*sizeof(complex_float);
    SR=(SR+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    SD=(SD+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    SY=(SY+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    szLoc=(SR+SD+SY);
    szScr=cholmxnnf_getScratchSize(M,N,L);
    szScr=XT_MAX(szScr,cholfwdmxnxpnf_getScratchSize(M,N,P,L));
    szScr=XT_MAX(szScr,cholbkwnxpnf_getScratchSize(N,P,L));
    return szLoc+szScr;
}
#else
DISCARD_FUN(void, cholmmse8x8x1nf,(
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L ))

size_t cholmmse8x8x1nf_getScratchSize   ( int M, int N, int P, int L )
{
  (void)M; (void)N; (void)P; (void)L;
  return 0;
}
#endif

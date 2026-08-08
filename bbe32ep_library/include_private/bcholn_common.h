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
/*===========================================================================
  Banded Cholesky backward recursion for pseudo-inversion API (complex data)
  bcholbkwwxnxpn    matrices (W+N-1)xNxP (band width W), P==1

  Integrit 2006-2016
===========================================================================*/
#ifndef BCHOLN_COMMON_H__
#define BCHOLN_COMMON_H__
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_types.h"
#include "common.h"

// bcholn_xxx require BBE_RSQRTLUNX40_0 and VSA arithmetic (BBE_SUBSR1SAVSN, etc.)
#define HAVE_BCHOLN (HAVE_RSQRT && HAVE_VSAMATH)

#if (HAVE_BCHOLN)
/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion for 
specific band width W. They use Cholesky decomposition
of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SR,SD,SY,SX are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SR=size(W*N)
SD=size(N)
SY=size(N*P)
SX=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR][2]  Cholesky upper-triangle matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagonal elements are in the last raw)
D[L][SD][2]   Sequence of L reciprocals of main diagonal R represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
Yt[L][SY][2]  Results of forward recursion stage. SY is computed as for complex 
              block ordered matrices of size N*P
qA,qX,qY      Fixed point representation of matrices A(or R which is the same), 
              x and y
Output:
Xt[L][SX][2]  Decision matrix x. SX is computed as for complex block ordered 
              matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
// W==4
int bcholbkw4xnx1n  (
            int16_t* restrict Xt, 
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict Yt, 
            int qXYA, int N, int L);

// W==8
int bcholbkw8xnx1n  (
            int16_t* restrict Xt, 
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict Yt, 
            int qXYA, int N, int L);
// W==12
int bcholbkw12xnx1n  (
            int16_t* restrict Xt, 
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict Yt, 
            int qXYA, int N, int L);

// W==16
int bcholbkw16xnx1n  (
            int16_t* restrict Xt, 
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict Yt, 
            int qXYA, int N, int L);

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion for 
specific band width W. They use Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD,SB,SY are selected as 
usual for complex block ordered matrix sequencies, i.e. total size is 
rounded up to the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)
SB=size((W+N-1)*P)
SY=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR][2]  Cholesky upper-triangle matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagonal elements are in the last raw)
D[L][SD][2]   Sequence of L reciprocals of main diagonal A represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
At[L][SA][2]  Original left-side matrices A represented in the compact 
              form (only band)
Bt[L][SB][2]  Original right-side matrices B. SB is computed as for complex 
              block ordered matrices of size (W+N-1)*P
qA,qB,qY      Fixed point representation of matrices A (or R which is the 
              same),B and y

Output:
Yt[L][SY][2]  Decision matrix y. SY is computed as for complex 
              block ordered matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
// W==4
int bcholfwd4xnx1n (
            int16_t* restrict Yt,
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict At, 
            const int16_t* restrict Bt, 
            int qYB,
            int N, int L);
// W==8
int bcholfwd8xnx1n (
            int16_t* restrict Yt,
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict At, 
            const int16_t* restrict Bt, 
            int qYB,
            int N, int L);
// W==12
int bcholfwd12xnx1n (
            int16_t* restrict Yt,
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict At, 
            const int16_t* restrict Bt, 
            int qYB,
            int N, int L);
// W==16
int bcholfwd16xnx1n (
            int16_t* restrict Yt,
            const int16_t* restrict Rt, 
            const int16_t* restrict D, 
            const int16_t* restrict At, 
            const int16_t* restrict Bt, 
            int qYB,
            int N, int L);

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right-hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxpn() and cholbkwnxpn(), respectively.

The code for banded matrices is intended for cases where matrix A contains 
W non-zero elements on the main diagonal and below. So, size M is N+W-1. 
Matrix A may be stored in the compact form of size WxN. Cholesky matrix R 
also has WxN non-zero elements

NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)


Input:
  W             Band width
  N             Dimensional parameters
  L             Number of matrices
  sigma2[L]     Regularization term; fixed point position is twice the
                number of fractional bits for matrices A, R
  At[L][SA][2]  Sequence of L complex matrices A represented in the 
                compact form (only band)
Output:
  Rt[L][SR][2]  Sequence of L upper triangular complex matrices R 
                represented in the compact form (saved only elements on 
                the main diagonal and above in such a way that diagonal
                elements are in the last raw)
  D[L][SD][2]   Sequence of L reciprocals of main diagonal A represented 
                in the  block floating point (mantissa and exponent).

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
Note:
  Function may speculatively read up to (W-2)*W complex elements
  beyond the upper bound of At and Rt.
---------------------------------------------------------------------------*/
// W==4
int bchol4xnn(
            int16_t * restrict Rt, 
            int16_t * restrict D, 
            const int16_t * restrict At, 
            const int32_t * restrict sigma2,
            int N,int L);
// W==8
int bchol8xnn(
            int16_t * restrict Rt, 
            int16_t * restrict D, 
            const int16_t * restrict At, 
            const int32_t * restrict sigma2,
            int N,int L);

#endif  // HAVE_BCHOLN

#endif

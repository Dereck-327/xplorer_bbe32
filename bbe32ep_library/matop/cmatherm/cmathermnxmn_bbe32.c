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
  NatureDSP_Baseband library. Matrix Operations
    Matrix Hermitian Product
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "cmathermnxmn_common.h"
#if (HAVE_MULPC && HAVE_PACKEDMUL && 1)
#define MAX(x,y) ((x)>(y)?(x):(y))
/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     M,N must be multiples of 4
*/
void cmathermnxmn ( void * pScr,
                    complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int M, int L, int Q )
{
  int _L;

  int MN = M*N, NN = N*N;

  NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);

  NASSERT(!(M & 3) && !(N & 3));
  if (L <= 0 || N <= 0 || M <= 0) return;
  if ((_L = (L&~7)))
  {
    cmathermnxmn_L8(pScr, y, x, N, M, _L, Q);
  }
  __Pragma("no_reorder")
  if ((_L = (L & 6)))
  {
    cmathermnxmn_L2(pScr,
      y + (L&~7)*NN,
      x + (L&~7)*MN,
      N, M, _L, Q);
  }
  __Pragma("no_reorder")
  if ((_L = (L & 1)))
  {
    cmathermnxmn_tail(pScr,
      y + (L&~1)*NN,
      x + (L&~1)*MN,
      N, M, Q);
  }
}

/* Return the scratch area size, in bytes. */
size_t cmathermnxmn_getScratchSize ( int N, int M, int L )
{
  size_t L2_scratch, L8_scratch, sz = 0;
  if (M <= 0 || N <= 0) return 0;
  (void)L;
  if (L >= 8) { L8_scratch = cmathermnxmn_L8_getScratchSize(N, M); sz = MAX(sz, L8_scratch); }

  if (L&6) { L2_scratch = cmathermnxmn_L2_getScratchSize(N, M); sz = MAX(sz, L2_scratch); }
  if (L & 1){ L2_scratch = cmathermnxmn_tail_getScratchSize(N, M); sz = MAX(sz, L2_scratch); }
  return sz;
} /* cmathermnxmn_getScratchSize() */
#else
DISCARD_FUN(void, cmathermnxmn, (void * pScr,
                   complex_fract16 * restrict y, 
             const complex_fract16 * restrict x, 
             int N, int M, int L, int Q ))
size_t cmathermnxmn_getScratchSize(int N, int M, int L) { (void)N; (void)M; (void) L; return 0; }
#endif

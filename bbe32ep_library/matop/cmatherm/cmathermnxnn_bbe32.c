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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmathermnxnn, (void * pScr,
                   complex_fract16 * restrict y, 
             const complex_fract16 * restrict x, 
             int N,int L, int Q ))
size_t cmathermnxnn_getScratchSize(int N, int L) { (void)N; (void) L; return 0; }
#else


/* Block Order, NxN*NxN->NxN, Sx=NxN, Sy=NxN
   Restrictions:
     N must be a multiple of 4
*/
void cmathermnxnn ( void * pScr,
                    complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int L, int Q )
{
  NASSERT_ALIGN(pScr, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);

  NASSERT(!(N & 3));

  cmathermnxmn(pScr, y, x, N, N, L, Q);

} /* cmathermnxnn() */

/* Return the scratch area size, in bytes. */
size_t cmathermnxnn_getScratchSize ( int N, int L )
{
  return cmathermnxmn_getScratchSize(N, N, L);
} /* cmathermnxnn_getScratchSize() */
#endif

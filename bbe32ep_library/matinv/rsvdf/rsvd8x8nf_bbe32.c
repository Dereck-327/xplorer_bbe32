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
/*
    NatureDSP_Baseband library. Singular Value Decomposition
    Real 8x8 block ordered matrices
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
/* SVD common declarations */
#include "svd_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
Thin SVD For Real/Complex Block Ordered Matrices

Description: compute the Thin Singular Value Decomposition of L complex 
(real) MxN matrices, with the number of rows greater than or equal to the
number of columns: M>=N. Input and output matrices are stored in block order.

Data format: IEEE-754 Std single precision floating-point

Storage sizes SA,Ss,SV denote the number of data elements required to store
a matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
Ss = storage_size(N)
SV = storage_size(N*N)

Notes:
1. SVD implementation may perform in-place transformations of input matrices,
   so INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once U or V matrix is not required, set the corresponding output pointer
   parameter to zero to allow for a lower complexity implementation of the
   SVD algorithm.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.

Temporary:
  pScr          Scratch area. Required size (in bytes) is defined by 
                functions [r]svd<size>nf_getScratchSize(N,L)
Input:
  M,N           Matrix dimensions
  L             Number of matrices
  A[L][SA]      MxN input matrices
Output:
  U[L][SA]      MxN matrices comprised of M left-singular column 
                vectors (optional)
  s[L][Ss]      Nx1 vectors of singular values in descending order. In an 
                exceptional case when the iterative algorithm fails to 
                converge for a particular matrix, all elements of the
                respective vector are set to NaN.
  V[L][SV]      NxN matrices comprised of N right-singular column
                vectors (optional)
Restrictions:
  pScr,U,s,V,A  Must not overlap and must be aligned on 32-byte boundary 
  M,N           Must be positive multiples of 4 such that M>=N
---------------------------------------------------------------------------*/
void rsvd8x8nf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L )
{
  NASSERT_ALIGN( U, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( s, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );

  rsvdmxnnf(pScr,U,s,V,A,8,8,L);

} /* rsvd8x8nf() */

size_t rsvd8x8nf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 8==M && 8==N );
  return ( rsvdmxnnf_getScratchSize(M,N,L) );
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, rsvd8x8nf, ( void * pScr,
                           float32_t * restrict U,
                           float32_t * restrict s,
                           float32_t * restrict V,
                           float32_t * restrict A,
                           int L ) )

size_t rsvd8x8nf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 8==M && 8==N );
  return (0);
}

#endif /* HAVE_VFPU */

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
    Real 4x4 stream ordered matrices
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
Thin SVD For Real/Complex Stream Ordered Matrices

Description: compute the Thin Singular Value Decomposition of L complex 
(real) MxN matrices, with the number of rows greater than or equal to the
number of columns: M>=N. Input and output matrices are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

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
  A[M*N][L]     MxN input matrices
Output:
  U[M*N][L]     MxN matrices comprised of M left-singular column 
                vectors If NULL, this agrument is assumed to be optional 
                and will not be computed
  s[N][L]       Nx1 vectors of singular values in descending order. In an 
                exceptional case when the iterative algorithm fails to 
                converge for a particular matrix, all elements of the
                respective vector are set to NaN.
  V[N*N][L]     NxN matrices comprised of N right-singular column
                vectors. If NULL, this agrument is assumed to be optional 
                and will not be computed
Restrictions:
  pScr,U,s,V,A  Must not overlap and must be aligned on 32-byte boundary 
  M,N           M>=N>1
  L             Must be a multiple of 8
---------------------------------------------------------------------------*/
void rsvd4x4sf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L )
{
  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( U   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( s   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  rsvdmxnsf(pScr,U,s,V,A,4,4,L);

} /* rsvd4x4sf() */

size_t rsvd4x4sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 4==M && 4==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return ( rsvdmxnsf_getScratchSize(M,N,L) );
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, rsvd4x4sf, ( void * pScr,
                       float32_t * restrict U,
                       float32_t * restrict s,
                       float32_t * restrict V,
                       float32_t * restrict A,
                       int L ) )

size_t rsvd4x4sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 4==M && 4==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (0);
}

#endif /* HAVE_VFPU */

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
    NatureDSP_Baseband library. Eigenvalues and eigenvectors
    Complex 16x16 block ordered matrices
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
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
Eigenvalues And Eigenvectors Of Real/Complex Block Ordered Matrices

Description: for each complex/real input matrix A of size NxN, compute N
(possibly repeated) eigenvalues s[N], and (optonally) N right eigenvectors
of size Nx1 V[N]. Input and output data are stored in block order.

Data format: IEEE-754 Std single precision floating-point

Storage sizes SA and Se denote the number of data elements required to store
a matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)
Ss = storage_size(N)

Notes:
1. Functions may perform in-place transformations of input matrices, so that
   INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once the eigenvectors are not required, set the corresponding output pointer
   V to zero, so that a lower complexity algorithm will be used.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.
4. In order to reduce the computational complexity, a preprocessing step known
   as "matrix balancing" is omitted from the implementation.

Temporary:
  pScr        Scratch area. Required size (in bytes) is defined by 
              functions [r]eigen<size>nf_getScratchSize(N,L)
Input:
  N           Matrix size
  L           Number of matrices
  A[L][SA]    NxN input matrices
Output:
  e[L][Se]    Nx1 vectors of eigenvalues. In an exceptional case when the
              iterative algorithm fails to converge for a particular matrix,
              all elements of the respective vector are set to NaN.
  V[L][SA]    NxN matrices comprised of N column eigenvectors (optional)
Restrictions:
  pScr,e,V,A  Must not overlap and must be aligned on 32-byte boundary 
  N           Must be a positive multipe of 4
---------------------------------------------------------------------------*/
void eigen16x16nf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L )
{
  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );

  eigennxnnf(pScr,e,V,A,16,L);

} /* eigen16x16nf() */

size_t eigen16x16nf_getScratchSize( int N, int L )
{
  NASSERT(16==N);
  return ( eigennxnnf_getScratchSize(N,L) );
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, eigen16x16nf, ( void * pScr,
                        complex_float * restrict e,
                        complex_float * restrict V,
                        complex_float * restrict A,
                        int L ) )

size_t eigen16x16nf_getScratchSize ( int N, int L ) 
{
  NASSERT(16==N);
  return (0);
}

#endif /* HAVE_VFPU */

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
    Real MxN block ordered matrices
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_Baseband_matop.h"
/* SVD common declarations */
#include "svd_common.h"

#define MIN(a,b)  ( (a)<(b) ? (a) : (b) )

#define sz_f32   sizeof(float32_t)

#if HAVE_VFPU

/* Calculate the number of data elements occupied by a matrix/vector stored
 * in block order. elemNum is the number of payload data elements, each of
 * elemSize bytes. */
static int getSpace( size_t elemSize, int elemNum )
{
  int stp, wid;
  /* Base-2 log of SIMD width for elemSize-byte elements. */
  wid = LOG2_BBE_SIMD_WIDTH - ( 29 - XT_NSA(elemSize) );
  /* At least one element must fit into a SIMD vector! */
  NASSERT(wid>=0);
  /* Select the allocation step from the number of elements: the next
    * power of two, not greater than the SIMD vector size. */
  stp = 30 - XT_NSA(elemNum);
  if (stp>wid) stp = wid;
  /* Allocation size is the storage size rounded up to the next
    * multiple of allocation step. */ 
  return ((1+((elemNum-1)>>stp))<<stp);
}

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
void rsvdmxnnf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int M, int N, int L )
{
  void *pScr_grsvd;
  float32_t *Ab,*As,*Us,*Vs,*Ds,*Fs;
  int _L,SA;

  NASSERT_ALIGN( U, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( s, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );
  NASSERT( M>=N && N>0 && 0==(M%4) && 0==(N%4) );

  /* Round the number of matrices to the next multiple of SIMD vector width. */
  _L = L + (-L & (BBE_SIMD_WIDTH/2-1));
  /* Space occupoed by single input matrix in block order. */
  SA = getSpace(sz_f32, M*N);

  {
    void * p = pScr;
    /* Partition the scratch area. */
    Ab = (float32_t*)p; p = Ab + SA*_L;
    As = (float32_t*)p; p = As + M*N*_L;
    Ds = (float32_t*)p; p = Ds + N*_L;
    Fs = (float32_t*)p; p = Fs + (N-1)*_L;
    pScr_grsvd = p; p = (uint8_t*)pScr_grsvd + rgrsvdsf_getScratchSize(M,N,_L);
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)rsvdmxnnf_getScratchSize(M,N,L) );
    /* Partially reuse the scratch area. */
    Us = ( U ? As : 0 );
    Vs = ( V ? Ab : 0 );
  }

  /* Copy input matrices to scratch and pad the storage with zeros to avoid 
   * SVD operation on phantom data. */
  memcpy(Ab, A, L*SA*sz_f32); memset(Ab+L*SA, 0, (_L-L)*SA*sz_f32);
  /* Convert input data to stream order. */
  rbsmxnf(As,Ab,M,N,_L);
  /* Reduce input matrixes to upper-bidiagonal form. Left-hand transform matrix
   * U replaces orginial matrix A. */
  rsvd_bidiag_mxnsf(Ds,Fs,Vs,As,M,N,_L,(U!=0));
  /* Apply the Golub-Reinsch SVD to bidiagonal matrices. */
  rgrsvdsf(pScr_grsvd,Ds,Fs,Us,Vs,M,N,_L);
  /* Convert results to block order. */
  svd_rsbmxnxsf(s,Ds,1,N,L,_L);
  if (U) svd_rsbmxnxsf(U,Us,M,N,L,_L);
  if (V) svd_rsbmxnxsf(V,Vs,N,N,L,_L);

} /* rsvdmxnnf() */

size_t rsvdmxnnf_getScratchSize ( int M, int N, int L )
{
  int _L,SA;
  NASSERT( M>=N && N>0 && 0==(M%4) && 0==(N%4) );
  /* Round the number of matrices to the next multiple of SIMD vector width. */
  _L = L + (-L & (BBE_SIMD_WIDTH/2-1));
  /* Space occupoed by single input matrix in block order. */
  SA = getSpace(sz_f32, M*N);
  return ( _L*SA   *sz_f32 + /* Ab: input matrices in block order padded with zeros        */
           _L*M*N  *sz_f32 + /* As: input matrix and left-singular vectors in stream order */
           _L*N    *sz_f32 + /* Ds: main diagonal in stream order                          */
           _L*(N-1)*sz_f32 + /* Fs: superdiagonal in stream order                          */
           rgrsvdsf_getScratchSize(M,N,_L) ); /* Golub-Reinsch SVD implementation req's    */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, rsvdmxnnf, ( void * pScr,
                           float32_t * restrict U,
                           float32_t * restrict s,
                           float32_t * restrict V,
                           float32_t * restrict A,
                           int M, int N, int L ) )

size_t rsvdmxnnf_getScratchSize ( int M, int N, int L )
{
  NASSERT( M>=N && N>0 && 0==(M%4) && 0==(N%4) );
  return (0);
}

#endif /* HAVE_VFPU */

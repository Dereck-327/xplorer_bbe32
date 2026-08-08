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
    Real NxN stream ordered matrices
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Fixed-point arithmetics. */
//#include "NatureDSP_Math.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_Baseband_matop.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

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

/* Scatter-gather copy. Both pointers, steps and element size arguments
 * must be aligned by SIMD vector size (in bytes). */
static void scatter_gather_copy( void * dst, const void * src, 
                                 int dstStep, int srcStep, 
                                 int elemSize, int elemNum )
{
#if 1
    const xb_vecNx16 * Pr = (xb_vecNx16*)src;
    xb_vecNx16 * Pw = (xb_vecNx16*)dst;
    xb_vecNx16 t;
    int n, N = elemSize / (2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(dst, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(src, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (dstStep % (2 * BBE_SIMD_WIDTH)));
    NASSERT(0 == (srcStep % (2 * BBE_SIMD_WIDTH)));
    NASSERT(0 == (elemSize % (2 * BBE_SIMD_WIDTH)));
    while (elemNum-->0) {
        for (n = 0; n<N; n++) {
            BBE_LVNX16_IP(t, Pr, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(t, Pw, 2 * BBE_SIMD_WIDTH);
        }
        Pr = (xb_vecNx16*)((uintptr_t)Pr - N*BBE_SIMD_WIDTH * 2 + srcStep);
        Pw = (xb_vecNx16*)((uintptr_t)Pw - N*BBE_SIMD_WIDTH * 2 + dstStep);
    }
#else
  NASSERT_ALIGN( dst, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( src, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(dstStep%(2*BBE_SIMD_WIDTH)) );
  NASSERT( 0==(srcStep%(2*BBE_SIMD_WIDTH)) );
  NASSERT( 0==(elemSize%(2*BBE_SIMD_WIDTH)) );
  while (elemNum-->0) {
    memcpy( dst, src, elemSize );
    dst = (uint8_t*)dst + dstStep;
    src = (uint8_t*)src + srcStep;
  }
#endif  
} /* scatter_gather_copy() */

#if 0
/* Convert L MxN real matrices from stream to block (packed) order. */
static void _rsbmxnf ( float32_t * restrict y, 
                 const float32_t * restrict x, 
                 int M, int N, int L )
{
  int n,k,S;
  NASSERT_ALIGN( y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( x, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  S = getSpace( sz_f32, M*N );
  for ( k=0; k<L; k++ ) {
    for ( n=0; n<M*N; n++ )
      y[n] = x[n*L];
    x++; y+=S;
  }

} /* _rsbmxnf() */

/* Convert L MxN complex matrices from block (packed) to stream order. */
static void _rbsmxnf ( float32_t * restrict y, 
                 const float32_t * restrict x, 
                 int M, int N, int L )
{
  int n,k,S;
  NASSERT_ALIGN( y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( x, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  S = getSpace( sz_f32, M*N );
  for ( k=0; k<L; k++ ) {
    for ( n=0; n<M*N; n++ )
      y[n*L] = x[n];
    x+=S; y++;
  }

} /* _rbsmxnf() */

/* Convert L MxN complex matrices from block (packed) to stream order. */
static void _cbsmxnf ( complex_float * restrict y, 
                 const complex_float * restrict x, 
                 int M, int N, int L )
{
  int n,k,S;
  NASSERT_ALIGN( y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( x, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  S = getSpace( sz_f32c, M*N );
  for ( k=0; k<L; k++ ) {
    for ( n=0; n<M*N; n++ )
      y[n*L] = x[n];
    x+=S; y++;
  }

} /* _cbsmxnf() */
#endif

/*-------------------------------------------------------------------------
Eigenvalues And Eigenvectors Of Real/Complex Stream Ordered Matrices

Description: for each complex/real input matrix A of size NxN, compute N
(possibly repeated) eigenvalues s[N], and (optonally) N right eigenvectors
of size Nx1 V[N]. Input and output data are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

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
              functions [r]eigen<size>sf_getScratchSize(N,L)
Input:
  N           Matrix size
  L           Number of matrices
  A[N*N][L]   NxN input matrices
Output:
  e[N][L]     Nx1 vectors of eigenvalues. In an exceptional case when the
              iterative algorithm fails to converge for a particular matrix,
              all elements of the respective vector are set to NaN.
  V[N*N][L]   NxN matrices comprised of N column eigenvectors (optional)
Restrictions:
  pScr,e,V,A  Must not overlap and must be aligned on 32-byte boundary 
  N           N>1
  L           Must be a multiple of 8 for real data, or a multiple of 4 for
              complex data
---------------------------------------------------------------------------*/
/* Number of matrices in a processing bunch. */
#define REIGSF_BLK  (BBE_SIMD_WIDTH/2)

void reigennxnsf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int N, int L )
{
  #define BLK  REIGSF_BLK

  float32_t *As,*Hb,*Hs,*Pb,*Ps,*Us;
  complex_float *eb,*es,*Vs;
  int SPb,SHb,Seb;
  int k,p;

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  NASSERT( N>1 );

  SPb = getSpace( sz_f32 , N*N         );
  SHb = getSpace( sz_f32 , N*(N+3)/2-1 );
  Seb = getSpace( sz_f32c, N           );

  {
    void * p = pScr;
    /* Partition the scratch area. */
    As = (float32_t    *)p; p = As + BLK*N*N;
    Hb = (float32_t    *)p; p = Hb + BLK*SHb;
    Pb = (float32_t    *)p; p = Pb + BLK*SPb;
    Vs = (complex_float*)p; p = Vs + BLK*N*N;
    es = (complex_float*)p; p = es + BLK*N;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)reigennxnsf_getScratchSize(N,L) );
    /* Partially reuse the scratch area. */
    Ps = (float32_t*)Vs; 
    eb = (complex_float*)As;
    Us = As; Hs = Pb;
  }

  /* Input matrices are separated into blocks for processing.  */
  for ( k=0; k<L; k+=BLK ) {
    /* Make a local copy of input block to use plain stream order functions */
    scatter_gather_copy( As, A, BLK*sz_f32, L*sz_f32, BLK*sz_f32, N*N );

    if (N>2) {
      /* Matrix size 3x3 and greater */
      if (V) {
        /* In-place reduction of input matrices to upper-Hessenberg form. Keep
         * orthogonal transformation matrices P */
        reigen_hess_nxnsf(Ps,As,N,BLK);
        /* Convert transformation matrices to block order */
        rsbmxnf(Pb,Ps,N,N,BLK);
      } else {
        /* Eigenvectors not required, thus transformation matrices may be dropped. */
        reigen_hess_nxnsf(0,As,N,BLK);
      }
      /* Convert upper-Hessenberg matrices from stream order to compact packed order,
       * with all elements below the 1st subdiagonal dropped. */
      reigen_s2hn_nxnf(Hb,As,N,BLK);
    } else {
      /* 2x2 matrices are already upper-Hessenberg, so simply convert them to
       * block order */
      rsbmxnf(Hb,As,2,2,BLK);
      if (V) {
        /* Initialize transformation matrices with 2x2 identity matrix. */
        for ( p=0; p<BLK; p++ ) {
          Pb[p*SPb+0] = 1.f; Pb[p*SPb+1] = 0.f;
          Pb[p*SPb+2] = 0.f; Pb[p*SPb+3] = 1.f;
        }
      }
    }

    /* Compute eigenpairs for each matrix in a block. */
    for ( p=0; p<BLK; p++ ) {
      if (V) {
        /* Apply the QR algorithm and update the transformation matrix. */
        reigen_hqr_f( eb+p*Seb, Hb+p*SHb, Pb+p*SPb, 0, N-1, N );
        /* Perfrom the backsubstitution to determine eigenvectors of quasi-
         * triangular form lying in Hb */
        reigen_bksubst_f( Hb+p*SHb, eb+p*Seb, N );
      } else {
        /* Use the QR algorithm to estimate eigenvalues; eigenvectors not needed. */
        reigen_hqr_f( eb+p*Seb, Hb+p*SHb, 0, 0, N-1, N );
      }
    }

    /* Convert eigenvalues to stream order. */
    cbsmxnf(es,eb,N,1,BLK);
    /* Copy eigenvalue results to output array. */
    scatter_gather_copy( e, es, L*sz_f32c, BLK*sz_f32c, BLK*sz_f32c, N );

    if (V) {
      /* Convert transformation matrices and eigenvectors of quasi-triangular
       * forms to stream order. */
      rbsmxnf(Ps,Pb,N,N,BLK);
      rbsmxnf(Hs,Hb,1,N*(N+3)/2-1,BLK);
      /* Left-multiply eigenvectors by transformation matrices to obtain re/im
       * components of eigenvectors for original input matrices. */
      reigen_hmulp_nxnsf(Us,Ps,Hs,N,BLK);
      /* Combine re/im components and rescale eigenvectors. */
      reigen_evcomb_nxnsf(Vs,Us,es,N,BLK);
      /* Copy eigenvector results to output array. */
      scatter_gather_copy( V, Vs, L*sz_f32c, BLK*sz_f32c, BLK*sz_f32c, N*N );
    }

    /* Proceed to the next block of input matrices. */
    e += BLK; A += BLK; if (V) V += BLK;
  }

  #undef BLK

} /* reigennxnsf() */

size_t reigennxnsf_getScratchSize ( int N, int L )
{
  NASSERT( N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return ( REIGSF_BLK*        (          N*N         )*sz_f32  +  /* As: a block of input matrices in stream order */
           REIGSF_BLK*getSpace( sz_f32 , N*(N+3)/2-1 )*sz_f32  +  /* Hb: upper-Hessenberg form in compact block order */
           REIGSF_BLK*getSpace( sz_f32 , N*N         )*sz_f32  +  /* Pb: transformation matrix in block order */
           REIGSF_BLK*        (          N*N         )*sz_f32c +  /* Vs: matrix of eigenvectors in stream order */
           REIGSF_BLK*        (          N           )*sz_f32c ); /* es: eigenvalues in stream order */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, reigennxnsf, ( void * pScr,
                         complex_float * restrict e,
                         complex_float * restrict V,
                         float32_t     * restrict A,
                         int N, int L ) )

size_t reigennxnsf_getScratchSize ( int N, int L )
{
  NASSERT( N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (0);
}

#endif /* HAVE_VFPU */

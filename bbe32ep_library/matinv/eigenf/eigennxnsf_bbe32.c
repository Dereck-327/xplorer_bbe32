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
    Complex NxN stream ordered matrices
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
/* Matrix Operations */
#include "NatureDSP_Baseband_matop.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if HAVE_VFPU

#define sz_f32c  sizeof(complex_float)

static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
}

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
  int n, N = elemSize/(2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN( dst, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( src, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(dstStep%(2*BBE_SIMD_WIDTH)) );
  NASSERT( 0==(srcStep%(2*BBE_SIMD_WIDTH)) );
  NASSERT( 0==(elemSize%(2*BBE_SIMD_WIDTH)) );
  while (elemNum-->0) {
    for ( n=0; n<N; n++ ) {
      BBE_LVNX16_IP(t,Pr,2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(t,Pw,2*BBE_SIMD_WIDTH);
    }
    Pr = (xb_vecNx16*)( (uintptr_t)Pr - N*BBE_SIMD_WIDTH*2 + srcStep );
    Pw = (xb_vecNx16*)( (uintptr_t)Pw - N*BBE_SIMD_WIDTH*2 + dstStep );
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
/* Convert L MxN complex matrices from stream to block (packed) order. */
static void _csbmxnf ( complex_float * restrict y, 
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
      y[n] = x[n*L];
    x++; y+=S;
  }

} /* _csbmxnf() */

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
#define EIGSF_BLK  (BBE_SIMD_WIDTH/4)

void eigennxnsf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int N, int L )
{
#if 1
  #define BLK  EIGSF_BLK

  complex_float *As,*Hb,*Ts,*Pb,*Ps;
  complex_float *eb,*es,*Vs;
  int SPb,SHb,Seb;
  int k,p;

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  NASSERT( N>1 );

  SPb = getSpace( sz_f32c, N*N         );
  SHb = getSpace( sz_f32c, N*(N+3)/2-1 );
  Seb = getSpace( sz_f32c, N           );

  {
    void * p = pScr;
    /* Partition the scratch area */
    As = (complex_float*)p; p = As + BLK*N*N;
    Hb = (complex_float*)p; p = Hb + BLK*SHb;
    Pb = (complex_float*)p; p = Pb + BLK*SPb;
    Ps = (complex_float*)p; p = Ps + BLK*N*N;
    eb = (complex_float*)p; p = eb + BLK*Seb;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)eigennxnsf_getScratchSize(N,L) );
    /* Partially reuse the scratch area. */
    es = Hb; Ts = Pb; Vs = As;
  }

  /* Input matrices are separated into blocks for processing.  */
  for ( k=0; k<L; k+=BLK ) {
    /* Make a local copy of input block to use plain stream order functions */
    scatter_gather_copy( As, A, BLK*sz_f32c, L*sz_f32c, BLK*sz_f32c, N*N );

    if (N>2) {
      /* Matrix size 3x3 and greater */
      if (V) {
        /* In-place reduction of input matrices to upper-Hessenberg form. Keep
         * unitary transformation matrices P */
        eigen_hess_nxnsf(Ps,As,N,BLK);
        /* Convert transformation matrices to block order */
        csbmxnf(Pb,Ps,N,N,BLK);
      } else {
        /* Eigenvectors not required, thus transformation matrices may be dropped. */
        eigen_hess_nxnsf(0,As,N,BLK);
      }
      /* Convert upper-Hessenberg matrices from stream order to compact packed order,
       * with all elements below the 1st subdiagonal dropped. */
      eigen_s2hn_nxnf(Hb,As,N,BLK);
    } else {
      /* 2x2 matrices are already upper-Hessenberg, so simply convert them to
       * block order */
      csbmxnf(Hb,As,2,2,BLK);
      if (V) {
        /* Initialize transformation matrices with 2x2 identity matrix. */
        const complex_float c0f = _makecomplexf(0.f,0.f);
        const complex_float c1f = _makecomplexf(1.f,0.f);
        for ( p=0; p<BLK; p++ ) {
          Pb[p*SPb+0] = c1f; Pb[p*SPb+1] = c0f;
          Pb[p*SPb+2] = c0f; Pb[p*SPb+3] = c1f;
        }
      }
    }

    /* Compute eigenpairs for each matrix in a block. */
    for ( p=0; p<BLK; p++ ) {
      if (V) {
        /* Apply the QR algorithm and update the transformation matrix. */
        eigen_hqr_f( eb+p*Seb, Hb+p*SHb, Pb+p*SPb, 0, N-1, N );
        /* Perfrom the backsubstitution to determine eigenvectors of 
         * triangular form lying in Hb */
        eigen_bksubst_f( Hb+p*SHb, N );
      } else {
        /* Use the QR algorithm to estimate eigenvalues; eigenvectors not needed. */
        eigen_hqr_f( eb+p*Seb, Hb+p*SHb, 0, 0, N-1, N );
      }
    }

    if (V) {
      /* Convert transformation matrices and eigenvectors of triangular forms
       * to stream order. */
      cbsmxnf(Ps,Pb,N,N,BLK);
      eigen_hn2ts_nxnf(Ts,Hb,N,BLK);
      /* Left-multiply eigenvectors by transformation matrices and rescale them to
       * obtain eigenvectors for original input matrices. */
      eigen_tmulp_nxnsf(Vs,Ps,Ts,N,BLK);
      /* Copy eigenvector results to output array. */
      scatter_gather_copy( V, Vs, L*sz_f32c, BLK*sz_f32c, BLK*sz_f32c, N*N );
    }

    /* Convert eigenvalues to stream order. */
    cbsmxnf(es,eb,N,1,BLK);
    /* Copy eigenvalue results to output array. */
    scatter_gather_copy( e, es, L*sz_f32c, BLK*sz_f32c, BLK*sz_f32c, N );

    /* Proceed to the next block of input matrices. */
    e += BLK; A += BLK; if (V) V += BLK;
  }

  #undef BLK
#else
  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
#endif
} /* eigennxnsf() */

size_t eigennxnsf_getScratchSize ( int N, int L )
{
  NASSERT( N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  return ( EIGSF_BLK*        (         N*N          )*sz_f32c +  /* As: a block of input matrices in stream order */
           EIGSF_BLK*getSpace( sz_f32c, N*(N+3)/2-1 )*sz_f32c +  /* Hb: upper-Hessenberg form in compact block order */
           EIGSF_BLK*getSpace( sz_f32c, N*N         )*sz_f32c +  /* Pb: transformation matrix in block order */
           EIGSF_BLK*        (          N*N         )*sz_f32c +  /* Ps: transformation matrix in strean order */
           EIGSF_BLK*getSpace( sz_f32c, N           )*sz_f32c ); /* eb: eigenvalues in block order */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, eigennxnsf, ( void * pScr,
                         complex_float * restrict e,
                         complex_float * restrict V,
                         complex_float * restrict A,
                         int N, int L ) )

size_t eigennxnsf_getScratchSize ( int N, int L ) 
{
  NASSERT( N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  return (0);
}

#endif /* HAVE_VFPU */

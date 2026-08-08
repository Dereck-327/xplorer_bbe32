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
    Real NxN block ordered matrices
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

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )

/* Calculate the number of data elements occupied by a matrix/vector stored
* in block order. elemNum is the number of payload data elements, each of
* elemSize bytes. */
static int getSpace(size_t elemSize, int elemNum)
{
    int stp, wid;
    /* Base-2 log of SIMD width for elemSize-byte elements. */
    wid = LOG2_BBE_SIMD_WIDTH - (29 - XT_NSA(elemSize));
    /* At least one element must fit into a SIMD vector! */
    NASSERT(wid >= 0);
    /* Select the allocation step from the number of elements: the next
    * power of two, not greater than the SIMD vector size. */
    stp = 30 - XT_NSA(elemNum);
    if (stp>wid) stp = wid;
    /* Allocation size is the storage size rounded up to the next
    * multiple of allocation step. */
    return ((1 + ((elemNum - 1) >> stp)) << stp);
}

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
#endif

/*
* Convert matrices from stream to block order. Note that the number of
* matrices to be converted and the total number of matrices in the
* input stream are specified through separate input arguments: L and
* stride, respectively.
* Input:
*   M,N             Matrix size
*   L               Number of matrices to be converted
*   stride          Number of matrices in the input stream
*   x[M*N][stride]  Input matrices in stream order
* Output:
*   y[L][SY]        Output matrices in block order
* Restrictions:
*   y,x             Must not overlap and must be aligned on
*                   2*BBE_SIMD_WIDTH-byte boundary
*   stride          Must be a multiple of BBE_SIMD_WIDTH/4
*   L<=stride       Number of converted matrices cannot exceed the number
*                   of matrices in the input stream
*   M,N             The product M*N must be a multiple of 4
* where SY = denotes the number of data entries to store M*N elements in
* block order with proper alignment.
*/

static void csbmxnxsf(complex_float * restrict y,
                const complex_float *          x,
                int M, int N, int L, int stride)
{
    xb_vecN_2xf32 * restrict Y0;
    xb_vecN_2xf32 * restrict Y1;
    xb_vecN_2xf32 * restrict Y2;
    xb_vecN_2xf32 * restrict Y3;
    const xb_vecN_2xf32 * X;
    xb_vecN_2xf32 a0, a1, a2, a3;
    xb_vecN_2xf32 b0, b1, b2, b3;
    vboolN_2 p1, p2, p3;
    int n, k, MN, S;

    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (M*N % (BBE_SIMD_WIDTH / 4)));
    NASSERT(0 == (stride % (BBE_SIMD_WIDTH / 4)));

    MN = M*N;
    S = (MN + BBE_SIMD_WIDTH / 4 - 1) / (BBE_SIMD_WIDTH / 4)*(BBE_SIMD_WIDTH / 4);
    Y0 = (xb_vecN_2xf32*)((uintptr_t)y + 0 * S*sz_f32c);
    Y1 = (xb_vecN_2xf32*)((uintptr_t)y + 1 * S*sz_f32c);
    Y2 = (xb_vecN_2xf32*)((uintptr_t)y + 2 * S*sz_f32c);
    Y3 = (xb_vecN_2xf32*)((uintptr_t)y + 3 * S*sz_f32c);
    for (k = 0; k<L; k += BBE_SIMD_WIDTH / 4) {
        p1 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k + 1), BBE_MOVVA16(L)));
        p2 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k + 2), BBE_MOVVA16(L)));
        p3 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k + 3), BBE_MOVVA16(L)));
        X = (xb_vecN_2xf32*)((uintptr_t)x + k*sz_f32c);
        for (n = 0; n<MN / (BBE_SIMD_WIDTH / 4); n++) {
            BBE_LVN_2XF32_XP(a0, X, stride*sz_f32c);
            BBE_LVN_2XF32_XP(a1, X, stride*sz_f32c);
            BBE_LVN_2XF32_XP(a2, X, stride*sz_f32c);
            BBE_LVN_2XF32_XP(a3, X, stride*sz_f32c);

            BBE_DSELN_2XF32I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_2XF32I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_2XF32I(a1, a0, b2, b0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_2XF32I(a3, a2, b3, b1, BBE_DSELI_INTERLEAVE_4);

            BBE_SVN_2XF32_IP(a0, Y0, BBE_SIMD_WIDTH / 4 * sz_f32c);
            BBE_SVN_2XF32T_IP(a1, Y1, BBE_SIMD_WIDTH / 4 * sz_f32c, p1);
            BBE_SVN_2XF32T_IP(a2, Y2, BBE_SIMD_WIDTH / 4 * sz_f32c, p2);
            BBE_SVN_2XF32T_IP(a3, Y3, BBE_SIMD_WIDTH / 4 * sz_f32c, p3);
        } /* n */

        Y0 = (xb_vecN_2xf32*)((uintptr_t)Y0 + 3 * S*sz_f32c);
        Y1 = (xb_vecN_2xf32*)((uintptr_t)Y1 + 3 * S*sz_f32c);
        Y2 = (xb_vecN_2xf32*)((uintptr_t)Y2 + 3 * S*sz_f32c);
        Y3 = (xb_vecN_2xf32*)((uintptr_t)Y3 + 3 * S*sz_f32c);
    } /* k */
} /* csbmxnxsf() */

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
/* Number of matrices in a processing bunch. */
#define REIGNF_BLK  (BBE_SIMD_WIDTH/2)

void reigennxnnf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int N, int L )
{
  #define BLK  REIGNF_BLK

  float32_t *Ab, *As, *Hb, *Hs, *Pb, *Ps, *Us;
  complex_float *eb, *es, *Vs;
  int SA,Se,SV,SH;
  int k,p,P;

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( N>0 && 0==(N%4) );

  SA = getSpace( sz_f32 , N*N         );
  Se = getSpace( sz_f32c, N           );
  SV = getSpace( sz_f32c, N*N         );
  SH = getSpace( sz_f32 , N*(N+3)/2-1 );

  {
    void * p = pScr;
    /* Partition the scratch area */
    As = (float32_t*)p; p = As + BLK*N*N;
    Pb = (float32_t*)p; p = Pb + BLK*SA;
    Vs = (complex_float*)p; p = Vs + BLK*N*N;
    es = (complex_float*)p; p = es + BLK*N;

    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)reigennxnnf_getScratchSize(N,L) );
    /* Partially reuse the scratch area. */
    Ab = Us = Pb; Hs = As;
    Ps = Hb = (float32_t*)Vs; eb = Vs;
  }

  /* Input matrices are separated into blocks for processing.  */
  for ( k=0; k<L; k+=BLK ) {
    P = MIN(L-k,BLK);
    if (P>=BLK) {
      /* Convert input matrices to stream order. */
      rbsmxnf(As,A,N,N,BLK);
    } else {
      /* Make a local copy of the last incomplete block of input matrices,
       * then safely convert it to stream order. */
      memcpy(Ab,A,P*SA*sz_f32);
      memset(Ab+P*SA,0,(BLK-P)*SA*sz_f32);
      rbsmxnf(As,Ab,N,N,BLK);
    }
    if (V) {
      /* In-place reduction of input matrices to upper-Hessenberg form. Keep
       * orthogonal transformation matrices P */
      reigen_hess_nxnsf(Ps,As,N,BLK);
      rsbmxnf(Pb,Ps,N,N,BLK);
    } else {
      /* Eigenvectors not required, thus transformation matrices may be dropped. */
      reigen_hess_nxnsf(0,As,N,BLK);
    }
    /* Convert upper-Hessenberg matrices from stream order to compact packed order,
     * with all elements below the 1st subdiagonal dropped. */
    reigen_s2hn_nxnf(Hb,As,N,BLK);

    /* Compute eigenpairs for each matrix in a block. */
    for ( p=0; p<P; p++ ) {
      if (V) {
        /* Apply the QR algorithm and update the transformation matrix. */
        reigen_hqr_f( e+p*Se, Hb+p*SH, Pb+p*SA, 0, N-1, N );
        /* Perfrom the backsubstitution to determine eigenvectors of quasi-
         * triangular form lying in Hb */
        reigen_bksubst_f( Hb+p*SH, e+p*Se, N );
      } else {
        /* Use the QR algorithm to estimate eigenvalues; eigenvectors not needed. */
        reigen_hqr_f( e+p*Se, Hb+p*SH, 0, 0, N-1, N );
      }
    }

    if (V) {
        /* Convert transformation matrices and eigenvectors of quasi-triangular
        * forms to stream order. */
        rbsmxnf(Hs, Hb, 1, N*(N + 3) / 2 - 1, BLK);
        rbsmxnf(Ps, Pb, N, N, BLK);
        /* Left-multiply eigenvectors by transformation matrices to obtain re/im
        * components of eigenvectors for original input matrices. */
        reigen_hmulp_nxnsf(Us, Ps, Hs, N, BLK);
        /* Convert eigenvalues to stream order. */
        if (P>=BLK){
          cbsmxnf(es, e, N, 1, BLK);
        } else {
          memcpy(eb, e, P*Se*sz_f32c);
          memset(eb+P*Se, 0, (BLK-P)*Se*sz_f32c);
          cbsmxnf(es, eb, N, 1, BLK);
        }
        /* Combine re/im components and rescale eigenvectors. */
        reigen_evcomb_nxnsf(Vs, Us, es, N, BLK);
        /* Copy eigenvector results to output array. */
        csbmxnxsf(V, Vs, N, N, P, BLK);
    }

    /* Proceed to the next block of input matrices. */
    A += BLK*SA; e += BLK*Se; if (V) V += BLK*SV;
  }

  #undef BLK

} /* reigennxnnf() */

size_t reigennxnnf_getScratchSize ( int N, int L )
{
  NASSERT( N>0 && 0==(N%4) );
  return (REIGNF_BLK*        (         N*N)*sz_f32  + /* As: input matrix in stream order */
          REIGNF_BLK*getSpace(sz_f32 , N*N)*sz_f32  + /* Pb: transformation matrix in block order */
          REIGNF_BLK*        (         N*N)*sz_f32c + /* Vs: matrix of eigenvectors in stream order */
          REIGNF_BLK*        (          N)*sz_f32c); /* es: eigenvalues in stream order */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, reigennxnnf, ( void * pScr,
                         complex_float * restrict e,
                         complex_float * restrict V,
                         float32_t     * restrict A,
                         int N, int L ) )

size_t reigennxnnf_getScratchSize ( int N, int L ) 
{
  NASSERT( N>0 && 0==(N%4) );
  return (0);
}

#endif /* HAVE_VFPU */

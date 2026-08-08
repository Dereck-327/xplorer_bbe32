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
    Extract the upper triangular part of a Hessenberg matrix stored in
    compact block order and convert it to compact stream order.
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Fixed-point arithmetics. */
//#include "NatureDSP_Math.h"
/* Matrix Operations */
#include "NatureDSP_Baseband_matop.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if HAVE_VFPU

#define sz_f32c  sizeof(complex_float)

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )
/* Index of (i,j)-th element of an NxN upper triangular matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * main diagonal aren't actually stored in memory. */
#define TIDX(i,j)   ( (i)*N - (i)*((i)+1)/2 + (j) )

#if 0
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
#endif

/*
 * Extract the upper triangular part of a Hessenberg matrix stored in
 * compact block order and convert it to compact stream order. 
 * Compactness impiles that elements below the first subdiagonal
 * (upper-Hessenberg) or below the main diagonal (upper triangular
 * matrix) are not stored in memory.
 * Number of payload elements in an NxN upper-Hessenberg matrix: N*(N+3)/2-1.
 * Number of elements in the upper triangular part of an NxN matrix, including
 * the main diagonal: N*(N+1)/2.
 * Input:
 *   N                  Matrix size
 *   x[L][SX]           Upper-Hessenberg matrices in block order
 * Output:
 *   y[N*(N+3)/2-1][L]  Upper triangular part of input matrices in stream order.
 *                      After the function completes, results are contained in
 *                      the first N*(N+1)/2*L elements of the array.
 * Restrictions:
 *   y,x  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *        boundary
 *   L    Must be a multiple of BBE_SIMD_WIDTH/4
 * where SX = S(sizeof(complex_float),N*(N+3)/2-1) denotes the number of
 * data elements needed to store all non-zerp elements of an NxN 
 * upper-Hessenberg matrix in block order with proper alignment.
 */

void eigen_hn2ts_nxnf ( complex_float * restrict y,
                  const complex_float * restrict x,
                  int N, int L )
#if 0
{
  int i,j,k,S;
  NASSERT_ALIGN( y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( x, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  S = getSpace( sz_f32c, N*(N+3)/2-1 );
  for ( k=0; k<L; k++ ) {
    for ( i=0; i<N; i++ )
      for ( j=i; j<N; j++ )
        y[TIDX(i,j)*L] = x[HIDX(i,j)];
    x+=S; y++;
  }

} /* eigen_hn2ts_nxnf() */
#else
{
    const xb_vecN_4xcf32 * restrict pY_r;
          xb_vecN_4xcf32 * restrict pY_w;

    const int NN = N*(N + 3) / 2 - 1;
    xb_vecN_4xcf32 val;
    int h_offset;
    int i, j, k;
    
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));
#if 0
    for (k = 0; k<L; k++) {
        for (i = 0; i < NN; i++){
            y[i*L+k] = x[i+S*k];
        }
    }
#else
    cbsmxnf(y, x, 1, NN, L);
#endif
#if 0
    for (k = 0; k<L; k++) {
        for (i = 0; i < N; i++){
            for (j = i; j < N; j++){
                y[TIDX(i, j)*L] = y[HIDX(i, j)*L];
            }
            y++;
        }
    }
#else    
    h_offset = 0;
    for (i = 0; i < N; i++){
        for (j = i; j < N; j++){
            pY_r = (xb_vecN_4xcf32 *)&y[(h_offset + j)*L];
            pY_w = (xb_vecN_4xcf32 *)&y[(h_offset + j - i)*L];
            for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++) {
                BBE_LVN_4XCF32_IP(val, pY_r, sz_f32c*(BBE_SIMD_WIDTH / 4));
                BBE_SVN_4XCF32_IP(val, pY_w, sz_f32c*(BBE_SIMD_WIDTH / 4));
            }
        }
        h_offset += (N - i);
    }
#endif

} /* eigen_hn2ts_nxnf() */
#endif

#endif /* HAVE_VFPU */

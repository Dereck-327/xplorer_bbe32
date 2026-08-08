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
    Left-multiply complex square upper triangular matrices T by unitary
    matrices P, and normalize columns of the product matrix.
    Complex Data, Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define sz_f32c   sizeof(complex_float)

#if 0
/* Index of (i,j)-th element of an NxN upper triangular matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * main diagonal aren't actually stored in memory. */
#define TIDX(i,j)   ( (i)*N - (i)*((i)+1)/2 + (j) )

static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
}

/* Complex floating-point multiplication, single precision */
static complex_float cmulf( complex_float x, complex_float y )
{
  return ( _makecomplexf( crealf(x)*crealf(y) - cimagf(x)*cimagf(y),
                          cimagf(x)*crealf(y) + crealf(x)*cimagf(y) ) );
}
#if 0
/* Multiply a complex floating-point number by an integral power
 * of two, single precision. */
static complex_float cldexpf( complex_float x, int y )
{
  return ( _makecomplexf(ldexpf(crealf(x),y), ldexpf(cimagf(x),y)) );
}
#endif
/* Complex by real floating-point division, single precision. */
static complex_float crdivf( complex_float x, float32_t y )
{
  return ( _makecomplexf( crealf(x)/y, cimagf(x)/y ) );
}

/* Complex floating-point addition, single precision */
static complex_float caddf( complex_float x, complex_float y )
{
  return ( _makecomplexf( crealf(x)+crealf(y), cimagf(x)+cimagf(y) ) );
}

/* Squared absolute of a complex floating-point number, single precision. */
static float32_t cabs2f( complex_float x )
{
  return ( crealf(x)*crealf(x) + cimagf(x)*cimagf(x) );
}
#if 0
/* Cheap estimation of complex number's magnitude, single precision
 * floating-point. */
static float32_t cmagf( complex_float x )
{
  return ( fabsf(crealf(x)) + fabsf(cimagf(x)) );
}
#endif
#endif
/*
 * Left-multiply complex square upper triangular matrices T by unitary
 * matrices P, then scale each column of a product matrix so that its
 * L2 norm is 1, and store results to output matrices: V <- P*T*S, where
 * S is a diagonal matrix with scaling factors on the main diagonal.
 * Note:
 *   TBD If balancing is implemented, the complexity of this functions may
 *       be reduced by taking into account the block diagonal structure of
 *       input matrices. See the MATLAB reference.
 * Input:
 *   N      Matrix size
 *   L      Number of matrices
 *   P[]    Unitary transformation matrices
 *   T[]    Upper triangular matrices. See function definitions for info
 *          on the storage format
 * Output
 *   V[]    Resulting matrices
 * Restrictions:
 *   V,P,H  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *          boundary
 *   Variant functions may impose additional restrictions
 */

/* Complex Data, Stream Order.
 * Elements of T below the MAIN DIAGONAL are not stored
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void eigen_tmulp_nxnsf( complex_float * restrict V, /* V[N*N][L]       */
                  const complex_float * restrict P, /* P[N*N][L]       */
                  const complex_float * restrict T, /* T[N*(N+1)/2][L] */
                  int N, int L )
{
#if 0
  /* A variant with exponent normalization. Use this if the eigenproblem solution
   * includes the balancing procedure. */
  int e,i,j,k,p;
  float32_t f,g;
  complex_float s;
  const complex_float c0f = _makecomplexf(0.f,0.f);
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( P, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( T, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  for ( k=0; k<L; k++ ) {
    for ( j=0; j<N; j++ ) {
      /* Compute the j-th column of P*T, and find the largest magnitude. */
      for ( f=0.f,i=0; i<N; i++ ) {
        for ( s=c0f,p=0; p<=j; p++ ) {
          s = caddf(s, cmulf(P[(i*N+p)*L], T[TIDX(p,j)*L]));
        }
        V[(i*N+j)*L] = s;
        g = cmagf(s); if (f<g) f = g;
      }
      /* Scale the j-th column by max exponent, and compute the L2 norm. */
      e = ilogbf(f); f = 0.f;
      for ( i=0; i<N; i++ ) {
        V[(i*N+j)*L] = s = cldexpf(V[(i*N+j)*L], -e);
        f += cabs2f(s);
      }
      g = sqrtf(f);
      /* Scale the j-th column by its L2 norm. */
      for ( i=0; i<N; i++ ) {
        V[(i*N+j)*L] = crdivf(V[(i*N+j)*L], g);
      }
    }
    V++; P++; T++;
  }
#else
#if 0
  /* A variant without exponent normalization. */
  int i,j,k,p;
  float32_t f,g;
  complex_float s;
  const complex_float c0f = _makecomplexf(0.f,0.f);
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( P, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( T, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  for ( k=0; k<L; k++ ) {
    for ( j=0; j<N; j++ ) {
      /* Compute the j-th column of P*T, and find the largest magnitude. */
      for ( i=0; i<N; i++ ) {
        for ( s=c0f,p=0; p<=j; p++ ) {
          s = caddf(s, cmulf(P[(i*N+p)*L], T[TIDX(p,j)*L]));
        }
        V[(i*N+j)*L] = s;
      }
      /* Scale the j-th column by max exponent, and compute the L2 norm. */
      for ( f=0.f,i=0; i<N; i++ ) {
        V[(i*N+j)*L] = s = V[(i*N+j)*L];
        f += cabs2f(s);
      }
      g = sqrtf(f);
      /* Scale the j-th column by its L2 norm. */
      for ( i=0; i<N; i++ ) {
        V[(i*N+j)*L] = crdivf(V[(i*N+j)*L], g);
      }
    }
    V++; P++; T++;
  }
#else
  /* A variant without exponent normalization. */
  const xb_vecN_4xcf32 * restrict P_r;
  const xb_vecN_4xcf32 * restrict T_r;
  const xb_vecN_4xcf32 * restrict V_r;
        xb_vecN_4xcf32 * V_w;

  xb_vecN_4xcf32 s, s0, s1, s2, s3, x, y;
  xb_vecN_2xf32 z, _f, _g;

  int i, j, k, p;  

  NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(T, 2 * BBE_SIMD_WIDTH);
  NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));

  for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
      for (j = 0; j < N; j++) {
        /* Compute the j-th column of P*T, and find the largest magnitude. */
        V_w = (xb_vecN_4xcf32 *)&V[j*L];
        for (i = 0; i < N; i++) {
            int step = N;
            
            P_r = (xb_vecN_4xcf32 *)&P[i*N*L];
            T_r = (xb_vecN_4xcf32 *)&T[j*L];
#if 0
            s = BBE_ZERON_4XCF32();
            for (p = 0; p <= j; p++) {
                BBE_LVN_4XCF32_XP(x, P_r, L*sz_f32c);
                BBE_LVN_4XCF32_XP(y, T_r, --step*L*sz_f32c);
                BBE_MULAN_4XCF32(s,x,y);
            }
#elif 0
            s0 = s1 = BBE_ZERON_4XCF32();
            for (p = 0; p <= j; p++) {
                BBE_LVN_4XCF32_XP(x, P_r, L*sz_f32c);
                BBE_LVN_4XCF32_XP(y, T_r, --step*L*sz_f32c);

                BBE_MULMASN_4XCF32(s0,x,y,0,4);
                BBE_MULMASN_4XCF32(s1,x,y,1,11);
            }
            s = BBE_ADDN_4XCF32(s0, s1);
#else
            s0 = s1 = s2 = s3 = BBE_ZERON_4XCF32();
            for (p = 0; p < ((j + 1)>>1); p++) {
                BBE_LVN_4XCF32_XP(x, P_r, L*sz_f32c);
                BBE_LVN_4XCF32_XP(y, T_r, --step*L*sz_f32c);

                BBE_MULMASN_4XCF32(s0, x, y, 0, 4);
                BBE_MULMASN_4XCF32(s1, x, y, 1, 11);

                BBE_LVN_4XCF32_XP(x, P_r, L*sz_f32c);
                BBE_LVN_4XCF32_XP(y, T_r, --step*L*sz_f32c);

                BBE_MULMASN_4XCF32(s2, x, y, 0, 4);
                BBE_MULMASN_4XCF32(s3, x, y, 1, 11);
            }

            if ((j+1)&1){
                BBE_LVN_4XCF32_XP(x, P_r, L*sz_f32c);
                BBE_LVN_4XCF32_XP(y, T_r, --step*L*sz_f32c);

                BBE_MULMASN_4XCF32(s0, x, y, 0, 4);
                BBE_MULMASN_4XCF32(s1, x, y, 1, 11);
            }
            s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(s0, s1), BBE_ADDN_4XCF32(s2, s3));
#endif
            BBE_SVN_4XCF32_XP(s, V_w, N*L*sz_f32c);
        }
             
        /* Scale the j-th column by max exponent, and compute the L2 norm. */
        _f = BBE_ZERON_2XF32();
        V_r = (xb_vecN_4xcf32 *)&V[j*L];
        for (i = 0; i < N; i++) {
            BBE_LVN_4XCF32_XP(x, V_r, N*L*sz_f32c);
            z = BBE_MOVN_2XF32_FROMN_4XCF32(x);
            BBE_MULAN_2XF32(_f, z, z);
        }

        _f = BBE_ADDN_2XF32(_f, BBE_SHFLN_2XF32I(_f, BBE_SHFLI_SWAP_2));
        _g = BBE_RSQRTN_2XF32(_f);

        /* Scale the j-th column by its L2 norm. */
        V_r = V_w = (xb_vecN_4xcf32 *)&V[j*L];
        for (i = 0; i < N; i++) {
            BBE_LVN_4XCF32_XP(x, V_r, N*L*sz_f32c);
            z = BBE_MOVN_2XF32_FROMN_4XCF32(x);
            x = BBE_MOVN_4XCF32_FROMN_2XF32( BBE_MULN_2XF32(z, _g) );
            BBE_SVN_4XCF32_XP(x, V_w, N*L*sz_f32c);
        }
    }
    V += (BBE_SIMD_WIDTH / 4); P += (BBE_SIMD_WIDTH / 4); T += (BBE_SIMD_WIDTH / 4);
  }
#endif
#endif
} /* eigen_tmulp_nxnsf() */

#endif /* HAVE_VFPU */

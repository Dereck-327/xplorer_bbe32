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
    Golub-Reinsch SVD algorithm for complex upper bidiagonal matrices stored in
    stream order.
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#if 0
#include <stdio.h>
#include <math.h>
#endif
#include <string.h>
#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Real and complex arithmetic primitives optimized for BBEN VFPU */
#include "vfpu_math.h"
/* NaN values for single precision routines */
#include "nanf_tbl.h"
/* SVD common declarations */
#include "svd_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define EPS          FLT_EPSILON
#define ITS_LIM      75 /* Iterations count limit */

#define MAX(a,b)     ((a)>(b)?(a):(b))
#define SIDX(i,k)    ((i)*L+(k))

#define sz_i16       sizeof(int16_t)
#define sz_f32       sizeof(float32_t)
#define sz_cf32      sizeof(complex_float)
#define sz_vbn4      sizeof(vboolN_4)

/* Internal function of Golub-Kahan SVD step implementation.
 * Summarize index data and fetch data for Wilkinson's shift computation. */
typedef void gks_fetch_fxn_t( 
                    int16_t       * p_l_lo,    int16_t       * p_k_up, 
                    int16_t       * a_its,     vboolN_4      * a_m,
                    complex_float * a_a,       complex_float * a_b, 
                    complex_float * a_x,       complex_float * a_y,
                    complex_float * a_c,       complex_float * a_s,
              const int16_t       * a_l, const int16_t       * a_k,
              const complex_float * D  , const complex_float * F,
              int itsLim, int N, int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Compute Wilkinson's shift. */
typedef void gks_wilkShift_fxn_t(
                    complex_float * a_a,       complex_float * a_b,
                    complex_float * a_c,       complex_float * a_s,
              const complex_float * a_x, const complex_float * a_y,
              int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Compute real Givens's rotator matrix: 
 *   G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) == [*,0], c*conj(c)+s*conj(s) == 1 */
typedef void gks_givens_fxn_t( 
                    complex_float * a_x,        complex_float * a_c,
                    complex_float * a_s,  const complex_float * a_a,
              const complex_float * a_b, 
              int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * 1st updating step of SVD step iteration: B <- B*G(a,b). */
typedef void gks_step1_fxn_t(
                    complex_float * D  ,       complex_float * F,
                    complex_float * a_a,       complex_float * a_b,
              const complex_float * a_c, const complex_float * a_s,
              const int16_t       * a_l, const vboolN_4  * a_m,
              int n, int L );

/*Internal function of Golub-Kahan SVD step implementation.
 * 2nd updating step of SVD step iteration: B <- G(a,b)'*B. */
typedef void gks_step2_fxn_t(
                    complex_float * D,         complex_float * F,
                    complex_float * a_a,       complex_float * a_b,
              const complex_float * a_c, const complex_float * a_s,
              const int16_t       * a_k, const vboolN_4      * a_m,
              int n, int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Accumulate transformation: U <- U*conj(G(a,b)) or V <- V*G(a,b). */
typedef void gks_accumLeft_fxn_t(
                    complex_float * W,   const complex_float * a_c,
              const complex_float * a_s, const vboolN_4      * a_m,
              int n, int M, int N, int L );
typedef void gks_accumRight_fxn_t(
                    complex_float * W,   const complex_float * a_c,
              const complex_float * a_s, const vboolN_4      * a_m,
              int n, int N, int L );

/* Golub-Reinsch SVD scratch area */
typedef struct tag_grsvdsf_scratch
{
  float32_t *tol;
  complex_float *a,*b,*c,*s,*x,*y;
  int16_t *its,*k,*l,*n;
  vboolN_4 *m;

} grsvdsf_scratch_t;

/* Golub-Kahan SVD Step internal functions */
typedef struct tag_grsvdsf_gksApi
{
  gks_fetch_fxn_t      * fetch;
  gks_wilkShift_fxn_t  * wilkShift;
  gks_givens_fxn_t     * givens;
  gks_step1_fxn_t      * step1;
  gks_step2_fxn_t      * step2;
  gks_accumLeft_fxn_t  * accumLeft;
  gks_accumRight_fxn_t * accumRight;

} grsvdsf_gksApi_t;

/* Golub-Kahan SVD Step internal functions set for generic matrix sizes. */
static const grsvdsf_gksApi_t grsvdsf_gksApi_mxn = {
  grsvdsf_gks_fetch    , grsvdsf_gks_wilkShift , grsvdsf_gks_givens,
  grsvdsf_gks_step1    , grsvdsf_gks_step2     ,
  grsvdsf_gks_accumLeft, grsvdsf_gks_accumRight
};

/* Golub-Reinsch SVD for complex upper bidiagonal matrices. Function processes
 * L matrices stored in stream order.
 * Based on [1] Algorithm 8.6.2. */
static void golubReinschSVD( 
                 const grsvdsf_scratch_t * scr,
                 const grsvdsf_gksApi_t  * gksApi,
                       float32_t     * restrict sv,    /* [N][L]   (out)    */
                       complex_float * restrict D,     /* [N][L]   (in/tmp) */
                       complex_float * restrict F,     /* [N-1][L] (in/tmp) */
                       complex_float * restrict U,     /* [M*N][L] (in/out) */
                       complex_float * restrict V,     /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Zero conditionally accessed scratch arrays. */
static void wipeScratch(
                 const grsvdsf_scratch_t * scr,
                 int L );

/* Derive tolerance level from maximum absolute value over a matrix. */
static void setTolLvl( float32_t     * restrict a_tol, /* [L]      (out) */
                 const complex_float * restrict D,     /* [N][L]   (in)  */
                 const complex_float * restrict F,     /* [N-1][L] (in(  */
                 int N, int L );

/* Look for a diagonal submatrix at the bottom and deflate the SVD. */
static void tryDeflate(
                       int16_t       * restrict a_k,   /* [L]      (in/out) */
                       int16_t       * restrict a_its, /* [L]      (in/out) */
                       float32_t     * restrict sv,    /* [N][L]   (in/out) */
                 const complex_float * restrict D,     /* [N][L]   (in)     */
                       complex_float * restrict F,     /* [N-1][L] (in/out) */
                       complex_float * restrict V,     /* [N*N][L] (in/out) */
                 int N, int L );

/* Search for a zero on the superdiagonal or on the main diagonal. If found, 
 * split the matrix. */
static void trySplit(  int16_t       * restrict a_n,   /* [L]      (tmp)    */
                       int16_t       * restrict a_l,   /* [L]      (out)    */
                 const int16_t       * restrict a_k,   /* [L]      (in)     */
                 const float32_t     * restrict a_tol, /* [L]      (in)     */
                       complex_float * restrict D,     /* [N][L]   (in/out) */
                       complex_float * restrict F,     /* [N-1][L] (in/out) */
                       complex_float * restrict U,     /* [M*N][L] (in/out) */
                 int k_up, int M, int N, int L );

/* Golub-Kahan SVD step for principal submatrices in rows/cols a_l[p]..a_k[p].
 * p=0..L-1. Function processes L matrices stored in stream order. 
 * Returns the highest index k among all matrices with a_its[]<ITS_LIM.
 * Based on [1] Algorithm 8.6.1. */
static int16_t golubKahanSVDStep( 
                 const grsvdsf_scratch_t * scr,
                 const grsvdsf_gksApi_t  * gksApi,
                       complex_float * restrict D, /* [N][L]   (in/out) */
                       complex_float * restrict F, /* [N-1][L] (in/out) */
                       complex_float * restrict U, /* [M*N][L] (in/out) */
                       complex_float * restrict V, /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Sort singular values in descending order, and permute left- and right-singular
 * vectors (columns of matrices U and V) accordingly. Data are reordered in-place. */
static void sortSingVal(
                       float32_t     * restrict sv,  /* [N][L]   (in/out) */
                       complex_float * restrict U,   /* [M*N][L] (in/out) */
                       complex_float * restrict V,   /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Compute Givens's rotation matrix.
 * Given scalars a and b, compute scalars s and c such that (in MATLAB notation):
 *   A) [a,b]*[c,conj(s);-s,conj(c)] == [*,0] 
 *   B) c*conj(c)+s*conj(s) == 1 */
static void givens( complex_float * restrict pc, 
                    complex_float * restrict ps, 
              const complex_float * restrict pa, 
              const complex_float * restrict pb );
#if 0
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

/* Complex floating-point multiplication with conjugation of the second
 * argument, single precision */
static complex_float cmuljf( complex_float x, complex_float y )
{
  return ( _makecomplexf( crealf(x)*crealf(y) + cimagf(x)*cimagf(y),
                          cimagf(x)*crealf(y) - crealf(x)*cimagf(y) ) );
}

/* Multiply the first complex argument by conjugate of the second argument,
 * and return the real part of the product. */
static float32_t cmulj_rf( complex_float x, complex_float y )
{
  return ( crealf(x)*crealf(y) + cimagf(x)*cimagf(y) );
}

/* Real by complex floating-point multiplication, single precision. */
static complex_float rcmulf( float32_t x, complex_float y )
{
  return ( _makecomplexf( x*crealf(y), x*cimagf(y) ) );
}

/* Complex floating-point division, single precision. */
static complex_float _cdivf( complex_float x, complex_float y )
{
  /*
   * Based on cdiv code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   * MATLAB reference code:
   *   function z = cdiv(x,y)
   *   xr = real(x); xi = imag(x);
   *   yr = real(y); yi = imag(y);
   *   if abs(yr)>abs(yi)
   *     h = yi/yr; yr = h*yi+yr;
   *     zr = (xr+h*xi)/yr;
   *     zi = (xi-h*xr)/yr;
   *   else
   *     h = yr/yi; yi = h*yr+yi;
   *     zr = (h*xr+xi)/yi;
   *     zi = (h*xi-xr)/yi;
   *   end
   *   z = complex(zr,zi);
   */

  float32_t xr,xi,yr,yi,zr,zi,h;
  xr = crealf(x); xi = cimagf(x);
  yr = crealf(y); yi = cimagf(y);
  if (fabsf(yr)>fabsf(yi)) {
    h = yi/yr; yr = h*yi+yr;
    zr = (xr+h*xi)/yr;
    zi = (xi-h*xr)/yr;
  } else {
    h = yr/yi; yi = h*yr+yi;
    zr = (h*xr+xi)/yi;
    zi = (h*xi-xr)/yi;
  }
  return ( _makecomplexf(zr,zi) );

} /* _cdivf() */

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

/* Complex floating-point subtraction, single precision */
static complex_float csubf( complex_float x, complex_float y )
{
  return ( _makecomplexf( crealf(x)-crealf(y), cimagf(x)-cimagf(y) ) );
}

/* Normalization of a complex number (x/|x|), single precision.
 * Returns 1 for a zero on input. */
static complex_float cnormf( complex_float x )
{
  /*
   * Based on csqrt code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   * MATLAB reference code:
   *   function z = cnorm(x)
   *   xr = abs(real(x));
   *   xi = abs(imag(x));
   *   zr = sign(real(x)); 
   *   zi = sign(imag(x)); 
   *   if xr+xi>0
   *     if xr>xi
   *       p = xi/xr;
   *       zi = zi*p;
   *     elseif xr<xi
   *       p = xr/xi;
   *       zr = zr*p;
   *     else
   *       p = 1; % Protect against Inf/Inf
   *     end
   *     q = 1/sqrt(1+p^2);
   *     z = complex(zr*q,zi*q);
   *   else
   *     z = 1;
   *   end
   */

  float32_t xr,xi,zr,zi,p,q;
  xr = fabsf(crealf(x)); zr = copysignf(1.f, crealf(x));
  xi = fabsf(cimagf(x)); zi = copysignf(1.f, cimagf(x));
  if ((xr+xi)==0) return ( _makecomplexf(1.f,0.f) );
  if (xr>xi) {
    p = xi/xr; zi *= p;
  } else if (xr<xi) {
    p = xr/xi; zr *= p;
  } else {
    p = ( xr==xi ? 1 : xr+xi ); /* Avoid Inf/Inf and arrange for NaN propagation. */
  }
  q = 1.f/sqrtf(1+p*p);
  return ( _makecomplexf(zr*q, zi*q) );

} /* cnormf() */

/* Complex floating-point negation, single precision. */
static complex_float cnegf( complex_float x )
{
  return ( _makecomplexf(-crealf(x), -cimagf(x)) );
}

/* Complex floating-point conjugate, single precision. */
static complex_float _conjf( complex_float x )
{
  return ( _makecomplexf(crealf(x), -cimagf(x)) );
}

/* Complex floating-point absolute value, single precision. */
static float32_t _cabsf( complex_float x )
{
  /*
   * Based on cabs code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   * MATLAB reference code:
   *   function r = cabs(x)
   *   xr = abs(real(x));
   *   xi = abs(imag(x));
   *   if xi>0
   *     if xi>xr, t = xr; xr = xi; xi = t; end;
   *     if xi==xr, t = 1; else t = xi/xr; end;
   *     r = xr*sqrt(1+t^2);
   *   else
   *     r = xr;
   *   end
   */

  float32_t t,xr,xi;
  xr = fabsf(crealf(x));
  xi = fabsf(cimagf(x));
  if (xi>0) {
    if (xi>xr) { t = xr; xr = xi; xi = t; }
    t = ( xi==xr ? 1 : xi/xr );  /* Avoid Inf/Inf */
    return ( xr*sqrtf(1+t*t) );
  } else {
    return ( xi==xi ? xr : xi );
  }

} /* _cabsf() */

/* Squared absolute of a complex floating-point number, single precision. */
static float32_t cabs2f( complex_float x )
{
  return ( crealf(x)*crealf(x) + cimagf(x)*cimagf(x) );
}

/* Cheap estimation of complex number's magnitude, single precision
 * floating-point. */
static float32_t cmagf( complex_float x )
{
  return ( fabsf(crealf(x)) + fabsf(cimagf(x)) );
}
#endif
/*
 * Golub-Reinsch SVD for real or complex upper bidiagonal matrices stored 
 * in stream order. Functions compute the SVD for L matrices, each of M rows
 * and N columns.
 *
 * These functions implement the Golub-Reinsch SVD algorithm as stated in:
 * [1] "Matrix Computations" by G.H. Golub and C.F. Van Loan, 4-th Edition
 * [2] "Singular Value Decomposition and Least Squares Solutions" by 
 *     G.H. Golub abd C. Reinsch, published in "Handbook for Automatie 
 *     Computation", Vol.II, Contribution I/10.
 *
 * Temporary:
 *   pScr            Scratch area.  Required size (in bytes) is defined by 
 *                   functions [r]grsvdsf_getScratchSize(M,N,L)
 * Input:
 *   M,N             Matrix dimensions
 *   L               Number of matrices
 *   F[N-1][L]       First superdiagonal of input matrices. Note that data
 *                   will be distorted by intermediate results.
 * Input/Output:
 *   D[N][L]         Main diagonal of input matrices (in); singular
 *                   values in descending order, or NaNs if failed to
 *                   converge (out, applies to the real-valued variant).
 *   U[M*N][L]       (OPTIONAL) left-hand transformation matrix of the 
 *                   bidiagonalization transform (in); matrix of left-
 *                   singular orthonormal vectors (out).
 *   V[N*N][L]       (OPTIONAL) right-hand transformation matrix of the
 *                   bidiagonalization transform (in); orthogonal (unitary)
 *                   matrix of right-singular vectors (out).
 * Output:
 *   s[N][L]         Singular values in descending order, or NaNs if failed to
 *                   converge (applies to the complex-vlaued variant)
 * Restrictiions:
 *   N>1             Input matrices must have at least 2 columns
 *   M>=N            Number of columns must not exceed the number of rows
 *   L               Must be a multiple of BBE_SIMD_WIDTH/2
 *   pScr,s,D,F,U,V  Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf ( void * pScr,
               float32_t     * restrict s,
               complex_float * restrict D,
               complex_float * restrict F,
               complex_float * restrict U,
               complex_float * restrict V,
               int M, int N, int L )
{
  const grsvdsf_gksApi_t * gksApi = &grsvdsf_gksApi_mxn;
  grsvdsf_scratch_t scr;

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(U, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(V, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));
  NASSERT(M>=N);
  NASSERT(N>1);

  {
    /* Partition the scratch area. */
    int Lix = (L+BBE_SIMD_WIDTH-1)/BBE_SIMD_WIDTH*BBE_SIMD_WIDTH;
    int LW4 = L/(BBE_SIMD_WIDTH/4);
    void * p = pScr;
    memset(&scr, 0, sizeof(scr));
    scr.tol = (float32_t    *)p; p = scr.tol + L;
    scr.a   = (complex_float*)p; p = scr.a + L;
    scr.b   = (complex_float*)p; p = scr.b + L;
    scr.x   = (complex_float*)p; p = scr.x + L;
    scr.y   = (complex_float*)p; p = scr.y + L;
    scr.c   = (complex_float*)p; p = scr.c + L;
    scr.s   = (complex_float*)p; p = scr.s + L;
    scr.its = (int16_t      *)p; p = scr.its + Lix;
    scr.k   = (int16_t      *)p; p = scr.k + Lix;
    scr.l   = (int16_t      *)p; p = scr.l + Lix;
    scr.n   = (int16_t      *)p; p = scr.n + Lix;
    scr.m   = (vboolN_4     *)p; p = scr.m + (N-1)*LW4;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)grsvdsf_getScratchSize(M,N,L) );
  }

  golubReinschSVD(&scr,gksApi,s,D,F,U,V,M,N,L);

} /* grsvdsf() */

size_t grsvdsf_getScratchSize ( int M, int N, int L )
{
  int Lix = (L+BBE_SIMD_WIDTH-1)/BBE_SIMD_WIDTH*BBE_SIMD_WIDTH;
  int LW4 = L/(BBE_SIMD_WIDTH/4);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));
  NASSERT(M>=N);
  NASSERT(N>1);
  return (   L  *sz_f32  +  /* Floating-point data         */
           6*L  *sz_cf32 +  /* Complex floating-point data */
           4*Lix*sz_i16  +  /* Indexing data               */
       (N-1)*LW4*sz_vbn4 ); /* Vector boolean data         */
}

/* Golub-Reinsch SVD for complex upper bidiagonal matrices. Function processes
 * L matrices stored in stream order.
 * Based on [1] Algorithm 8.6.2. */
void golubReinschSVD( const grsvdsf_scratch_t * scr,
                      const grsvdsf_gksApi_t  * gksApi,
                      float32_t     * restrict sv, /* [N][L]   (out)    */
                      complex_float * restrict D,  /* [N][L]   (in/tmp) */
                      complex_float * restrict F,  /* [N-1][L] (in/tmp) */
                      complex_float * restrict U,  /* [M*N][L] (in/out) */
                      complex_float * restrict V,  /* [N*N][L] (in/out) */
                      int M, int N, int L )
{
  /*
   * MATLAB reference code:
   *
   *   function varargout = csvd(A)
   *   nout = max(1,nargout);
   *   needU = nout>1; needV = nout>2;
   *   n = size(A,2);
   *   EPS = eps(cast(1,'like',A));
   *   % Convert input matrix A to bidiagonal form: A = U*B*V'.
   *   if nout>1
   *     if needV, [U,B,V] = cbidiag(A); else [U,B] = cbidiag(A); end;
   *     d = B(1:n+1:n*n); f = B(n+1:n+1:n*n);
   *   else
   *     b = cbidiag(A);
   *     d = b(1,:); f = b(2,2:n);
   *   end
   *   % Preserve the data type for empty matrices because it impacts other
   *   % arguments when passed to a MATLAB function.
   *   if ~needU, U = cast([],'like',A); end;
   *   if ~needV, V = cast([],'like',A); end;
   *   tol = EPS*max(cmag([d,f]));
   *   iterCnt = zeros(1,n);
   *   k = n;
   *   % This loop implements the Golub-Reinsch SVD, [1] Algorithm 8.6.2.
   *   while k>1
   *     % Look for a diagonal submatrix at the bottom.
   *     while k==1 || (k>1 && cmag(f(k-1))<=EPS*(cmag(d(k-1))+cmag(d(k))))
   *       % Tiny elements on the superdiagonal are flushed to zero.
   *       if k>1, f(k-1) = 0; end;
   *       % Normalize the result to obtain a conventional singular value (real,
   *       % non-negative), and adjust the matrix V accordingly.
   *       if ~isempty(V), V(:,k) = V(:,k)*conj(cnorm(d(k))); end;
   *       d(k) = cabs(d(k));
   *       % Deflate the problem
   *       k = k-1;
   *     end
   *     % Search for a zero on the superdiagonal. If found, split the matrix.
   *     l = k;
   *     while l>1
   *       if cmag(f(l-1))<=EPS*(cmag(d(l-1))+cmag(d(l)))
   *         % Set a tiny superdiagonal element to zero and decouple the problem
   *         f(l-1) = 0; break;
   *       elseif cmag(d(l-1))<=tol
   *         % If d(l-1) is zero, set f(l-1) to zero through plane rotations 
   *         % (l-1,l),(l-1,l+1),...,(l-1,k).
   *         c = 1; s = 0;
   *         for p=l:k
   *           a = f(p-1); f(p-1) = -s*a;
   *           a = conj(c)*a; b = d(p);
   *           % G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
   *           % [conj(s),conj(c);c,-s]*[a;b] = [0;*]
   *           [c,s] = givens(a,b);
   *           d(p) = c*a-s*b;
   *           if ~isempty(U)
   *             u = U(:,l-1); v = U(:,p);
   *             % [u,v] <- [u,v]*[conj(s),conj(c);c,-s]'
   *             U(:,l-1) = s*u+c*v;
   *             U(:,p) = conj(c)*u-conj(s)*v;
   *           end
   *         end
   *         % Restart the search.
   *         l = k; continue;
   *       end
   *       l = l-1;
   *     end
   *     if k>l
   *       if iterCnt(k)<75
   *         % Perform a single Golub-Kahan SVD step for B(l:k,l:k)
   *         [d,f,U,V] = GolubKahanStep(d,f,U,V,l,k);
   *         iterCnt(k) = iterCnt(k)+1;
   *       else
   *         fprintf('Solution will not converge\n');
   *         d = []; U = []; V = []; break;
   *       end
   *     end
   *   end
   *   % Sort singular values in descending order. We use the bubble sort due to
   *   % simplicity. Its miserable running time of O(n^2) is completely buried
   *   % under the overall complexity of O(n^3).
   *   if ~isempty(d)
   *     ord = 1:n;
   *     for p=1:n
   *       for q=1:n-p
   *         if d(ord(q))<d(ord(q+1))
   *           t = ord(q); ord(q) = ord(q+1); ord(q+1) = t;
   *         end
   *       end
   *     end
   *     d = d(ord);
   *     if ~isempty(U), U = U(:,ord); end;
   *     if ~isempty(V), V = V(:,ord); end;
   *   end
   *   % Set output argument(s).
   *   if nout>1
   *     if needV, varargout = {U,diag(d),V}; else varargout = {U,diag(d)}; end;
   *   else
   *     varargout = {d'};
   *   end
   */

  float32_t * restrict a_tol = scr->tol;
  int16_t   * restrict a_its = scr->its;
  int16_t   * restrict a_l   = scr->l;
  int16_t   * restrict a_k   = scr->k;
  int16_t   * restrict a_n   = scr->n;

  int k_up,i,p;

#if 0
  static FILE * f_a_l = 0;
  static FILE * f_a_k = 0;
  static FILE * f_a_its = 0;
  static int cnt = 0;
  if (!f_a_l) f_a_l = fopen("a_l.txt", "wt");
  if (!f_a_k) f_a_k = fopen("a_k.txt", "wt");
  if (!f_a_its) f_a_its = fopen("a_its.txt", "wt");
  fprintf(f_a_l, "#%d\n", cnt);
  fprintf(f_a_k, "#%d\n", cnt);
  fprintf(f_a_its, "#%d\n", cnt);
  cnt++;
#endif

  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));
  NASSERT(M>=N);
  NASSERT(N>1);

  /* Zero conditionally accessed scratch arrays. This will prevent alerts
   * from bounds checking tools on reading uninitialized memory locations. */
  wipeScratch(scr,L);
  /* Derive tolerance level from maximum absolute value over a matrix. */
  setTolLvl(a_tol,D,F,N,L);
  /* Start from bottom-right. */
  for ( p=0; p<L; p++ ) {
    a_k[p] = N-1; a_its[p] = 0;
  } /* p */
  /* Iterate until all singular values are revealed, or failed to converge. */
  k_up = N-1;
  while (k_up>0) {
    /* Look for a diagonal submatrix at the bottom. */
    tryDeflate(a_k,a_its,sv,D,F,V,N,L);
    /* Search for a zero on the superdiagonal or on the main diagonal. 
     * If found, split the matrix. */
    trySplit(a_n,a_l,a_k,a_tol,D,F,U,k_up,M,N,L);
#if 0
    for ( p=0; p<L; p++ ) {
      fprintf(f_a_l, "%3d ", (int)scr->l[p]);
      fprintf(f_a_k, "%3d ", (int)scr->k[p]);
      fprintf(f_a_its, "%3d ", (int)scr->its[p]); 
    }
    fprintf(f_a_l, "\n");   fflush(f_a_l);
    fprintf(f_a_k, "\n");   fflush(f_a_k);
    fprintf(f_a_its, "\n"); fflush(f_a_its);
#endif
    /* Perform a single Golub-Kahan SVD step for each matrix */
    k_up = golubKahanSVDStep(scr,gksApi,D,F,U,V,M,N,L);
  } /* k_up */
  /* Check if the algorithm failed to converge for any input matrix. */
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>0) {
      NASSERT(a_its[p]>=ITS_LIM);
      for ( i=0; i<N; i++ ) sv[SIDX(i,p)] = qNaNf.f;
    }
  }
  /* Sort singular values in descending order. */
  sortSingVal(sv,U,V,M,N,L);

} /* golubReinschSVD() */

/* Zero conditionally accessed scratch arrays. */
void wipeScratch(
          const grsvdsf_scratch_t * scr,
          int L )
{
#if 1
  xb_vecN_2xf32 * restrict A_w = (xb_vecN_2xf32*)scr->a;
  xb_vecN_2xf32 * restrict B_w = (xb_vecN_2xf32*)scr->b;
  xb_vecN_2xf32 * restrict C_w = (xb_vecN_2xf32*)scr->c;
  xb_vecN_2xf32 * restrict S_w = (xb_vecN_2xf32*)scr->s;
  xb_vecN_2xf32 * restrict X_w = (xb_vecN_2xf32*)scr->x;
  xb_vecN_2xf32 * restrict Y_w = (xb_vecN_2xf32*)scr->y;
  const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
  int p;
  for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
    BBE_SVN_2XF32_IP(c0f, A_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, B_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, C_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, S_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, X_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, Y_w, 2*BBE_SIMD_WIDTH);
  }
#else
  const complex_float c0f = _makecomplexf(0.f,0.f);
  int p;
  for ( p=0; p<L; p++ ) {
    scr->a[p] = scr->b[p] =
    scr->c[p] = scr->s[p] =
    scr->x[p] = scr->y[p] = c0f;
  }
#endif
} /* wipeScratch() */

/* Derive tolerance level from maximum absolute value over a matrix. */
void setTolLvl( float32_t     * restrict a_tol, /* [L]      (out) */
          const complex_float * restrict D,     /* [N][L]   (in)  */
          const complex_float * restrict F,     /* [N-1][L] (in(  */
          int N, int L )
{
  /*
   * MATLAB outline:
   *   tol = EPS*max(abs(real([d,f]))+abs(imag([d,f])));
   */
#if 1
  const xb_vecN_2xf32 * restrict D_r;
  const xb_vecN_2xf32 * restrict F_r;
  const xb_vecN_2xf32 * restrict TOL_r;
        xb_vecN_2xf32 * restrict TOL_w;

  xb_vecN_2xf32 d0,d1,f0,f1,s,t;
  int i,p;

#if 1
  D_r = (xb_vecN_2xf32*)&D[SIDX(N-1,0)];
  TOL_w = (xb_vecN_2xf32*)a_tol;
  for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
    BBE_LVN_2XF32_IP(d0, D_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(d1, D_r, 2*BBE_SIMD_WIDTH);
    d0 = BBE_ABSN_2XF32(d0); d1 = BBE_ABSN_2XF32(d1);
    BBE_DSELN_2XF32I(d1, d0, d1, d0, BBE_DSELI_DEINTERLEAVE_2);
    t = BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(d1,d0));
    BBE_SVN_2XF32_IP(t, TOL_w, 2*BBE_SIMD_WIDTH);
  } /* p */
  __Pragma("no_reorder");
#else
  for ( p=0; p<L; p++ ) {
    a_tol[p] = cmagf(D[SIDX(N-1,p)])*EPS;
  } /* p */
#endif
#if 1
  D_r = (xb_vecN_2xf32*)&D[SIDX(0,0)];
  F_r = (xb_vecN_2xf32*)&F[SIDX(0,0)];
  for ( i=0; i<N-1; i++ ) {
    TOL_r = TOL_w = (xb_vecN_2xf32*)a_tol;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(d0, D_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(d1, D_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(f0, F_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(f1, F_r, 2*BBE_SIMD_WIDTH);
      d0 = BBE_ABSN_2XF32(d0); d1 = BBE_ABSN_2XF32(d1);
      f0 = BBE_ABSN_2XF32(f0); f1 = BBE_ABSN_2XF32(f1);
      BBE_DSELN_2XF32I(d1, d0, d1, d0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELN_2XF32I(f1, f0, f1, f0, BBE_DSELI_DEINTERLEAVE_2);
      s = BBE_MAXNUMN_2XF32(BBE_ADDN_2XF32(d1,d0), BBE_ADDN_2XF32(f1,f0));
      BBE_LVN_2XF32_IP(t, TOL_r, 2*BBE_SIMD_WIDTH);
      t = BBE_MAXNUMN_2XF32(t, BBE_MULN_2XF32(EPS,s));
      BBE_SVN_2XF32_IP(t, TOL_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  } /* i */
#else
  for ( i=0; i<N-1; i++ ) {
    for ( p=0; p<L; p++ ) {
      a_tol[p] = MAX(a_tol[p], MAX(cmagf(D[SIDX(i,p)]), cmagf(F[SIDX(i,p)]))*EPS);
    } /* p */
  } /* i */
#endif
#else
  int i,p;
  for ( p=0; p<L; p++ ) {
    a_tol[p] = cmagf(D[SIDX(N-1,p)])*EPS;
    for ( i=N-2; i>=0; i-- ) {
      a_tol[p] = MAX(a_tol[p], MAX(cmagf(D[SIDX(i,p)]), cmagf(F[SIDX(i,p)]))*EPS);
    }
  } /* p */
#endif
#if 0 /* !!!! */
  {
    static FILE * f = 0;
    int p;
    if (!f) {
      f = fopen("setTolLvl.txt","wt");
      NASSERT(f);
    }
    for ( p=0; p<L; p++ ) {
      if (isnan(a_tol[p])) {
        fprintf(f,"NaN ");
      } else if (isinf(a_tol[p])) {
        fprintf(f,"Inf ");
      } else {
        fprintf(f, "%.3e ", a_tol[p]);
      }
    }
    fprintf(f,"\n");
    fflush(f);
  }
#endif
} /* setTolLvl() */

/* Look for a diagonal submatrix at the bottom and deflate the SVD. */
void tryDeflate(
                int16_t       * restrict a_k,   /* [L]      (in/out) */
                int16_t       * restrict a_its, /* [L]      (in/out) */
                float32_t     * restrict sv,    /* [N][L]   (in/out) */
          const complex_float * restrict D,     /* [N][L]   (in)     */
                complex_float * restrict F,     /* [N-1][L] (in/out) */
                complex_float * restrict V,     /* [N*N][L] (in/out) */
          int N, int L )
{
  /*
   * MATLAN outline:
   *   % Look for a diagonal submatrix at the bottom.
   *   while k==1 || (k>1 && cmag(f(k-1))<=EPS*(cmag(d(k-1))+cmag(d(k))))
   *     % Tiny elements on the superdiagonal are flushed to zero.
   *     if k>1, f(k-1) = 0; end;
   *     % Normalize the result to obtain a conventional singular value (real,
   *     % non-negative), and adjust the matrix V accordingly.
   *     if ~isempty(V), V(:,k) = V(:,k)*conj(cnorm(d(k))); end;
   *     d(k) = cabs(d(k));
   *     % Deflate the problem
   *     k = k-1;
   *   end
   */
#if 1
  const xtcomplexfloat * restrict _D = (xtcomplexfloat*)D;
        xtcomplexfloat * restrict _F = (xtcomplexfloat*)F;
  const xtcomplexfloat * restrict V_r;
        xtcomplexfloat * restrict V_w;

  float32_t f,g,h;
  xtcomplexfloat s,t;
  vbool1 bnz;
  int i,p;

  const xtcomplexfloat c1f = BBE_CMPLXF32(XT_CONST_S(0), XT_CONST_S(1));

  for ( p=0; p<L; p++ ) {
    while (a_k[p]>=0) {
      if (a_k[p]>0) {
        f = IT_CMAGCF32(_F[SIDX(a_k[p]-1,p)]);
        g = IT_CMAGCF32(_D[SIDX(a_k[p]-1,p)]);
        h = IT_CMAGCF32(_D[SIDX(a_k[p],p)]);
        if (f>EPS*(g+h)) break;
        /* Tiny elements on the superdiagonal are flushed to zero. */
        _F[SIDX(a_k[p]-1,p)] = BBE_CONSTCF32(0);
      }
      /* Normalize the result to obtain a conventional singular value (real,
       * non-negative), and adjust the matrix V accordingly. */
      sv[SIDX(a_k[p],p)] = f = IT_ABSCF32(_D[SIDX(a_k[p],p)]);
      if (V) {
        bnz = XT_UNEQ_S(f, XT_CONST_S(0));
        s = IT_CRDIVCF32(_D[SIDX(a_k[p],p)],f);
        s = BBE_MOVCF32T(s,c1f,bnz);
        V_r = V_w = (xtcomplexfloat*)&V[SIDX(a_k[p],p)];
        for ( i=0; i<N; i++ ) {
          xtcomplexfloat_loadxp(t, V_r, N*L*sz_cf32);
          t = BBE_MULJCF32(t,s);
          xtcomplexfloat_storexp(t, V_w, N*L*sz_cf32);
        }
      } /* V */
      a_its[p] = 0; a_k[p]--;
    } /* a_k[p] */
  } /* p */
#else
  float32_t f,g,h;
  complex_float s;
  int i,p;

  const complex_float c0f = _makecomplexf(0.f,0.f);

  for ( p=0; p<L; p++ ) {
    while (a_k[p]>=0) {
      if (a_k[p]>0) {
        f = cmagf(F[SIDX(a_k[p]-1,p)]); g = cmagf(D[SIDX(a_k[p]-1,p)]); h = cmagf(D[SIDX(a_k[p],p)]);
        if (f>EPS*(g+h)) break;
        /* Tiny elements on the superdiagonal are flushed to zero. */
        F[SIDX(a_k[p]-1,p)] = c0f;
      }
      /* Normalize the result to obtain a conventional singular value (real,
       * non-negative), and adjust the matrix V accordingly. */
      if (V) {
        s = cnormf(D[SIDX(a_k[p],p)]);
        for ( i=0; i<N; i++ ) V[SIDX(i*N+a_k[p],p)] = cmuljf(V[SIDX(i*N+a_k[p],p)], s);
      }
      sv[SIDX(a_k[p],p)] = _cabsf(D[SIDX(a_k[p],p)]);
      a_its[p] = 0; a_k[p]--;
    } /* a_k[p] */
  } /* p */
#endif
} /* tryDeflate() */

/* Search for a zero on the superdiagonal or on the main diagonal. If found, 
 * split the matrix. */
void trySplit(  int16_t       * restrict a_n,   /* [L]      (tmp)    */
                int16_t       * restrict a_l,   /* [L]      (out)    */
          const int16_t       * restrict a_k,   /* [L]      (in)     */
          const float32_t     * restrict a_tol, /* [L]      (in)     */
                complex_float * restrict D,     /* [N][L]   (in/out) */
                complex_float * restrict F,     /* [N-1][L] (in/out) */
                complex_float * restrict U,     /* [M*N][L] (in/out) */
          int k_up, int M, int N, int L )
{
  /*
   * MATLAB outline:
   *   % Search for a zero on the superdiagonal. If found, split the matrix.
   *   l = k;
   *   while l>1
   *     if cmag(f(l-1))<=EPS*(cmag(d(l-1))+cmag(d(l)))
   *       % Set a tiny superdiagonal element to zero and decouple the problem
   *       f(l-1) = 0; break;
   *     elseif cmag(d(l-1))<=tol
   *       % If d(l-1) is zero, set f(l-1) to zero through plane rotations 
   *       % (l-1,l),(l-1,l+1),...,(l-1,k).
   *       c = 1; s = 0;
   *       for p=l:k
   *         a = f(p-1); f(p-1) = -s*a;
   *         a = conj(c)*a; b = d(p);
   *         % G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
   *         % [conj(s),conj(c);c,-s]*[a;b] = [0;*]
   *         [c,s] = givens(a,b);
   *         d(p) = c*a-s*b;
   *         if ~isempty(U)
   *           u = U(:,l-1); v = U(:,p);
   *           % [u,v] <- [u,v]*[conj(s),conj(c);c,-s]'
   *           U(:,l-1) = s*u+c*v;
   *           U(:,p) = conj(c)*u-conj(s)*v;
   *         end
   *       end
   *       % Restart the search.
   *       l = k; continue;
   *     end
   *     l = l-1;
   *   end
   */
#if 1
        xb_vecNx16     * restrict N_w;
        xb_vecNx16     * restrict L_w;
  const xb_vecNx16     * restrict K_r;
  const xb_vecN_2xf32  * restrict TOL_r;
        xtcomplexfloat * restrict _D = (xtcomplexfloat*)D;
  const xb_vecN_2xf32  * restrict D0_r;
  const xb_vecN_2xf32  * restrict D1_r;
        xtcomplexfloat * restrict _F = (xtcomplexfloat*)F;
  const xb_vecN_2xf32  * restrict F_r;
        xb_vecN_2xf32  * restrict F_w;
  const xtcomplexfloat * restrict U0_r;
  const xtcomplexfloat * restrict U1_r;
        xtcomplexfloat * restrict U0_w;
        xtcomplexfloat * restrict U1_w;

  int n_lo,i,j,p;

  const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);

  for ( p=0; p<L; p++ ) {
    a_l[p] = a_n[p] = 0;
  } /* p */
  n_lo = 0;
  do { 
#if 1
    D0_r = (xb_vecN_2xf32*)&D[SIDX(n_lo,0)];
    D1_r = (xb_vecN_2xf32*)&D[SIDX(n_lo+1,0)];
    F_r = F_w = (xb_vecN_2xf32*)&F[SIDX(n_lo,0)];
    for ( i=n_lo+1; i<=k_up; i++ ) {
      xb_vecN_2xf32 a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3;
      xb_vecN_2xf32 f0,f1,g0,g1,h0,h1,t0,t1;
      xb_vecNx16 vi,vk;
      vboolN bz0,bz1,bz2,bz3;
      vboolN_2 bl0,bl1,bn0,bn1;
      vboolN bk,bl,bn;
      vi = BBE_MOVVA16(i);
      N_w = (xb_vecNx16*)a_n;
      L_w = (xb_vecNx16*)a_l;
      K_r = (xb_vecNx16*)a_k;
      TOL_r = (xb_vecN_2xf32*)a_tol;
      for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
        /* i<=a_k[p] */
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        bk = BBE_LENX16(vi,vk);
        /* F[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(a0, F_r, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XF32_IP(a1, F_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a2, F_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a3, F_r, 2*BBE_SIMD_WIDTH);
        /* f <- cmagf(a); */
        a0 = BBE_ABSN_2XF32(a0); a1 = BBE_ABSN_2XF32(a1);
        a2 = BBE_ABSN_2XF32(a2); a3 = BBE_ABSN_2XF32(a3);
        BBE_DSELN_2XF32I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2);
        f0 = BBE_ADDN_2XF32(a0,a1); f1 = BBE_ADDN_2XF32(a2,a3);
        /* D[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(b0, D0_r, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XF32_IP(b1, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(b2, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(b3, D0_r, 2*BBE_SIMD_WIDTH);
        /* g <- cmagf(b); */
        b0 = BBE_ABSN_2XF32(b0); b1 = BBE_ABSN_2XF32(b1);
        b2 = BBE_ABSN_2XF32(b2); b3 = BBE_ABSN_2XF32(b3);
        BBE_DSELN_2XF32I(b1, b0, b1, b0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(b3, b2, b3, b2, BBE_DSELI_DEINTERLEAVE_2);
        g0 = BBE_ADDN_2XF32(b0,b1); g1 = BBE_ADDN_2XF32(b2,b3);
        /* D[SIDX(i,p)] */
        BBE_LVN_2XF32_IP(c0, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(c1, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(c2, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(c3, D1_r, 2*BBE_SIMD_WIDTH);
        /* h <- cmagf(c); */
        c0 = BBE_ABSN_2XF32(c0); c1 = BBE_ABSN_2XF32(c1);
        c2 = BBE_ABSN_2XF32(c2); c3 = BBE_ABSN_2XF32(c3);
        BBE_DSELN_2XF32I(c1, c0, c1, c0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(c3, c2, c3, c2, BBE_DSELI_DEINTERLEAVE_2);
        h0 = BBE_ADDN_2XF32(c0,c1); h1 = BBE_ADDN_2XF32(c2,c3);
        /* f<=EPS*(g+h): zero on the superdiagonal? */
        bl0 = BBE_OLEN_2XF32(f0, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(g0,h0)));
        bl1 = BBE_OLEN_2XF32(f1, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(g1,h1)));
        bl = BBE_ANDBN(bk, BBE_JOINBN_2(bl1,bl0));
        BBE_SVNX16T_IP(vi, L_w, 2*BBE_SIMD_WIDTH, bl);
        /* Tiny superdiagonal elements are flushed to zero. Condition of
         * i<=k does not matter here! */
        BBE_EXTRACTB(bz1, bz0, BBE_MOVN_FROMN_2(bl0));
        BBE_EXTRACTB(bz3, bz2, BBE_MOVN_FROMN_2(bl1));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz0));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz1));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz2));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz3));
        /* g<=a_tol[p]: zero on the main diagonal? */
        BBE_LVN_2XF32_IP(t0, TOL_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t1, TOL_r, 2*BBE_SIMD_WIDTH);
        bn0 = BBE_OLEN_2XF32(g0,t0); bn1 = BBE_OLEN_2XF32(g1,t1);
        bn = BBE_ANDBN(bk, BBE_JOINBN_2(bn1,bn0));
        BBE_SVNX16T_IP(vi, N_w, 2*BBE_SIMD_WIDTH, bn);
      } /* p */
      if ((L&(BBE_SIMD_WIDTH/2))!=0) {
        /* i<=a_k[p] */
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        bk = BBE_LENX16(vi,vk);
        /* F[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(a0, F_r, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XF32_IP(a1, F_r, 2*BBE_SIMD_WIDTH);
        /* D[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(b0, D0_r, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XF32_IP(b1, D0_r, 2*BBE_SIMD_WIDTH);
        /* D[SIDX(i,p)] */
        BBE_LVN_2XF32_IP(c0, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(c1, D1_r, 2*BBE_SIMD_WIDTH);
        /* f <- cmagf(a); */
        a0 = BBE_ABSN_2XF32(a0); a1 = BBE_ABSN_2XF32(a1);
        BBE_DSELN_2XF32I(a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2);
        f0 = BBE_ADDN_2XF32(a0,a1);
        /* g <- cmagf(b); */
        b0 = BBE_ABSN_2XF32(b0); b1 = BBE_ABSN_2XF32(b1);
        BBE_DSELN_2XF32I(b1, b0, b1, b0, BBE_DSELI_DEINTERLEAVE_2);
        g0 = BBE_ADDN_2XF32(b0,b1);
        /* h <- cmagf(c); */
        c0 = BBE_ABSN_2XF32(c0); c1 = BBE_ABSN_2XF32(c1);
        BBE_DSELN_2XF32I(c1, c0, c1, c0, BBE_DSELI_DEINTERLEAVE_2);
        h0 = BBE_ADDN_2XF32(c0,c1);
        /* f<=EPS*(g+h): zero on the superdiagonal? */
        bl0 = BBE_OLEN_2XF32(f0, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(g0,h0)));
        bl1 = BBE_MOVN_2_FROMN(BBE_NEQNX16(vi,vi));
        bl = BBE_ANDBN(bk, BBE_JOINBN_2(bl1,bl0));
        BBE_SVNX16T_IP(vi, L_w, 2*BBE_SIMD_WIDTH, bl);
        /* Tiny superdiagonal elements are flushed to zero. Condition of
         * i<=k does not matter here! */
        BBE_EXTRACTB(bz1, bz0, BBE_MOVN_FROMN_2(bl0));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz0));
        BBE_SVN_2XF32T_IP(c0f, F_w, 2*BBE_SIMD_WIDTH, BBE_MOVN_2_FROMN(bz1));
        /* g<=a_tol[p]: zero on the main diagonal? */
        BBE_LVN_2XF32_IP(t0, TOL_r, 2*BBE_SIMD_WIDTH);
        bn0 = BBE_OLEN_2XF32(g0,t0);
        bn1 = BBE_MOVN_2_FROMN(BBE_NEQNX16(vi,vi));
        bn = BBE_ANDBN(bk, BBE_JOINBN_2(bn1,bn0));
        BBE_SVNX16T_IP(vi, N_w, 2*BBE_SIMD_WIDTH, bn);
      } /* L */
    } /* i */
#else
    for ( i=n_lo+1; i<=k_up; i++ ) {
      float32_t f,g,h;
      const complex_float c0f = _makecomplexf(0.f,0.f);
      for ( p=0; p<L; p++ ) {
        if (i<=a_k[p]) {
          f = cmagf(F[SIDX(i-1,p)]); g = cmagf(D[SIDX(i-1,p)]); h = cmagf(D[SIDX(i,p)]);
          /* Zero on the superdiagonal? */
          if (f<=EPS*(g+h)) {
            F[SIDX(i-1,p)] = c0f; a_l[p] = i;
          }
          /* Zero on the main diagonal? */
          if (g<=a_tol[p]) {
            a_n[p] = i;
          }
        } /* i<=a_k[p] */
      } /* p */
    } /* i */
#endif
    n_lo = k_up+1;
#if 1
    for ( p=0; p<L; p++ ) {
      xtcomplexfloat a,b,c,s,t,u,v,w,z;
      const xtcomplexfloat c0f = BBE_CONSTCF32(0);
      const xtcomplexfloat c1f = BBE_CMPLXF32(XT_CONST_S(0), XT_CONST_S(1));
      if (a_n[p]>a_l[p]) {
        /* If D(n-1) is zero, set F(n-1) to zero through plane rotations
         * (n-1,n),(n-1,n+1),...,(n-1,k). */
        c = c1f; s = c0f;
        for ( j=a_n[p]; j<=a_k[p]; j++ ) {
          a = _F[SIDX(j-1,p)]; 
          _F[SIDX(j-1,p)] = BBE_NEGCF32(BBE_MULCF32(s,a));
          a = BBE_MULCF32(a,c); b = _D[SIDX(j,p)];
          /* G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
           * [conj(s),conj(c);c,-s]*[a;b] = [0;*] */
          givens((complex_float*)&c, (complex_float*)&s, 
                 (complex_float*)&a, (complex_float*)&b);
          t = BBE_MULCF32(c,a); BBE_MULSCF32(t,s,b); _D[SIDX(j,p)] = t;
          if (U) {
            U0_r = U0_w = (xtcomplexfloat*)&U[SIDX(a_n[p]-1,p)];
            U1_r = U1_w = (xtcomplexfloat*)&U[SIDX(j,p)];
            for ( i=0; i<M; i++ ) {
              /* [u,v] <- [u,v]*[conj(s),conj(c);c,-s]' */
              xtcomplexfloat_loadxp(u, U0_r, N*L*sz_cf32);
              xtcomplexfloat_loadxp(v, U1_r, N*L*sz_cf32);
              w = BBE_MULCF32(u,s); BBE_MULACF32(w,v,c);
              z = BBE_MULJCF32(u,c); BBE_MULJSCF32(z,v,s);
              xtcomplexfloat_storexp(w, U0_w, N*L*sz_cf32);
              xtcomplexfloat_storexp(z, U1_w, N*L*sz_cf32);
            }
          } /* U */
        } /* j */
        /* Now F(n-1) is also zero. */
        a_l[p] = a_n[p];
        /* Main diagonal and superdiagonal should be re-examined starting 
         * from position n. */
        if (n_lo>a_n[p]) n_lo = a_n[p];
      } /* a_n[p], a_l[p] */
    } /* p */
#else
    for ( p=0; p<L; p++ ) {
      complex_float c,s,a,b,u,v;
      const complex_float c0f = _makecomplexf(0.f,0.f);
      const complex_float c1f = _makecomplexf(1.f,0.f);
      if (a_n[p]>a_l[p]) {
        /* If D(n-1) is zero, set F(n-1) to zero through plane rotations
         * (n-1,n),(n-1,n+1),...,(n-1,k). */
        c = c1f; s = c0f;
        for ( j=a_n[p]; j<=a_k[p]; j++ ) {
          a = F[SIDX(j-1,p)]; F[SIDX(j-1,p)] = cnegf(cmulf(s,a));
          a = cmulf(a,c); b = D[SIDX(j,p)];
          /* G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
            * [conj(s),conj(c);c,-s]*[a;b] = [0;*] */
          givens(&c,&s,&a,&b);
          D[SIDX(j,p)] = csubf(cmulf(c,a), cmulf(s,b));
          if (U) {
            for ( i=0; i<M; i++ ) {
              /* [u,v] <- [u,v]*[conj(s),conj(c);c,-s]' */
              u = U[SIDX(i*N+a_n[p]-1,p)]; v = U[SIDX(i*N+j,p)];
              U[SIDX(i*N+a_n[p]-1,p)] = caddf(cmulf(s,u), cmulf(c,v));
              U[SIDX(i*N+j,p)] = csubf(cmuljf(u,c), cmuljf(v,s));
            }
          } /* U */
        } /* j */
        /* Now F(n-1) is also zero. */
        a_l[p] = a_n[p];
        /* Main diagonal and superdiagonal should be re-examined starting 
         * from position n. */
        if (n_lo>a_n[p]) n_lo = a_n[p];
      } /* a_n[p], a_l[p] */
    } /* p */
#endif
  } while (n_lo<=k_up);
#else
  complex_float c,s,a,b,u,v;
  float32_t f,g,h;
  int n_lo,i,j,p;

  const complex_float c0f = _makecomplexf(0.f,0.f);
  const complex_float c1f = _makecomplexf(1.f,0.f);

  for ( p=0; p<L; p++ ) {
    a_l[p] = a_n[p] = 0;
  } /* p */
  n_lo = 0;
  do {
    for ( i=n_lo+1; i<=k_up; i++ ) {
      for ( p=0; p<L; p++ ) {
        if (i<=a_k[p]) {
          f = cmagf(F[SIDX(i-1,p)]); g = cmagf(D[SIDX(i-1,p)]); h = cmagf(D[SIDX(i,p)]);
          /* Zero on the superdiagonal? */
          if (f<=EPS*(g+h)) {
            F[SIDX(i-1,p)] = c0f; a_l[p] = i;
          }
          /* Zero on the main diagonal? */
          if (g<=a_tol[p]) {
            a_n[p] = i;
          }
        } /* i<=a_k[p] */
      } /* p */
    } /* i */
    n_lo = k_up+1;
    for ( p=0; p<L; p++ ) {
      if (a_n[p]>a_l[p]) {
        /* If D(n-1) is zero, set F(n-1) to zero through plane rotations
          * (n-1,n),(n-1,n+1),...,(n-1,k). */
        c = c1f; s = c0f;
        for ( j=a_n[p]; j<=a_k[p]; j++ ) {
          a = F[SIDX(j-1,p)]; F[SIDX(j-1,p)] = cnegf(cmulf(s,a));
          a = cmulf(a,c); b = D[SIDX(j,p)];
          /* G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
            * [conj(s),conj(c);c,-s]*[a;b] = [0;*] */
          givens(&c,&s,&a,&b);
          D[SIDX(j,p)] = csubf(cmulf(c,a), cmulf(s,b));
          if (U) {
            for ( i=0; i<M; i++ ) {
              /* [u,v] <- [u,v]*[conj(s),conj(c);c,-s]' */
              u = U[SIDX(i*N+a_n[p]-1,p)]; v = U[SIDX(i*N+j,p)];
              U[SIDX(i*N+a_n[p]-1,p)] = caddf(cmulf(s,u), cmulf(c,v));
              U[SIDX(i*N+j,p)] = csubf(cmuljf(u,c), cmuljf(v,s));
            }
          } /* U */
        } /* j */
        /* Now F(n-1) is also zero. */
        a_l[p] = a_n[p];
        /* Main diagonal and superdiagonal should be re-examined starting 
          * from position n. */
        if (n_lo>a_n[p]) n_lo = a_n[p];
      } /* a_n[p], a_l[p] */
    } /* p */
  } while (n_lo<=k_up);
#endif
} /* trySplit() */

/* Golub-Kahan SVD step for principal submatrices in rows/cols a_l[p]..a_k[p].
 * p=0..L-1. Function processes L matrices stored in stream order. 
 * Returns the highest index k among all matrices with a_its[]<ITS_LIM.
 * Based on [1] Algorithm 8.6.1. */
int16_t golubKahanSVDStep( 
                      const grsvdsf_scratch_t * scr,
                      const grsvdsf_gksApi_t  * gksApi,
                      complex_float * restrict D, /* [N][L]   (in/out) */
                      complex_float * restrict F, /* [N-1][L] (in/out) */
                      complex_float * restrict U, /* [M*N][L] (in/out) */
                      complex_float * restrict V, /* [N*N][L] (in/out) */
                      int M, int N, int L )
{
  /*
   * MATLAB reference code:
   *
   *   function [d,f,U,V] = GolubKahanStep(d,f,U,V,l,k)
   *   % Wilkinson's shift: for the bottom-right 2x2 subblock of B'*B, find an
   *   % eigenvalue which is closer to the last diagonal element. Derived from
   *   % Algol code provided in [2].
   *   x = d(k-1); y = d(k); z = d(l); v = f(k-1);
   *   if k>2, w = f(k-2); else w = 0; end;
   *   t = real((x-y)*conj(x+y)+(w-v)*conj(w+v))/(2*cabs(x)*cabs(v));
   *   r = sqrt(t*t+1); if t<0, r = t-r; else r = t+r; end;
   *   % a <- t11-mu; b <- t12
   *   a = real((z-y)*conj(z+y))+cabs(v)*(cabs(x)/r-cabs(v));
   *   b = conj(d(l))*f(l);
   *   % QR algorithm iteration
   *   x = d(l); y = f(l); 
   *   for p=l:k-1
   *     q = p+1;
   *     % G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0]
   *     [c,s] = givens(a,b);
   *     % B <- B*G(a,b)
   *     if p>l, f(p-1) = c*a-s*b; end;
   *     z = d(q);
   *     a = c*x-s*y; b = -s*z;
   *     y = conj(s)*x+conj(c)*y; z = conj(c)*z;
   *     % V <- V*G(a,b)
   *     if ~isempty(V);
   *       u = V(:,p); v = V(:,q);
   *       V(:,p) = c*u-s*v;
   *       V(:,q) = conj(s)*u+conj(c)*v;
   *     end
   *     % G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) = [*,0] ->
   *     % [c,-s;conj(s),conj(c)]*[a;b] = [*;0]
   *     [c,s] = givens(a,b);
   *     % B <- G(a,b).'*B 
   *     d(p) = c*a-s*b;
   *     a = c*y-s*z; x = conj(s)*y+conj(c)*z;
   *     if p<k-1, y = f(q); b = -s*y; y = conj(c)*y; end;
   *     % U <- U*conj(G(a,b))
   *     if ~isempty(U);
   *       u = U(:,p); v = U(:,q);
   *       U(:,p) = conj(c)*u-conj(s)*v;
   *       U(:,q) = s*u+c*v;
   *     end
   *   end
   *   f(k-1) = a;
   *   d(k) = x;
   */

  int16_t l_lo,k_up;
  int n;
  gksApi->fetch(&l_lo, &k_up, scr->its, scr->m, scr->a, scr->b, 
                scr->x, scr->y, scr->c, scr->s, scr->l, scr->k, 
                D, F, ITS_LIM, N, L);
  if (l_lo<k_up) {
    gksApi->wilkShift(scr->a, scr->b, scr->c, scr->s, scr->x, scr->y, L);
    for ( n=l_lo; n<k_up; n++ )      {
      gksApi->givens            (scr->x, scr->c, scr->s, scr->a, scr->b                 , L);
      gksApi->step1             (D, F, scr->a, scr->b, scr->c, scr->s, scr->l, scr->m, n, L);
      if (V) gksApi->accumRight (V, scr->c, scr->s, scr->m, n, N                        , L);
      gksApi->givens            (scr->x, scr->c, scr->s, scr->a, scr->b                 , L);
      gksApi->step2             (D, F, scr->a, scr->b, scr->c, scr->s, scr->k, scr->m, n, L);
      if (U) gksApi->accumLeft  (U, scr->c, scr->s, scr->m, n, M, N                     , L);
    } /* n */
  } /* l_lo, k_up */
  return (k_up);

} /* golubKahanSVDStep() */

/* Sort singular values in descending order, and permute left- and right-singular
 * vectors (columns of matrices U and V) accordingly. Data are reordered in-place. */
void sortSingVal(
                float32_t     * restrict sv,  /* [N][L]   (in/out) */
                complex_float * restrict U,   /* [M*N][L] (in/out) */
                complex_float * restrict V,   /* [N*N][L] (in/out) */
                int M, int N, int L )
{
  float32_t f;
  complex_float s;
  int k,l,n,p;
  for ( p=0; p<L; p++ ) {
    /* We use a kind of the bubble sort due to its simplicity. Normally, the number
     * of singular N is too small for an advanced sorting algorithm to demonstrate
     * its power. */
    for ( k=0; k<N-1; k++ ) {
      l = k; f = sv[SIDX(k,p)];
      for ( n=k+1; n<N; n++ ) {
        if (sv[SIDX(n,p)]>f) { l = n; f = sv[SIDX(n,p)]; }
      }
      if (k!=l) {
        f = sv[SIDX(l,p)]; sv[SIDX(l,p)] = sv[SIDX(k,p)]; sv[SIDX(k,p)] = f;
        if (U) {
          for ( n=0; n<M; n++ ) {
            s = U[SIDX(n*N+l,p)]; U[SIDX(n*N+l,p)] = U[SIDX(n*N+k,p)]; U[SIDX(n*N+k,p)] = s;
          }
        }
        if (V) {
          for ( n=0; n<N; n++ ) {
            s = V[SIDX(n*N+l,p)]; V[SIDX(n*N+l,p)] = V[SIDX(n*N+k,p)]; V[SIDX(n*N+k,p)] = s;
          }
        }
      } /* k!=l */
    } /* k */
  } /* p */

} /* sortSingVal() */

/* Compute Givens's rotation matrix.
 * Given scalars a and b, compute scalars s and c such that (in MATLAB notation):
 *   A) [a,b]*[c,conj(s);-s,conj(c)] == [*,0] and 
 *   B) c*conj(c)+s*conj(s) == 1 */
void givens( complex_float * restrict pc, 
             complex_float * restrict ps, 
       const complex_float * restrict pa, 
       const complex_float * restrict pb )
{
  /*
   * function [c,s] = givens(a,b)
   * if b~=0
   *   if b*conj(b)>a*conj(a)
   *     t = cdiv(a,b);
   *     s = 1/sqrt(1+t*conj(t));
   *     c = -conj(t)*s;
   *   else
   *     t = cdiv(b,a); 
   *     c = 1/sqrt(1+t*conj(t));
   *     s = -conj(t)*c;
   *   end
   * else
   *   c = 1; s = 0;
   * end;
   */
#if 1
  xtcomplexfloat _c,_s;
  xb_vecN_2xf32 a,b,c,r,s,t,u,v;
  vboolN_2 bnz,ble;
  vboolN bn0,bn1;
  vselN sel;
  int dum;

  const xb_vecN_2xf32 c1f = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0), 
                                            BBE_CONSTN_2XF32(1),
                                            BBE_SELI_INTERLEAVE_2_LO);
  const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);

  u = BBE_MOVN_2XF32_FROMCF32(MOV_CF32_FROM_COMPLEX_FLOAT(*pa));
  v = BBE_MOVN_2XF32_FROMCF32(MOV_CF32_FROM_COMPLEX_FLOAT(*pb));
  a = BBE_SELN_2XF32I(v,u,BBE_SELI_INTERLEAVE_4_LO);
  b = BBE_SELN_2XF32I(u,v,BBE_SELI_INTERLEAVE_4_LO);

  u = BBE_MULN_2XF32(a,a); v = BBE_MULN_2XF32(b,b);
  u = BBE_ADDN_2XF32(u, BBE_SHFLN_2XF32I(u, BBE_SHFLI_SWAP_2));
  v = BBE_ADDN_2XF32(v, BBE_SHFLN_2XF32I(v, BBE_SHFLI_SWAP_2));
  ble = BBE_ULEN_2XF32(u,v);
  u = BBE_ABSN_2XF32(b);
  u = BBE_ADDN_2XF32(u, BBE_SHFLN_2XF32I(u, BBE_SHFLI_SWAP_2));
  v = BBE_SHFLN_2XF32I(a, BBE_SHFLI_SWAP_2);
  bnz = BBE_ORBN_2(BBE_UNEQN_2XF32(u,c0f), BBE_UNN_2XF32(a,v));
  BBE_EXTRACTB(bn1, bn0, BBE_MOVN_FROMN_2(bnz));
  bnz = BBE_MOVN_2_FROMN(bn0);

  t = IT_CDIVN_2XF32(a,b,0);
  r = c1f; BBE_MULAN_2XF32(r,t,t);
  r = BBE_ADDN_2XF32(r, BBE_SHFLN_2XF32I(r, BBE_SHFLI_SWAP_2));
  r = BBE_RSQRTN_2XF32(r);
  s = BBE_SELN_2XF32I(c0f, r, BBE_SELI_INTERLEAVE_2_EVEN);
  c = BBE_MULN_2XF32(BBE_CONJN_2XF32(BBE_NEGN_2XF32(t)),r);

  BBE_SQZN(sel, dum, BBE_MOVN_FROMN_2(ble)); (void)dum;
  u = BBE_SHFLN_2XF32(c, BBE_MOVVSELN_2_FROMVSELN(sel));
  v = BBE_SHFLN_2XF32(s, BBE_MOVVSELN_2_FROMVSELN(sel));
  c = BBE_MOVN_2XF32T(u,v,ble);
  s = BBE_MOVN_2XF32T(v,u,ble);
  c = BBE_MOVN_2XF32T(c,c1f,bnz);
  s = BBE_MOVN_2XF32T(s,c0f,bnz);

  _c = BBE_MOVCF32_FROMN_2XF32(c);
  _s = BBE_MOVCF32_FROMN_2XF32(s);
  *pc = MOV_COMPLEX_FLOAT_FROM_CF32(_c);
  *ps = MOV_COMPLEX_FLOAT_FROM_CF32(_s);
#else
  complex_float a,b,c,s,t;
  float32_t r;
  a = *pa; b = *pb;
  /* Use inverted condition to arrange for NaN propagation. */
  if (!(cmagf(b)==0) || !isfinite(cmagf(a))) {
    if (cabs2f(b)>cabs2f(a)) {
      t = _cdivf(a,b);
      r = 1.f/sqrtf(1.f + cabs2f(t));
      s = _makecomplexf(r,0.f);
      c = _conjf(rcmulf(-r,t));
    } else {
      t = _cdivf(b,a);
      r = 1.f/sqrtf(1.f + cabs2f(t));
      c = _makecomplexf(r,0.f);
      s = _conjf(rcmulf(-r,t));
    }
  } else {
    c = _makecomplexf(1.f,0.f);
    s = _makecomplexf(0.f,0.f);
  }
  *pc = c; *ps = s;
#endif
} /* givens() */

#if 0
/* Sanity check of the givens() function against NaNs on input. */
int main(void)
{
  float32_t cr,ci,sr,si;
  #define GIVENS(ar,ai,br,bi) \
          { \
            xtcomplexfloat a,b,c,s; \
            a = BBE_CMPLXF32(ai,ar); b = BBE_CMPLXF32(bi,br); \
            givens((complex_float*)&c,(complex_float*)&s, \
                   (complex_float*)&a,(complex_float*)&b); \
            cr = BBE_CREALCF32(c); ci = BBE_CIMAGCF32(c); \
            sr = BBE_CREALCF32(s); si = BBE_CIMAGCF32(s); \
          }
  GIVENS(qNaNf.f,  -2.5f,  +0.6f,  +3.8f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +1.7f,qNaNf.f,  +0.6f,  +3.8f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +1.7f,  -2.5f,qNaNf.f,  +3.8f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +1.7f,  -2.5f,  +0.6f,qNaNf.f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(qNaNf.f,qNaNf.f,  +1.7f,  -2.5f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(qNaNf.f,  +3.8f,qNaNf.f,  -2.5f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(qNaNf.f,  +3.8f,  +1.7f,qNaNf.f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.6f,qNaNf.f,qNaNf.f,  -2.5f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.6f,qNaNf.f,  +1.7f,qNaNf.f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.6f,  +3.8f,qNaNf.f,qNaNf.f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.0f,  +0.0f,  +0.0f,   0.0f); NASSERT((cr==1.f)&&(ci==0.f)&&(sr==0.f)&&(si==0.f));
  GIVENS(qNaNf.f,  +0.0f,  +0.0f,   0.0f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.0f,qNaNf.f,  +0.0f,   0.0f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.0f,  +0.0f,qNaNf.f,   0.0f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  GIVENS(  +0.0f,  +0.0f,  +0.0f,qNaNf.f); NASSERT((isnan(cr)||isnan(ci))&&(isnan(sr)||isnan(si)));
  #undef GIVENS
  return (0);
}
#endif

#endif /* HAVE_VFPU */

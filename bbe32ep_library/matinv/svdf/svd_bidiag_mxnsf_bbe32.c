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
    Reduction of generic matrices to upper bidiagonal form
    Complex Data, Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <float.h>
#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Real and complex arithmetic primitives optimized for BBEN VFPU */
#include "vfpu_math.h"
/* SVD common declarations */
#include "svd_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define sz_f32    sizeof(float32_t)
#define sz_f32c   sizeof(complex_float)

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
 * Reduce MxN matrices with M>=N to upper bidiagonal form by left- and
 * right-hand Householder transformations.
 * Original matrix A, upper bidiagonal form B, left- and right-hand 
 * transformation matrices U and V relate to each other through the 
 * following identities: A == U*B*V', U'*A*V == B.
 * Matrix U is comprised of orthonormal columns, matrix V is orthogonal
*  (unitary for complex data).
 * Input:
 *   M,N         Matrix dimensions
 *   L           Number of matrices
 *   needU       Set to non-zero value if left-hand transformation matrices
 *               U are required (applies to stream order functions)
 * Input/Output:
 *   A[M*N]xL    Original matrices (in); left-hand transformation matrices
 *               U (out, optional, stream order functions). Even when 
 *               not utilized to keep the resulting matrix U, A is still used
 *               as an intermediate storage, so input data are damaged in any
 *               case.
 * Output:
 *   D[N]xL      Main diagonal of reduced matrices
 *   F[N-1]xL    First superdiagonal of reduced matrices
 *   U[M*N]xL    Left-hand transformation matrix (block order functions). Set
 *               to zero if not needed.
 *   V[N*N]xL    Orthogonal matrix of accumulated right-hand transformations,
 *               optional. Set to zero if not needed.
 * Restrictiions:
 *   D,F,U,V,A   Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *               boundary
 *   N>1         Input matrices must have at least 2 columns
 *   M>=N        Number of columns must not exceed the number of rows
 *   M,N,L       Subject to additional limitations exposed by a particular
 *               function
 */

/* Complex-valued, Stream Order
 * Note that output matrices U replace input matrices A.
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void svd_bidiag_mxnsf (
                complex_float * restrict D, /* D[N][L]   */
                complex_float * restrict F, /* F[N-1][L] */
                complex_float * restrict V, /* V[N*N][L] */
                complex_float * restrict A, /* A[M*N][L] */
                int M, int N, int L, int needU )
{
  /*
   * This implementation is based on the code from "Singular Value 
   * Decomposition and Least Squares Solutions" by G.H. Golub and 
   * C. Reinsch, Handbook for Automatie Computation, Vol.II Linear
   * Algebra, Contribution I/10.
   *
   * MATLAB reference code with a simple test:
   *
   *   M = 7; N = 5;
   *   A = single(complex(1-2*rand(M,N),1-2*rand(M,N))); R = A;
   *   % Perform the bidiagonalization
   *   tol = realmin(class(A));
   *   D = cast(zeros(1,N),'like',A);
   *   F = cast(zeros(1,N-1),'like',A);
   *   for i=1:N
   *     l = i+1;
   *     % Left-hand transformation
   *     t = A(i:M,i)'*A(i:M,i);
   *     if t>=tol
   *       r = sqrt(t);
   *       f = A(i,i); w = abs(f);
   *       if w>0, g = -r*(f/w); else g = -r; end;
   *       h = -t-w*r; A(i,i) = f-g;
   *       for j=l:N
   *         s = A(i:M,i)'*A(i:M,j)/h;
   *         A(i:M,j) = A(i:M,j)+s*A(i:M,i);
   *       end
   *     else
   *       g = 0;
   *     end
   *     D(i) = g;
   *     % Right-hand transformation
   *     if i<N-1
   *       t = A(i,l:N)*A(i,l:N)';
   *       if t>=tol
   *         r = sqrt(t);
   *         f = A(i,l); w = abs(f);
   *         if w>0, g = -r*(f/w); else g = -r; end;
   *         h = -t-w*r;  A(i,l) = f-g;
   *         for j=l:M
   *           s = A(j,l:N)*A(i,l:N)'/h;
   *           A(j,l:N) = A(j,l:N)+s*A(i,l:N);
   *         end
   *       else
   *         g = 0;
   *       end
   *       F(i) = g; 
   *     elseif i==N-1
   *       F(i) = A(i,l);
   *     end
   *   end
   *   % Accumulate right-hand transformations from reflection axes stored in 
   *   % the work matrix A.
   *   V = cast(complex(1-2*rand(N),1-2*rand(N)),'like',A);
   *   V(N-1,N-1) = 1; V(N,N) = 1;
   *   for i=N-2:-1:1
   *     l = i+1; g = F(i);
   *     if g~=0
   *       h = real(A(i,l)*conj(g));
   *       for j=l+1:N
   *         s = A(i,l+1:N)*V(l+1:N,j)/h;
   *         V(l,j) = s*A(i,l)';
   *         V(l+1:N,j) = V(l+1:N,j)+s*A(i,l+1:N)';
   *       end
   *       s = A(i,l)/h;
   *       V(1:N,l) = s*A(i,1:N)';
   *     else
   *       V(l:N,l) = 0; V(l,l+1:N) = 0;
   *     end
   *     V(l,l) = V(l,l)+1;
   *   end
   *   V(1,1) = 1; V(1,2:N) = 0; V(2:N,1) = 0;
   *   % Accumulate left-hand transformations. Here we build matrix U over 
   *   % intermediate data stored in A.
   *   for i=N:-1:1
   *     l = i+1; g = D(i);
   *     if g~=0
   *       h = real(A(i,i)*conj(g));
   *       for j=l:N
   *         s = A(l:M,i)'*A(l:M,j)/h;
   *         A(i,j) = s*A(i,i);
   *         A(l:M,j) = A(l:M,j)+s*A(l:M,i);
   *       end
   *       s = conj(A(i,i))/h;
   *       A(i:M,i) = s*A(i:M,i);
   *     else
   *       A(i:M,i) = 0; A(i,i+1:N) = 0;
   *     end
   *     A(i,i) = A(i,i)+1;
   *   end
   *   % Verify the resulting decomposition of the test matrix.
   *   B = diag(D)+diag(F,1);
   *   R = R*V-A*B;
   *   fprintf('M(R) = %.1e\n',max(abs(R(:))));
   */

        xb_vecN_4xcf32 * restrict D_rw = (xb_vecN_4xcf32*)D;
        xb_vecN_4xcf32 * restrict F_rw = (xb_vecN_4xcf32*)F;
        xb_vecN_4xcf32 * restrict V_rw = (xb_vecN_4xcf32*)V;
  const xb_vecN_4xcf32 *          V_r;
        xb_vecN_4xcf32 * restrict V_w0;
        xb_vecN_4xcf32 * restrict V_w1;
  const xb_vecN_4xcf32 *          A_r0;
  const xb_vecN_4xcf32 *          A_r1;
        xb_vecN_4xcf32 * restrict A_w;
        xb_vecN_4xcf32 * restrict A_rw = (xb_vecN_4xcf32*)A;

  const xb_vecN_2xf32  tol  = FLT_MIN;
  const xb_vecN_2xf32  c0f  = BBE_CONSTN_2XF32(0);
  const xb_vecN_4xcf32 c0cf = BBE_CONSTN_4XCF32(0);
  const xb_vecN_4xcf32 c1cf = BBE_SELN_4XCF32I(c0cf, BBE_CONSTN_4XCF32(1), BBE_SELI_INTERLEAVE_2_LO);

  xb_vecN_4xcf32 f,g,s,u,v,x,y;
  xb_vecN_2xf32 h;
  vboolN_4 btol,bgtz;
  int i,j,k,l,p;

  #define SIDX(i,j,k)  (((i)*N+(j))*L+(k))

  NASSERT_ALIGN( D, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( F, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );
  NASSERT( M>=N && N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );

  for ( k=0; k<L; k+=BBE_SIMD_WIDTH/4 ) {
    /*
     * Perform the bidiagonalization.
     */
    for ( i=0; i<N; i++ ) {
      l = i+1;
      /*======================================================================*
       *                      Left-hand transformation                        *
       *======================================================================*/
#if 1
      {
        xb_vecN_2xf32 r,t,w;
        /* for ( t=0.f, j=i; j<M; j++ ) {
         *   t += cabs2f(A[SIDX(j,i,k)]);
         * } */
        A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,i,k)];
        __Pragma("loop_count min=1");
        for ( t=c0f, j=i; j<M; j++ ) {
          BBE_LVN_4XCF32_XP(f, A_r0, N*L*sz_f32c);
          r = BBE_MOVN_2XF32_FROMN_4XCF32(f);
          BBE_MULAN_2XF32(t,r,r);
        }
        t = BBE_ADDN_2XF32(t, BBE_SHFLN_2XF32I(t, BBE_SHFLI_SWAP_2));
        /* Mask out columns whose L2 norm is finite, but so small that
         * it cannot be evaluated to working precision. */
        btol = BBE_MOVN_4_FROMN_2(BBE_ULEN_2XF32(tol,t));
        /*---------------------------------------------------------------------* 
         * Construct a Householder vector                                      */
        /* r = sqrtf(t); */
        r = BBE_NEGN_2XF32(BBE_FSQRTN_2XF32(t)); /* !!!! */
        /* f = A[SIDX(i,i,k)]; w = _cabsf(f); */
        f = BBE_LVN_4XCF32_X(A_rw, SIDX(i,i,k)*sz_f32c);
        w = IT_ABSN_4XCF32(f);
        /*  g = (w>0 ? rcmulf(-r, crdivf(f,w)) : _makecomplexf(-r,0) ); */
        u = IT_RCMULN_4XCF32(r, IT_CRDIVN_4XCF32(f,w));
        bgtz = BBE_MOVN_4_FROMN_2(BBE_OLTN_2XF32(c0f, w));
        v = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f,r,BBE_SELI_INTERLEAVE_2_EVEN));
        g = BBE_MOVN_4XCF32T(u,v,bgtz);
        g = BBE_MOVN_4XCF32T(g,c0cf,btol);
        BBE_SVN_4XCF32_X(g, D_rw, SIDX(0,i,k)*sz_f32c);
        /* h = -t-w*r; A[SIDX(i,i,k)] = csubf(f,g); */
        h = BBE_NEGN_2XF32(t); BBE_MULAN_2XF32(h,w,r);
        BBE_SVN_4XCF32_X(BBE_SUBN_4XCF32(f,g), A_rw, SIDX(i,i,k)*sz_f32c);
        /*---------------------------------------------------------------------* 
         * Apply the Householder reflector from the left to rows i..M-1,       *
         * columns i+1..N-1                                                    */
        for ( j=l; j<N; j++ ) {
          /* for ( s=c0f, p=i; p<M; p++ ) {
           *   s = caddf(s, cmuljf(A[SIDX(p,j,k)], A[SIDX(p,i,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,j,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(i,i,k)];
          for ( f=g=u=v=c0cf, p=0; p<(M-i)/2; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULMASN_4XCF32(f,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(g,x,y,0x2,0x7);
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULMASN_4XCF32(u,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(v,x,y,0x2,0x7);
          }
          if ((M-i)&1) {
            x = BBE_LVN_4XCF32_I(A_r0,0);
            y = BBE_LVN_4XCF32_I(A_r1,0);
            BBE_MULMASN_4XCF32(f,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(g,x,y,0x2,0x7);
          }
          s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(f,g), BBE_ADDN_4XCF32(u,v));
          /* s = crdivf(s,h); */
          s = IT_CRDIVN_4XCF32(s,h);
          /* for ( p=i; p<M; p++ ) {
           *   A[SIDX(p,j,k)] = caddf(A[SIDX(p,j,k)], cmulf(s, A[SIDX(p,i,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,j,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(i,i,k)];
          A_w  = (xb_vecN_4xcf32*)&A[SIDX(i,j,k)];
          __Pragma("loop_count min=1");
          __Pragma("super_swp ii=12, unroll=4");
          for ( p=i; p<M; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULAN_4XCF32(x,y,s);
            BBE_SVN_4XCF32T_XP(x, A_w, N*L*sz_f32c,btol);
          }
         } /* j */
        __Pragma("no_reorder");
      }
#else
      {
        const float32_t tol = FLT_MIN;
        complex_float f,g,s;
        float32_t h,t,r,w;
        int m;
        const complex_float c0f = _makecomplexf(0.f,0.f);
        const complex_float c1f = _makecomplexf(1.f,0.f);
        for ( m=k; m<k+BBE_SIMD_WIDTH/4; m++ ) {
          for ( t=0.f, j=i; j<M; j++ ) {
            t += cabs2f(A[SIDX(j,i,m)]);
          }
          /* Use inverted condition for robust NaN propagation. */
          if (!(t<tol)) {
            r = sqrtf(t);
            f = A[SIDX(i,i,m)]; w = _cabsf(f);
            g = (w>0 ? rcmulf(-r, crdivf(f,w)) : _makecomplexf(-r,0) );
            h = -t-w*r; A[SIDX(i,i,m)] = csubf(f,g);
            for ( j=l; j<N; j++ ) {
              for ( s=c0f, p=i; p<M; p++ ) {
                s = caddf(s, cmuljf(A[SIDX(p,j,m)], A[SIDX(p,i,m)]));
              }
              s = crdivf(s,h);
              for ( p=i; p<M; p++ ) {
                A[SIDX(p,j,m)] = caddf(A[SIDX(p,j,m)], cmulf(s, A[SIDX(p,i,m)]));
              }
            }
          } else {
            g = c0f;
          }
          D[SIDX(0,i,m)] = g;
        } /* m */
      }
#endif
      /*======================================================================*
       *                      Right-hand transformation                        *
       *======================================================================*/
#if 1
      if (l<N-1) {
        xb_vecN_2xf32 r,t,w;
        /* for ( t=0.f, j=l; j<N; j++ ) {
         *   t += cabs2f(A[SIDX(i,j,k)]);
         * } */
        A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,l,k)];
        __Pragma("loop_count min=1");
        for ( t=c0f, j=l; j<N; j++ ) {
          BBE_LVN_4XCF32_XP(f, A_r0, L*sz_f32c);
          r = BBE_MOVN_2XF32_FROMN_4XCF32(f);
          BBE_MULAN_2XF32(t,r,r);
        }
        t = BBE_ADDN_2XF32(t, BBE_SHFLN_2XF32I(t, BBE_SHFLI_SWAP_2));
        /* Mask out columns whose L2 norm is finite, but so small that
         * it cannot be evaluated to working precision. */
        btol = BBE_MOVN_4_FROMN_2(BBE_ULEN_2XF32(tol,t));
        /*---------------------------------------------------------------------* 
         * Construct a Householder vector                                      */
        /* r = sqrtf(t); */
        r = BBE_NEGN_2XF32(BBE_FSQRTN_2XF32(t)); /* !!!! */
        /* f = A[SIDX(i,l,k)]; w = _cabsf(f); */
        f = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
        w = IT_ABSN_4XCF32(f);
        /* g = (w>0 ? rcmulf(-r, crdivf(f,w)) : _makecomplexf(-r,0) ); */
        u = IT_RCMULN_4XCF32(r, IT_CRDIVN_4XCF32(f,w));
        bgtz = BBE_MOVN_4_FROMN_2(BBE_OLTN_2XF32(c0f, w));
        v = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f,r,BBE_SELI_INTERLEAVE_2_EVEN));
        g = BBE_MOVN_4XCF32T(u,v,bgtz);
        g = BBE_MOVN_4XCF32T(g,c0cf,btol);
        BBE_SVN_4XCF32_X(g, F_rw, SIDX(0,i,k)*sz_f32c);
        /* h = -t-w*r; A[SIDX(i,l,k)] = csubf(f,g); */
        h = BBE_NEGN_2XF32(t); BBE_MULAN_2XF32(h,w,r);
        BBE_SVN_4XCF32_X(BBE_SUBN_4XCF32(f,g), A_rw, SIDX(i,l,k)*sz_f32c);
        /*---------------------------------------------------------------------* 
         * Apply the Householder reflector from the right to rows i+1..M-1,    *
         * columns i+1..N-1                                                    */
        for ( j=l; j<M; j++ ) {
          /* for ( s=c0f, p=l; p<N; p++ ) {
           *   s = caddf(s, cmuljf(A[SIDX(j,p,k)], A[SIDX(i,p,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(j,l,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(i,l,k)];
          for ( f=g=u=v=c0cf, p=0; p<(N-l)/2; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, L*sz_f32c);
            BBE_MULMASN_4XCF32(f,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(g,x,y,0x2,0x7);
            BBE_LVN_4XCF32_XP(x, A_r0, L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, L*sz_f32c);
            BBE_MULMASN_4XCF32(u,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(v,x,y,0x2,0x7);
          }
          if ((N-l)&1) {
            x = BBE_LVN_4XCF32_I(A_r0,0);
            y = BBE_LVN_4XCF32_I(A_r1,0);
            BBE_MULMASN_4XCF32(f,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(g,x,y,0x2,0x7);
          }
          s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(f,g), BBE_ADDN_4XCF32(u,v));
          /* s = crdivf(s,h); */
          s = IT_CRDIVN_4XCF32(s,h);
          /* for ( p=l; p<N; p++ ) {
           *   A[SIDX(j,p,k)] = caddf(A[SIDX(j,p,k)], cmulf(s, A[SIDX(i,p,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(j,l,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(i,l,k)];
          A_w  = (xb_vecN_4xcf32*)&A[SIDX(j,l,k)];
          __Pragma("loop_count min=1");
          __Pragma("super_swp ii=12, unroll=4");
          for ( p=l; p<N; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, L*sz_f32c);
            BBE_MULAN_4XCF32(x,y,s);
            BBE_SVN_4XCF32T_XP(x, A_w, L*sz_f32c,btol);
          }
        } /* j */
      } else if (l<N) {
        xb_vecN_4xcf32 g;
        /* F[SIDX(0,i,k)] = A[SIDX(i,l,k)]; */
        g = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
        BBE_SVN_4XCF32_X(g, F_rw, SIDX(0,i,k)*sz_f32c);
      }
      __Pragma("no_reorder");
#else
      {
        const float32_t tol = FLT_MIN;
        complex_float f,g,s;
        float32_t h,t,r,w;
        int m;
        const complex_float c0f = _makecomplexf(0.f,0.f);
        const complex_float c1f = _makecomplexf(1.f,0.f);
        for ( m=k; m<k+BBE_SIMD_WIDTH/4; m++ ) {
          /* Right-hand transformation */
          if (l<N-1) {
            for ( t=0.f, j=l; j<N; j++ ) {
              t += cabs2f(A[SIDX(i,j,m)]);
            }
            if (!(t<tol)) {
              r = sqrtf(t);
              f = A[SIDX(i,l,m)]; w = _cabsf(f);
              g = (w>0 ? rcmulf(-r, crdivf(f,w)) : _makecomplexf(-r,0) );
              h = -t-w*r; A[SIDX(i,l,m)] = csubf(f,g);
              for ( j=l; j<M; j++ ) {
                for ( s=c0f, p=l; p<N; p++ ) {
                  s = caddf(s, cmuljf(A[SIDX(j,p,m)], A[SIDX(i,p,m)]));
                }
                s = crdivf(s,h);
                for ( p=l; p<N; p++ ) {
                  A[SIDX(j,p,m)] = caddf(A[SIDX(j,p,m)], cmulf(s, A[SIDX(i,p,m)]));
                }
              }
            } else {
              g = c0f;
            }
            F[SIDX(0,i,m)] = g;
          } else if (l<N) {
            F[SIDX(0,i,m)] = A[SIDX(i,l,m)];
          } /* l */
        } /* m */
      }
#endif
    } /* i */
    /*
     * If needed, accumulate right-hand transformations from reflection axes 
     * stored in  the work matrix A.
     */
#if 1
    if (V) {
      xb_vecN_4xcf32 r,t;
      /* V[SIDX(N-2,N-2,k)] = V[SIDX(N-1,N-1,k)] = c1f; */
      BBE_SVN_4XCF32_X(c1cf, V_rw, SIDX(N-2,N-2,k)*sz_f32c);
      BBE_SVN_4XCF32_X(c1cf, V_rw, SIDX(N-1,N-1,k)*sz_f32c);
      for ( i=N-3; i>=0; i-- ) {
        /* l = i+1; g = F[SIDX(0,i,k)]; */
        l = i+1;
        g = BBE_LVN_4XCF32_X(F_rw, SIDX(0,i,k)*sz_f32c);
        /* bgtz = (cmagf(g)>0); */
        h = BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(g));
        h = BBE_ADDN_2XF32(h, BBE_SHFLN_2XF32I(h, BBE_SHFLI_SWAP_2));
        bgtz = BBE_MOVN_4_FROMN_2(BBE_OLTN_2XF32(c0f,h));
        /* f = A[SIDX(i,l,k)];  */
        f = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
        /* h = crealf(f)*crealf(g) + cimagf(f)*cimagf(g); */
        h = BBE_MULN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(f),
                           BBE_MOVN_2XF32_FROMN_4XCF32(g));
        h = BBE_ADDN_2XF32(h, BBE_SHFLN_2XF32I(h, BBE_SHFLI_SWAP_2));
        for ( j=l+1; j<N; j++ ) {
          /* for ( s=c0f, p=l+1; p<N; p++ ) {
           *   s = caddf(s, cmulf(A[SIDX(i,p,k)], V[SIDX(p,j,k)]));
           * } */
          V_r  = (xb_vecN_4xcf32*)&V[SIDX(l+1,j,k)];
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,l+1,k)];
          for ( r=t=u=v=c0cf, p=0; p<(N-(l+1))/2; p++ ) {
            BBE_LVN_4XCF32_XP(x, V_r, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r0, L*sz_f32c);
            BBE_MULMASN_4XCF32(r,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(t,x,y,0x1,0x7);
            BBE_LVN_4XCF32_XP(x, V_r, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r0, L*sz_f32c);
            BBE_MULMASN_4XCF32(u,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(v,x,y,0x1,0x7);
          }
          if ((N-(l+1))&1) {
            x = BBE_LVN_4XCF32_I(V_r,0);
            y = BBE_LVN_4XCF32_I(A_r0,0);
            BBE_MULMASN_4XCF32(r,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(t,x,y,0x1,0x7);
          }
          s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(r,t), BBE_ADDN_4XCF32(u,v));
          /* s = crdivf(s,h); */
          s = IT_CRDIVN_4XCF32(s,h);
          /* V[SIDX(l,j,k)] = cmuljf(s, A[SIDX(i,l,k)]); */
          y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
          x = c0cf; BBE_MULJN_4XCF32T(x,s,y,bgtz);
          BBE_SVN_4XCF32_X(x, V_rw, SIDX(l,j,k)*sz_f32c);
          /* for ( p=l+1; p<N; p++ ) {
           *   V[SIDX(p,j,k)] = caddf(V[SIDX(p,j,k)], cmuljf(s, A[SIDX(i,p,k)]));
           * } */
          V_r  = (xb_vecN_4xcf32*)&V[SIDX(l+1,j,k)];
          V_w0 = (xb_vecN_4xcf32*)&V[SIDX(l+1,j,k)];
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,l+1,k)];
          __Pragma("super_swp ii=12, unroll=4");
          __Pragma("loop_count min=1");
          for ( p=l+1; p<N; p++ ) {
            BBE_LVN_4XCF32_XP(x, V_r, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r0, L*sz_f32c);
            BBE_MULJAN_4XCF32T(x,s,y,bgtz);
            BBE_SVN_4XCF32_XP(x, V_w0, N*L*sz_f32c);
          }
        } /* j */
        /* s = crdivf(A[SIDX(i,l,k)], h); */
        y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
        s = IT_CRDIVN_4XCF32(y,h);
        /* V[SIDX(l,l,k)] = caddf(c1f, cmuljf(s,A[SIDX(i,l,k)])); */
        y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,l,k)*sz_f32c);
        x = c1cf; BBE_MULJAN_4XCF32T(x,s,y,bgtz);
        BBE_SVN_4XCF32_X(x, V_rw, SIDX(l,l,k)*sz_f32c);
        /* for ( j=l+1; j<N; j++ ) {
         *   V[SIDX(j,l,k)] = cmuljf(s,A[SIDX(i,j,k)]);
         * } */
        V_w0 = (xb_vecN_4xcf32*)&V[SIDX(l+1,l,k)];
        A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i,l+1,k)];
        __Pragma("loop_count min=1");
        for ( j=l+1; j<N; j++ ) {
          BBE_LVN_4XCF32_XP(y, A_r0, L*sz_f32c);
          x = BBE_MULJN_4XCF32(s,y);
          x = BBE_MOVN_4XCF32T(x,c0cf,bgtz);
          BBE_SVN_4XCF32_XP(x, V_w0, N*L*sz_f32c);
        }
        __Pragma("no_reorder");
      } /* i */
      /* V[SIDX(0,0,k)] = c1f; */
      /* for ( j=1; j<N; j++ ) {
       *   V[SIDX(0,j,k)] = V[SIDX(j,0,k)] = c0f;
       * } */
      BBE_SVN_4XCF32_X(c1cf, V_rw, SIDX(0,0,k)*sz_f32c);
      V_w0 = (xb_vecN_4xcf32*)&V[SIDX(1,0,k)];
      V_w1 = (xb_vecN_4xcf32*)&V[SIDX(0,1,k)];
      __Pragma("loop_count min=1");
      for ( i=1; i<N; i++ ) {
        BBE_SVN_4XCF32_XP(c0cf, V_w0, N*L*sz_f32c);
        BBE_SVN_4XCF32_XP(c0cf, V_w1, L*sz_f32c);
      }
    } /* B */
#else
    if (V) {
      complex_float f,g,s;
      float32_t h;
      int m;
      const complex_float c0f = _makecomplexf(0.f,0.f);
      const complex_float c1f = _makecomplexf(1.f,0.f);
      for ( m=k; m<k+BBE_SIMD_WIDTH/4; m++ ) {
        V[SIDX(N-2,N-2,m)] = V[SIDX(N-1,N-1,m)] = c1f;
        for ( i=N-3; i>=0; i-- ) {
          l = i+1; g = F[SIDX(0,i,m)];
          if (cmagf(g)>0) {
            f = A[SIDX(i,l,m)]; 
            h = crealf(f)*crealf(g) + cimagf(f)*cimagf(g);
            for ( j=l+1; j<N; j++ ) {
              for ( s=c0f, p=l+1; p<N; p++ ) {
                s = caddf(s, cmulf(A[SIDX(i,p,m)], V[SIDX(p,j,m)]));
              }
              s = crdivf(s,h);
              V[SIDX(l,j,m)] = cmuljf(s, A[SIDX(i,l,m)]);
              for ( p=l+1; p<N; p++ ) {
                V[SIDX(p,j,m)] = caddf(V[SIDX(p,j,m)], cmuljf(s, A[SIDX(i,p,m)]));
              }
            }
            s = crdivf(A[SIDX(i,l,m)], h);
            for ( j=l; j<N; j++ ) {
              V[SIDX(j,l,m)] = cmuljf(s,A[SIDX(i,j,m)]);
            }
          } else {
            V[SIDX(l,l,m)] = c0f;
            for ( j=l+1; j<N; j++ ) {
              V[SIDX(l,j,m)] = V[SIDX(j,l,m)] = c0f;
            }
          }
          V[SIDX(l,l,m)] = caddf(V[SIDX(l,l,m)], c1f);
        } /* i */
        V[SIDX(0,0,m)] = c1f;
        for ( j=1; j<N; j++ ) {
          V[SIDX(0,j,m)] = V[SIDX(j,0,m)] = c0f;
        }
      } /* m */
    } /* V */
#endif
    /*
     * If needed, accumulate the left-hand transformation. Matrix U overrides 
     * intermediate data stored in A. 
     */
#if 1
    if (needU) {
      xb_vecN_4xcf32 r,t;
      for ( i=N-1; i>=0; i-- ) {
        l = i+1; 
        /* g = D[SIDX(0,i,k)]; */
        g = BBE_LVN_4XCF32_X(D_rw, SIDX(0,i,k)*sz_f32c);
        /* bgtz = (cmagf(g)>0); */
        h = BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(g));
        h = BBE_ADDN_2XF32(h, BBE_SHFLN_2XF32I(h, BBE_SHFLI_SWAP_2));
        bgtz = BBE_MOVN_4_FROMN_2(BBE_OLTN_2XF32(c0f,h));
        /* f = A[SIDX(i,i,k)]; */
        f = BBE_LVN_4XCF32_X(A_rw, SIDX(i,i,k)*sz_f32c);
        /* h = crealf(f)*crealf(g) + cimagf(f)*cimagf(g); */
        h = BBE_MULN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(f),
                           BBE_MOVN_2XF32_FROMN_4XCF32(g));
        h = BBE_ADDN_2XF32(h, BBE_SHFLN_2XF32I(h, BBE_SHFLI_SWAP_2));
        for ( j=l; j<N; j++ ) {
          /* for ( s=c0f, p=l; p<M; p++ ) {
           *   s = caddf(s, cmuljf(A[SIDX(p,j,k)], A[SIDX(p,i,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(l,j,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(l,i,k)];
          for ( r=t=u=v=c0cf, p=0; p<(M-l)/2; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULMASN_4XCF32(r,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(t,x,y,0x2,0x7);
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULMASN_4XCF32(u,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(v,x,y,0x2,0x7);
          }
          if ((M-l)&1) {
            x = BBE_LVN_4XCF32_I(A_r0,0);
            y = BBE_LVN_4XCF32_I(A_r1,0);
            BBE_MULMASN_4XCF32(r,x,y,0x0,0x8);
            BBE_MULMASN_4XCF32(t,x,y,0x2,0x7);
          }
          s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(r,t), BBE_ADDN_4XCF32(u,v));
          /* s = crdivf(s,h); */
          s = IT_CRDIVN_4XCF32(s,h);
          /* A[SIDX(i,j,k)] = cmulf(s, A[SIDX(i,i,k)]); */
          y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,i,k)*sz_f32c);
          x = c0cf; BBE_MULN_4XCF32T(x,s,y,bgtz);
          BBE_SVN_4XCF32_X(x, A_rw, SIDX(i,j,k)*sz_f32c);
          /* for ( p=l; p<M; p++ ) {
           *   A[SIDX(p,j,k)] = caddf(A[SIDX(p,j,k)], cmulf(s, A[SIDX(p,i,k)]));
           * } */
          A_r0 = (xb_vecN_4xcf32*)&A[SIDX(l,j,k)];
          A_r1 = (xb_vecN_4xcf32*)&A[SIDX(l,i,k)];
          A_w  = (xb_vecN_4xcf32*)&A[SIDX(l,j,k)];
          __Pragma("super_swp ii=12, unroll=4");
          for ( p=l; p<M; p++ ) {
            BBE_LVN_4XCF32_XP(x, A_r0, N*L*sz_f32c);
            BBE_LVN_4XCF32_XP(y, A_r1, N*L*sz_f32c);
            BBE_MULAN_4XCF32T(x,s,y,bgtz);
            BBE_SVN_4XCF32_XP(x, A_w, N*L*sz_f32c);
          }
        } /* j */
        /* s = _conjf(crdivf(A[SIDX(i,i,k)], h)); */
        y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,i,k)*sz_f32c);
        s = BBE_CONJN_4XCF32(IT_CRDIVN_4XCF32(y,h));
        /* A[SIDX(i,i,k)] = caddf(c1f, cmulf(s, A[SIDX(i,i,k)])); */
        y = BBE_LVN_4XCF32_X(A_rw, SIDX(i,i,k)*sz_f32c);
        x = c1cf; BBE_MULAN_4XCF32T(x,s,y,bgtz);
        BBE_SVN_4XCF32_X(x, A_rw, SIDX(i,i,k)*sz_f32c);
        /* for ( j=i+1; j<M; j++ ) {
         *   A[SIDX(j,i,k)] = cmulf(s, A[SIDX(j,i,k)]);
         * } */
        A_r0 = (xb_vecN_4xcf32*)&A[SIDX(i+1,i,k)];
        A_w  = (xb_vecN_4xcf32*)&A[SIDX(i+1,i,k)];
        for ( j=i+1; j<M; j++ ) {
          BBE_LVN_4XCF32_XP(y, A_r0, N*L*sz_f32c);
          x = BBE_MULN_4XCF32(s,y);
          x = BBE_MOVN_4XCF32T(x,c0cf,bgtz);
          BBE_SVN_4XCF32_XP(x, A_w, N*L*sz_f32c);
        }
      } /* i */
    } /* needU */
#else
    if (needU) {
      complex_float f,g,s;
      float32_t h;
      int m;
      const complex_float c0f = _makecomplexf(0.f,0.f);
      const complex_float c1f = _makecomplexf(1.f,0.f);
      for ( m=k; m<k+BBE_SIMD_WIDTH/4; m++ ) {
        for ( i=N-1; i>=0; i-- ) {
          l = i+1; g = D[SIDX(0,i,m)];
          if (cmagf(g)>0) {
            f = A[SIDX(i,i,m)];
            h = crealf(f)*crealf(g) + cimagf(f)*cimagf(g);
            for ( j=l; j<N; j++ ) {
              for ( s=c0f, p=l; p<M; p++ ) {
                s = caddf(s, cmuljf(A[SIDX(p,j,m)], A[SIDX(p,i,m)]));
              }
              s = crdivf(s,h);
              A[SIDX(i,j,m)] = cmulf(s, A[SIDX(i,i,m)]);
              for ( p=l; p<M; p++ ) {
                A[SIDX(p,j,m)] = caddf(A[SIDX(p,j,m)], cmulf(s, A[SIDX(p,i,m)]));
              }
            }
            s = _conjf(crdivf(A[SIDX(i,i,m)], h));
            A[SIDX(i,i,m)] = caddf(c1f, cmulf(s, A[SIDX(i,i,m)]));
            for ( j=i+1; j<M; j++ ) {
              A[SIDX(j,i,m)] = cmulf(s, A[SIDX(j,i,m)]);
            }
          } else {
            A[SIDX(i,i,m)] = c1f;;
            for ( j=i+1; j<M; j++ ) {
              A[SIDX(j,i,m)] = c0f;
              A[SIDX(i,j,m)] = c0f;
            }
          }
        } /* i */
      } /* m */
    } /* needU */
#endif
  } /* k */

  #undef SIDX

} /* svd_bidiag_mxnsf() */

#endif /* HAVE_VFPU */

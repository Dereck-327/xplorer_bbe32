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
    Hessenberg Form Of Complex Matrices, Block and Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <float.h>
#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Fixed-point arithmetics. */
//#include "NatureDSP_Math.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#include "vfpu_math.h"

#if HAVE_VFPU

#define sz_f32c   sizeof(complex_float)

#if 0
#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

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
#endif

#if 0
/* Reduction of a single matrix to upper-Hessenberg form. Reduced matrix overwrites 
 * input data. Both the block and stream order are supported.
 * This function is based on orthes() algorithm from
 *   "Similarity Reduction of a General Matrix to Hessenberg Form" by R.S. Martin and
 *   J.H. Wilkinson, Handbook for Automatie Computation, Vol.II Linear Algebra, 
 *   Contribution II/13. */
static void chessf ( complex_float * restrict P, 
                     complex_float * restrict A, 
                     int N, int stride )
{
  /*
   * MATLAB reference code and a simple test:
   *
   *   n = 5; A = complex(1-2*rand(n),1-2*rand(n)); B = A;
   *   % Reduce the matrix to upper-Hessenberg form
   *   P = eye(n);
   *   tol = realmin(class(A));
   *   for m=2:n-1
   *     u = A(m:n,m-1);
   *     h = u'*u;
   *     if h>tol
   *       % Construct and apply a Householder reflector.
   *       f = sqrt(h);
   *       % Note that u(1) may be subnormal or even zero!
   *       t = cabs(u(1)); 
   *       if t>tol, g = -f*(u(1)/t); else g = -f; end;
   *       h = h+f*t;
   *       u(1) = u(1)-g;
   *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
   *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
   *       % Accumulate transformations
   *       P(2:n,m:n) = P(2:n,m:n)-(P(2:n,m:n)*u/h)*u';
   *     else
   *       g = 0;
   *     end
   *     A(m,m-1) = g;
   *     A(m+1:n,m-1) = 0;
   *   end
   *   % Validate the result
   *   R = B-P*A*P';
   *   r = sqrt(R(:)'*R(:));
   *   fprintf('r = %.1e\n',r);
   */

  #define IDX(i,j)  (((i)*N+(j))*stride)

  float32_t f,h,t;
  complex_float g,s;
  int m,n,p;
  const complex_float c0f = _makecomplexf(0.f,0.f);
  /* Initialize P with identity matrix. */
  if (P) {
    for ( m=0; m<N; m++ ) {
      for ( n=0; n<N; n++ ) {
        P[IDX(m,n)] = _makecomplexf( m==n ? 1.f : 0.f, 0.f );
      }
    }
  }
  /* Perform N-2 Householder reflections. */
  for ( m=1; m<N-1; m++ ) {
    /* Squared norm of the pivot column */
    for ( h=0.f, n=m; n<N; n++ ) h += cabs2f(A[IDX(n,m-1)]);
    /* Skip the transformation if the norm is too small to 
     * be evaluated to working precision. */
    if (h<=FLT_MIN) { A[IDX(m,m-1)] = c0f; continue; }
    /* In-place construction of a Householder vector. Note that
     * A(m,m-1) may be subnormal or even zero. */
    f = sqrtf(h); t = _cabsf(A[IDX(m,m-1)]);
    if (t>FLT_MIN) {
      g = rcmulf(-f, crdivf(A[IDX(m,m-1)],t));
    } else {
      g = _makecomplexf(-f,0.f);
    }
    h += f*t;
    A[IDX(m,m-1)] = csubf(A[IDX(m,m-1)], g);
    /* Left-hand transformation */
    for ( n=m; n<N; n++ ) {
      for ( s=c0f, p=m; p<N; p++ ) {
        s = caddf(s, cmuljf(A[IDX(p,n)], A[IDX(p,m-1)]));
      }
      s = crdivf(s,h);
      for ( p=m; p<N; p++ ) {
        A[IDX(p,n)] = csubf(A[IDX(p,n)], cmulf(s, A[IDX(p,m-1)]));
      }
    }
    /* Right-hand transformation */
    for ( n=0; n<N; n++ ) {
      for ( s=c0f, p=m; p<N; p++ ) {
        s = caddf(s, cmulf(A[IDX(n,p)], A[IDX(p,m-1)]));
      }
      s = crdivf(s,h);
      for ( p=m; p<N; p++ ) {
        A[IDX(n,p)] = csubf(A[IDX(n,p)], cmuljf(s, A[IDX(p,m-1)]));
      }
    }
    /* Accumulate transformations */
    if (P) {
      for ( n=0; n<N; n++ ) {
        for ( s=c0f, p=m; p<N; p++ ) {
          s = caddf(s, cmulf(P[IDX(n,p)], A[IDX(p,m-1)]));
        }
        s = crdivf(s,h);
        for ( p=m; p<N; p++ ) {
          P[IDX(n,p)] = csubf(P[IDX(n,p)], cmuljf(s, A[IDX(p,m-1)]));
        }
      }
    }
    A[IDX(m,m-1)] = g;
  } /* m */

  #undef IDX

} /* chessf() */
#endif

/*
 * Reduce square input matrices to upper-Hessenberg form by a unitary 
 * similarity transformation: H <- P'*A*P. Resulting matrices H replace
 * input matrices A.
 * Input:
 *   N    Matrix size
 *   L    Number of matrices
 * Input/Output:
 *   A[]  General square matrices (in), reduced matrices (out). 
 *        Values below the first subdiagonal of resulting 
 *        matrices are not defined
 * Output:
 *   P[]  Unitary/orthogonal trasformation matrix, This argument is 
 *        OPTIONAL, set it to zero if matrix P is not required
 * Restrictions:
 *   A,P  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *        boundary
 *   N>2  Minimum maxtrix size is 3x3
 *   N,L  Subject to additional limitations exposed by a particular function
 */

/* Complex-valued, Stream Order
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void eigen_hess_nxnsf ( 
                complex_float * restrict P, /* P[N*N][L] */
                complex_float * restrict A, /* A[N*N][L] */
                int N, int L )
#if 0
{
  int k;
  NASSERT_ALIGN( P, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  for ( k=0; k<L; k++ ) {
    chessf( ( P ? P+k : 0 ), A+k, N, L );
  }
} /* eigen_hess_nxnsf() */
#elif 0
{
    /*
    * MATLAB reference code and a simple test:
    *
    *   n = 5; A = complex(1-2*rand(n),1-2*rand(n)); B = A;
    *   % Reduce the matrix to upper-Hessenberg form
    *   P = eye(n);
    *   tol = realmin(class(A));
    *   for m=2:n-1
    *     u = A(m:n,m-1);
    *     h = u'*u;
    *     if h>tol
    *       % Construct and apply a Householder reflector.
    *       f = sqrt(h);
    *       % Note that u(1) may be subnormal or even zero!
    *       t = cabs(u(1));
    *       if t>tol, g = -f*(u(1)/t); else g = -f; end;
    *       h = h+f*t;
    *       u(1) = u(1)-g;
    *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
    *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
    *       % Accumulate transformations
    *       P(2:n,m:n) = P(2:n,m:n)-(P(2:n,m:n)*u/h)*u';
    *     else
    *       g = 0;
    *     end
    *     A(m,m-1) = g;
    *     A(m+1:n,m-1) = 0;
    *   end
    *   % Validate the result
    *   R = B-P*A*P';
    *   r = sqrt(R(:)'*R(:));
    *   fprintf('r = %.1e\n',r);
    */
#define IDX(i,j)  (((i)*N+(j))*stride)
    const complex_float c0f = _makecomplexf(0.f, 0.f);
    const int stride = L;

    float32_t f, h, t;
    complex_float g, s;
    int k, m, n, p;

    NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));

    for (k = 0; k<L; k++) {
        /* Initialize P with identity matrix. */
        if (P) {
            for (m = 0; m<N; m++) {
                for (n = 0; n<N; n++) {
                    P[IDX(m, n)] = _makecomplexf(m == n ? 1.f : 0.f, 0.f);
                }
            }
        }
        /* Perform N-2 Householder reflections. */
        for (m = 1; m<N - 1; m++) {
            /* Squared norm of the pivot column */
            for (h = 0.f, n = m; n<N; n++) h += cabs2f(A[IDX(n, m - 1)]);
            /* Skip the transformation if the norm is too small to
            * be evaluated to working precision. */
            if (h <= FLT_MIN) { A[IDX(m, m - 1)] = c0f; continue; }
            /* In-place construction of a Householder vector. Note that
            * A(m,m-1) may be subnormal or even zero. */
            f = sqrtf(h); t = _cabsf(A[IDX(m, m - 1)]);
            if (t>FLT_MIN) {
                g = rcmulf(-f, crdivf(A[IDX(m, m - 1)], t));
            } else {
                g = _makecomplexf(-f, 0.f);
            }
            h += f*t;
            A[IDX(m, m - 1)] = csubf(A[IDX(m, m - 1)], g);
            /* Left-hand transformation */
            for (n = m; n<N; n++) {
                for (s = c0f, p = m; p<N; p++) {
                    s = caddf(s, cmuljf(A[IDX(p, n)], A[IDX(p, m - 1)]));
                }
                s = crdivf(s, h);
                for (p = m; p<N; p++) {
                    A[IDX(p, n)] = csubf(A[IDX(p, n)], cmulf(s, A[IDX(p, m - 1)]));
                }
            }
            /* Right-hand transformation */
            for (n = 0; n<N; n++) {
                for (s = c0f, p = m; p<N; p++) {
                    s = caddf(s, cmulf(A[IDX(n, p)], A[IDX(p, m - 1)]));
                }
                s = crdivf(s, h);
                for (p = m; p<N; p++) {
                    A[IDX(n, p)] = csubf(A[IDX(n, p)], cmuljf(s, A[IDX(p, m - 1)]));
                }
            }
            /* Accumulate transformations */
            if (P) {
                for (n = 0; n<N; n++) {
                    for (s = c0f, p = m; p<N; p++) {
                        s = caddf(s, cmulf(P[IDX(n, p)], A[IDX(p, m - 1)]));
                    }
                    s = crdivf(s, h);
                    for (p = m; p<N; p++) {
                        P[IDX(n, p)] = csubf(P[IDX(n, p)], cmuljf(s, A[IDX(p, m - 1)]));
                    }
                }
            }
            A[IDX(m, m - 1)] = g;
        } /* m */

        A++;
        if (P) P ++;
    }
} /* eigen_hess_nxnsf() */

#undef IDX
#else
{
    /*
    * MATLAB reference code and a simple test:
    *
    *   n = 5; A = complex(1-2*rand(n),1-2*rand(n)); B = A;
    *   % Reduce the matrix to upper-Hessenberg form
    *   P = eye(n);
    *   tol = realmin(class(A));
    *   for m=2:n-1
    *     u = A(m:n,m-1);
    *     h = u'*u;
    *     if h>tol
    *       % Construct and apply a Householder reflector.
    *       f = sqrt(h);
    *       % Note that u(1) may be subnormal or even zero!
    *       t = cabs(u(1));
    *       if t>tol, g = -f*(u(1)/t); else g = -f; end;
    *       h = h+f*t;
    *       u(1) = u(1)-g;
    *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
    *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
    *       % Accumulate transformations
    *       P(2:n,m:n) = P(2:n,m:n)-(P(2:n,m:n)*u/h)*u';
    *     else
    *       g = 0;
    *     end
    *     A(m,m-1) = g;
    *     A(m+1:n,m-1) = 0;
    *   end
    *   % Validate the result
    *   R = B-P*A*P';
    *   r = sqrt(R(:)'*R(:));
    *   fprintf('r = %.1e\n',r);
    */
#define IDX(i,j)  (((i)*N+(j))*stride)
    const xb_vecN_4xcf32 * restrict A_r0;
    const xb_vecN_4xcf32 * restrict A_r1;
    const xb_vecN_4xcf32 * restrict P_r;
          xb_vecN_4xcf32 * restrict A_w;
          xb_vecN_4xcf32 * restrict P_w;

    const int stride = L;
    const xb_vecN_4xcf32 _c0f = BBE_CONSTN_4XCF32(0);
    const xb_vecN_4xcf32 _c1f = BBE_SELN_4XCF32I(_c0f, BBE_CONSTN_4XCF32(1), BBE_SELI_INTERLEAVE_2_LO);
    const xb_vecN_2xf32 F_MIN = FLT_MIN;

    xb_vecN_2xf32 _f, _fneg, _h, _hrecip, _t, x;
    xb_vecN_4xcf32 a, a0, a1, _g;
    xb_vecN_4xcf32 s, s0, s1, s2, s3;
    vboolN_2 bskip_2;
    vboolN_4 bskip_4, bless_4;

    int k, m, n, p;
    
    NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));

    if (P){
        /* Initialize P with identity matrix. */
        const int step_norm = (N + 1)*L*sz_f32c;
        const int step_back = (BBE_SIMD_WIDTH / 4) * sz_f32c - (N - 1)*(N + 1)*L*sz_f32c;
        int step;
        P_w = (xb_vecN_4xcf32*)P;
        for (n = 0; n < (N*N*L / (BBE_SIMD_WIDTH / 4)); n++) {
            BBE_SVN_4XCF32_IP(_c0f, P_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        }
        P_w = (xb_vecN_4xcf32*)P;
        for (n = 0, k = 0; k < (L / (BBE_SIMD_WIDTH / 4))*N; k++) {
            n = BBE_ADDMOD16U(n, (N << 16) | 1);
            step = step_norm;
            XT_MOVEQZ(step, step_back, n);
            BBE_SVN_4XCF32_XP(_c1f, P_w, step);
        }
    }

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++) {
        /* Perform N-2 Householder reflections. */
        for (m = 1; m < N - 1; m++) {
            _h = BBE_ZERON_2XF32();
            /* Squared norm of the pivot column */
            A_r0 = (xb_vecN_4xcf32 *)&A[(m*N + (m-1))*L];
            for (n = m; n < N; n++){
                BBE_LVN_4XCF32_XP(a, A_r0, N*L*sz_f32c);
                x = BBE_MOVN_2XF32_FROMN_4XCF32(a);
                BBE_MULAN_2XF32(_h, x, x);
            }
            _h = BBE_ADDN_2XF32(_h, BBE_SHFLN_2XF32I(_h, BBE_SHFLI_SWAP_2));
            /* Skip the transformation if the norm is too small to
            * be evaluated to working precision. */
            A_r0 = A_w = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
            bskip_2 = BBE_OLEN_2XF32(_h, F_MIN);
            bskip_4 = BBE_MOVN_4_FROMN_2(bskip_2);
            BBE_SVN_4XCF32T_I(_c0f, A_w, 0, bskip_4);
            /* In-place construction of a Householder vector. Note that
            * A(m,m-1) may be subnormal or even zero. */
            //f = sqrtf(h);
            _f = BBE_SQRTN_2XF32(_h);
            _fneg = BBE_NEGN_2XF32(_f);
            // t = _cabsf(A[IDX(m, m - 1)]);
            a = BBE_LVN_4XCF32_I(A_r0, 0);
            _t = IT_ABSN_4XCF32(a);
            bless_4 = BBE_MOVN_4_FROMN_2(BBE_OLTN_2XF32(F_MIN, _t));
            //g = rcmulf(-f, crdivf(A[IDX(m, m - 1)], t));
            a = IT_CRDIVN_4XCF32(a, _t);
            x = BBE_MULN_2XF32(_fneg, BBE_MOVN_2XF32_FROMN_4XCF32(a));
            _g = BBE_MOVN_4XCF32_FROMN_2XF32(x);
            //g = _makecomplexf(-f, 0.f);
            x = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0), _fneg, BBE_SELI_INTERLEAVE_2_EVEN);
            a = BBE_MOVN_4XCF32_FROMN_2XF32(x);
            _g = BBE_MOVN_4XCF32T(_g, a, bless_4);
            //h += f*t;
            BBE_MULAN_2XF32(_h, _f, _t);
            //A[IDX(m, m - 1)] = csubf(A[IDX(m, m - 1)], g);
            a = BBE_LVN_4XCF32F_I(A_r0, 0, bskip_4);
            a = BBE_SUBN_4XCF32(a,_g);
            BBE_SVN_4XCF32F_I(a, A_w, 0, bskip_4);
            // Calculate reciprocal
            _hrecip = BBE_RECIP0N_2XF32(_h);
            _t = BBE_CONSTN_2XF32(1);
            BBE_MULSN_2XF32(_t, _hrecip, _h);
            BBE_MULAN_2XF32(_hrecip, _t, _hrecip);

            /* Left-hand transformation */
            for (n = m; n < N; n++) {
                A_r0 = (xb_vecN_4xcf32 *)&A[IDX(m, n)];
                A_r1 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
#if 0
                s = BBE_ZERON_4XCF32();
                for (p = m; p < N; p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0, N*L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULJAN_4XCF32(s, a0, a1);
                }
#else
                s0 = s1 = s2 = s3 = BBE_ZERON_4XCF32();
                for (p = 0; p < ((N - m)>>1); p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0, N*L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s0, a0, a1, 2, 4);
                    BBE_MULMASN_4XCF32(s1, a0, a1, 0, 11);
                    BBE_LVN_4XCF32_XP(a0, A_r0, N*L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s2, a0, a1, 2, 4);
                    BBE_MULMASN_4XCF32(s3, a0, a1, 0, 11);
                }
                if ((N - m) & 1){
                    BBE_LVN_4XCF32_XP(a0, A_r0, N*L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s0, a0, a1, 2, 4);
                    BBE_MULMASN_4XCF32(s1, a0, a1, 0, 11);
                }
                s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(s0,s1), BBE_ADDN_4XCF32(s2,s3));
#endif
#if 0
                s = IT_CRDIVN_4XCF32(s,_h);
#else
                x = BBE_MOVN_2XF32_FROMN_4XCF32(s);
                _t = BBE_MULN_2XF32(_hrecip, x);
                BBE_MULSN_2XF32(x, _h, _t); BBE_MULAN_2XF32(_t, _hrecip, x);
                s = BBE_MOVN_4XCF32_FROMN_2XF32(_t);
#endif
                A_r0 = A_w = (xb_vecN_4xcf32 *)&A[IDX(m, n)];
                A_r1 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
                for (p = m; p < N; p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0, N*L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULSN_4XCF32(a0, s, a1);
                    BBE_SVN_4XCF32F_XP(a0, A_w, N*L*sz_f32c, bskip_4);
                }
                __Pragma("no_reorder");
            }
            /* Right-hand transformation */
            for (n = 0; n < N; n++) {
                A_r0 = (xb_vecN_4xcf32 *)&A[IDX(n, m)];
                A_r1 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
#if 0
                s = BBE_ZERON_4XCF32();
                for (p = m; p < N; p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0,   L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULAN_4XCF32(s, a0, a1);
                }
#else
                s0 = s1 = s2 = s3 = BBE_ZERON_4XCF32();
                for (p = 0; p < (N-m)>>1; p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0, L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s0, a0, a1, 0, 4);
                    BBE_MULMASN_4XCF32(s1, a0, a1, 1, 11);
                    BBE_LVN_4XCF32_XP(a0, A_r0, L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s2, a0, a1, 0, 4);
                    BBE_MULMASN_4XCF32(s3, a0, a1, 1, 11);
                }
                if ((N - m) & 1){
                    BBE_LVN_4XCF32_XP(a0, A_r0, L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULMASN_4XCF32(s0, a0, a1, 0, 4);
                    BBE_MULMASN_4XCF32(s1, a0, a1, 1, 11);
                }
                s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(s0, s1), BBE_ADDN_4XCF32(s2, s3));
#endif
                //s = crdivf(s, h);
#if 0
                s = IT_CRDIVN_4XCF32(s, _h);
#else
                x = BBE_MOVN_2XF32_FROMN_4XCF32(s);
                _t = BBE_MULN_2XF32(_hrecip, x);
                BBE_MULSN_2XF32(x, _h, _t); BBE_MULAN_2XF32(_t, _hrecip, x);
                s = BBE_MOVN_4XCF32_FROMN_2XF32(_t);
#endif
                A_r0 = A_w = (xb_vecN_4xcf32 *)&A[IDX(n, m)];
                A_r1 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
                for (p = m; p < N; p++) {
                    BBE_LVN_4XCF32_XP(a0, A_r0,   L*sz_f32c);
                    BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                    BBE_MULJSN_4XCF32(a0, s, a1);
                    BBE_SVN_4XCF32F_XP(a0, A_w, L*sz_f32c, bskip_4);
                }
                __Pragma("no_reorder");
            }
            /* Accumulate transformations */
            if (P) {
                for (n = 0; n < N; n++) {
                    P_r = (xb_vecN_4xcf32 *)&P[IDX(n, m)];
                    A_r1 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
#if 0
                    s = BBE_ZERON_4XCF32();
                    for (p = m; p < N; p++) {
                        BBE_LVN_4XCF32_XP(a0, P_r,    L*sz_f32c);
                        BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                        BBE_MULAN_4XCF32(s, a0, a1);
                    }
#else
                    s0 = s1 = s2 = s3 = BBE_ZERON_4XCF32();
                    for (p = 0; p < (N-m)>>1; p++) {
                        BBE_LVN_4XCF32_XP(a0, P_r, L*sz_f32c);
                        BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                        BBE_MULMASN_4XCF32(s0, a0, a1, 0, 4);
                        BBE_MULMASN_4XCF32(s1, a0, a1, 1, 11);
                        BBE_LVN_4XCF32_XP(a0, P_r, L*sz_f32c);
                        BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                        BBE_MULMASN_4XCF32(s2, a0, a1, 0, 4);
                        BBE_MULMASN_4XCF32(s3, a0, a1, 1, 11);
                    }
                    if ((N - m) & 1){
                        BBE_LVN_4XCF32_XP(a0, P_r, L*sz_f32c);
                        BBE_LVN_4XCF32_XP(a1, A_r1, N*L*sz_f32c);
                        BBE_MULMASN_4XCF32(s0, a0, a1, 0, 4);
                        BBE_MULMASN_4XCF32(s1, a0, a1, 1, 11);
                    }
                    s = BBE_ADDN_4XCF32(BBE_ADDN_4XCF32(s0, s1), BBE_ADDN_4XCF32(s2, s3));
#endif
                    //s = crdivf(s, h);
#if 0
                    s = IT_CRDIVN_4XCF32(s, _h);
#else
                    x = BBE_MOVN_2XF32_FROMN_4XCF32(s);
                    _t = BBE_MULN_2XF32(_hrecip, x);
                    BBE_MULSN_2XF32(x, _h, _t); BBE_MULAN_2XF32(_t, _hrecip, x);
                    s = BBE_MOVN_4XCF32_FROMN_2XF32(_t);
#endif
                    P_r = P_w = (xb_vecN_4xcf32 *)&P[IDX(n, m)];
                    A_r0 = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
                    for (p = m; p < N; p++) {
                        BBE_LVN_4XCF32_XP(a0, P_r, L*sz_f32c);
                        BBE_LVN_4XCF32_XP(a1, A_r0, N*L*sz_f32c);
                        BBE_MULJSN_4XCF32(a0, s, a1);
                        BBE_SVN_4XCF32F_XP(a0, P_w, L*sz_f32c, bskip_4);
                    }
                }
            }
            //A[IDX(m, m - 1) + k] = g;
            A_w = (xb_vecN_4xcf32 *)&A[IDX(m, m - 1)];
            BBE_SVN_4XCF32F_I(_g, A_w, 0, bskip_4);
        } /* m */
        A += (BBE_SIMD_WIDTH / 4);
        if (P) P += (BBE_SIMD_WIDTH / 4);
    }
} /* eigen_hess_nxnsf() */
#undef IDX
#endif

#endif /* HAVE_VFPU */

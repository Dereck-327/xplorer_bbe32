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
    Francis QR algorithm for complex upper-Hessenberg matrix
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* !!!! *
#include <stdio.h>
* !!!! */

#include <math.h>
#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Real and complex arithmetic primitives optimized for BBEN VFPU */
#include "vfpu_math.h"
/* NaN values for single precision routines */
#include "nanf_tbl.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define sz_f32c   sizeof(complex_float)

#define EPS       FLT_EPSILON
#define REALMAX   FLT_MAX
#define ITS_LIM   40 /* Iterations count limit for the QR algorithm */

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )
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

/* Complex floating-point square root, single precision. Returns the
 * square root that lies in the right half of the complex plane. */
static complex_float _csqrtf( complex_float x )
{
  /*
   * Based on csqrt code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   * MATLAB reference code:
   *   function z = csqrt(x)
   *   xr = real(x); xi = imag(x);
   *   h = sqrt((abs(xr)+cabs(x))/2);
   *   if xi~=0, xi = xi/(2*h); end;
   *   if xr<0
   *     if xi<0
   *       xr = -xi; xi = -h;
   *     else
   *       xr = xi; xi = h;
   *     end
   *   else
   *     xr = h;
   *   end
   *   z = complex(xr,xi);
   */

  float32_t xr,xi,h;
  xr = crealf(x); xi = cimagf(x);
  h = sqrtf((fabsf(xr) + _cabsf(x))/2);
  if (xi!=0) xi /= 2*h;
  if (xr<0) {
    xr = ( xi<0 ? -xi : xi );
    xi = ( xi<0 ? -h : h );
  } else {
    xr = h;
  }
  return ( _makecomplexf(xr,xi) );

} /* _csqrtf() */
#endif
/* !!!! */
#if 0
void print_csqrt( float32_t xr, float32_t xi )
{
  xtcomplexfloat x,y0,y1;
  union ufloat32uint32 _xr,_xi,_y0r,_y0i,_y1r,_y1i;
  x = BBE_CMPLXF32(xi,xr);
  *(complex_float*)&y0 = _csqrtf(*(complex_float*)&x);
  y1 = IT_SQRTCF32(x);
  _xr.f = xr; _xi.f = xi;
  _y0r.f = BBE_CREALCF32(y0); _y0i.f = BBE_CIMAGCF32(y0);
  _y1r.f = BBE_CREALCF32(y1); _y1i.f = BBE_CIMAGCF32(y1);
  printf("x: %+f 0x%08x %+f 0x%08x;\n  y0: %+f 0x%08x %+f 0x%08x;\n  y1: %+f 0x%08x %+f 0x%08x;\n\n",
         (float64_t)_xr.f, _xr.u, (float64_t)_xi.f, _xi.u,
         (float64_t)_y0r.f, _y0r.u, (float64_t)_y0i.f, _y0i.u,
         (float64_t)_y1r.f, _y1r.u, (float64_t)_y1i.f, _y1i.u );
}

int main(void)
{
  print_csqrt(0.f,0.f);
  print_csqrt(1.f,1.f);
  print_csqrt(0.f,1.f);
  print_csqrt(-1.f,1.f);
  print_csqrt(-1.f,0.f);
  print_csqrt(-1.f,-1.f);
  print_csqrt(0.f,-1.f);
  print_csqrt(1.f,-1.f);
  print_csqrt(1e-40f,1e-40f);
  print_csqrt(-1e-40f,0.f);
  print_csqrt(qNaNf.f,0.f);
  print_csqrt(qNaNf.f,1.f);
  print_csqrt(0.f,qNaNf.f);
  print_csqrt(1.f,qNaNf.f);
  return (0);
}
#endif
/* !!!! */


/* Divide 2 complex floating-point numbers by a real number and
 * overwrite input values with the results. */
inline_ void crdiv_2xcf32( xtcomplexfloat * pa0, 
                           xtcomplexfloat * pa1, xtfloat b )
{
  xb_vecN_2xf32 a0,a1,r0,r1,_b,s,t;
  a0 = BBE_MOVN_2XF32_FROMCF32(*pa0);
  a1 = BBE_MOVN_2XF32_FROMCF32(*pa1);
  _b = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(b), 0);
  s = BBE_RECIPN_2XF32(_b);
  t = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(t,s,_b); BBE_MULAN_2XF32(s,t,s);
  r0 = BBE_MULN_2XF32(a0,s); BBE_MULSN_2XF32(a0,r0,_b); BBE_MULAN_2XF32(r0,a0,s);
  r1 = BBE_MULN_2XF32(a1,s); BBE_MULSN_2XF32(a1,r1,_b); BBE_MULAN_2XF32(r1,a1,s);
  *pa0 = BBE_MOVCF32_FROMN_2XF32(r0);
  *pa1 = BBE_MOVCF32_FROMN_2XF32(r1);

} /* crdiv_2xcf32() */

/* Divide 3 complex floating-point numbers by a real number and
 * overwrite input values with the results. */
inline_ void crdiv_3xcf32( xtcomplexfloat * pa0, 
                           xtcomplexfloat * pa1,
                           xtcomplexfloat * pa2, xtfloat b )
{
  xb_vecN_2xf32 a0,a1,a2,r0,r1,r2,_b,s,t;
  a0 = BBE_MOVN_2XF32_FROMCF32(*pa0);
  a1 = BBE_MOVN_2XF32_FROMCF32(*pa1);
  a2 = BBE_MOVN_2XF32_FROMCF32(*pa2);
  _b = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(b), 0);
  s = BBE_RECIPN_2XF32(_b);
  t = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(t,s,_b); BBE_MULAN_2XF32(s,t,s);
  r0 = BBE_MULN_2XF32(a0,s); BBE_MULSN_2XF32(a0,r0,_b); BBE_MULAN_2XF32(r0,a0,s);
  r1 = BBE_MULN_2XF32(a1,s); BBE_MULSN_2XF32(a1,r1,_b); BBE_MULAN_2XF32(r1,a1,s);
  r2 = BBE_MULN_2XF32(a2,s); BBE_MULSN_2XF32(a2,r2,_b); BBE_MULAN_2XF32(r2,a2,s);
  *pa0 = BBE_MOVCF32_FROMN_2XF32(r0);
  *pa1 = BBE_MOVCF32_FROMN_2XF32(r1);
  *pa2 = BBE_MOVCF32_FROMN_2XF32(r2);

} /* crdiv_3xcf32() */

/* Perform a Francis QR step for complex upper-Hessenberg matrix, and
 * optionally update the transform matrix. */
static void francisQRStep ( complex_float * restrict H,
                            complex_float * restrict P,
                            complex_float p, complex_float q, complex_float r,
                            int l, int m, int en, int low, int upp, int N );

/*
 * Francis QR algorithm for complex upper-Hessenberg matrix. 
 * For NxN matrix A, the function computes the full set of N eigenvalues e[N]
 * from its upper-Hessenberg form H. It also outputs a triangular form T and an
 * updated unitary matrix P (optional) such that A == P*T*P'. 
 * Input:
 *   N               Matrix size
 *   low,upp         Lower/upper boundaries of block diagonal structure. 
 *                   Set to 0,N-1 if QR should be applied to full matrix H.
 * Input/Output:
 *   H[N*(N+3)/2-1]  Upper-Hessenberg form (in), triangular form (out). Zero
 *                   elements below the first subdiagonal are not stored.
 *                   Values on the first subdiagonal of triangular output
 *                   form are not specified.
 *   P[N*N]          Transformation matrix reducing original matrix to 
 *                   upper-Hessenberg form (in), or to triangular form T (out).
 *                   This argument is optional, set to zero if transformation
 *                   matrix is not required
 * Output:
 *   e[N]            Eigenvalues of original matrix. Complex eigenvalues
 *                   come in conjugate pairs. If the QR algorithm fails to 
 *                   converge, real and imaginary components of ev[0..N-1]
 *                   are set to NaN.
 * Return Value:
 *   Non-zero if successfull, or zero if the QR algorithm failed to converge.
 * Restrictions:
 *   0<=low<=upp<N
 */

int eigen_hqr_f ( complex_float * restrict e, 
                  complex_float * restrict H, 
                  complex_float * restrict P, 
                  int low, int upp, int N )
{
  /*
   * This implementation is based on the following articles:
   * [1] "The QR Algorithm for Real Hessenberg Matrices" by R.S. Martin, 
   *     G. Petern and J.H. Wilkinson, Handbook for Automatie Computation,
   *     Vol.II Linear Algebra, Contribution II/14.
   * [2] "Eigenvectors of Real and Complex Matrices by LR and QR
   *     triangularizations" by G. Petern and J.H. Wilkinson, Handbook for
   *     Automatie Computation, Vol.II Linear Algebra, Contribution II/15.
   *
   * MATLAB reference code:
   *
   *   function [ev,V,T,S,cnt] = hqr2(H,P,low,upp)
   *   n = length(H); 
   *   if ~exist('P','var'), P = eye(n); end;
   *   if ~exist('low','var'), low = 1; end;
   *   if ~exist('upp','var'), upp = n; end;
   *   % Valudate input matrix, it must be square upper-Hessenberg
   *   msg = 'Argument H must be a square upper-Hessenberg matrix!';
   *   assert(n^2==numel(H),msg);
   *   for i=3:n, assert(all(H(i,1:i-2)==0),msg); end;
   *   assert((1<=low)&&(upp<=n));
   *   EPS = eps(cast(1,'like',H));
   *   % Prepare a vector for eigenvalues.
   *   ev = cast(zeros(n,1),'like',H);
   *   cnt = zeros(n,1);
   *   % Copy the isolated roots 
   *   for i=[1:low-1,upp+1:n], ev(i) = H(i,i); end;
   *   % Loop until all eogenvalues are found, or failed to converge.
   *   its = 0;
   *   t = cast(0,'like',H);
   *   T = H; S = P;
   *   en = upp;
   *   na = en-1;
   *   while en>=low
   *     % Look for a negligible element on the subdiagonal, starting from
   *     % (en,en-1) upwards. If found at l, then the iteration will be performed
   *     % for submatrix T(l:en,l:en), otherwise l=low. Negligibility criterion
   *     % defined by [2] (23).
   *     l = en;
   *     while l>low
   *       if cmag(T(l,l-1))<=EPS*(cmag(T(l-1,l-1))+cmag(T(l,l))), break; end;
   *       l = l-1;
   *     end
   *     % Analyse the current state of the working matrix:
   *     %  A) no eigenvalue revealed. We should either:
   *     %     A1) perform a QR iteration with regular double shift, or
   *     %     A2) perform a QR iteration with exceptional shift, or
   *     %     A3) declare a failure!
   *     %  B) single eigenvalue revealed, or
   *     %  C) a pair of eigenvalues revealed
   *     if l<en-1
   *       % Select a double shift and perform a QR iteration.
   *       if its>0 && mod(its,10)==0
   *         %  Try exceptional shift [2] (19)
   *         t = t+T(en,en);
   *         for i=low:en, T(i,i) = T(i,i)-T(en,en); end;
   *         s = cmag(T(en,na))+cmag(T(na,en-2));
   *         x = 0.75*s; y = 0.75*s; w = -0.4375*s^2;
   *       elseif its<45
   *         % Regular shift
   *         x = T(en,en); y = T(na,na); w = T(na,en)*T(en,na);
   *       else
   *         % Declare a failure
   *         fprintf('Failed to converge at en=%d\n',en);
   *         ev = []; V = []; T = []; S = []; cnt = [];
   *         return;
   *       end
   *       % Look for two small consecutive elements on the subdiagonal. If found,
   *       % the QR iteration should work on a smaller matrix.
   *       for m = en-2:-1:l
   *         z = T(m,m); r = x-z; s = y-z;
   *         p = cdiv(r*s-w,T(m+1,m))+T(m,m+1);
   *         q = T(m+1,m+1)-z-r-s;
   *         r = T(m+2,m+1);
   *         s = cmag(p)+cmag(q)+cmag(r);
   *         p = p/s; q = q/s; r = r/s;
   *         if m>l
   *           a = cmag(T(m,m-1))*(cmag(q)+cmag(r));
   *           b = cmag(p)*(cmag(T(m-1,m-1))+cmag(z)+cmag(T(m+1,m+1)));
   *           if a<=EPS*b, break; end;
   *         end
   *       end
   *       % Perform a double QR step involving rows 1 to n and columns m to n.
   *       %fprintf('%d %d %d %d %g %g %g\n',its+1+sum(cnt(cnt>0)),l,m,en,p,q,r);
   *       [T,S] = FrancisQRStep(T,S,l,m,en,na,n,low,upp,p,q,r);
   *       its = its+1;
   *     elseif l==en
   *       % One eigenvalie revealed
   *       cnt(en) = its;
   *       T(en,en) = t+T(en,en);
   *       ev(en) = T(en,en);
   *       en = na; na = en-1;
   *       its = 0;
   *     else
   *       % Revealed a pair of eigenvalues, which are the two roots of 2x2 
   *       % complex block at bottom right.
   *       assert(l==na);
   *       cnt(en) = its; cnt(na) = -its;
   *       x = T(en,en); y = T(na,na); w = T(na,en)*T(en,na);
   *       p = (y-x)/2; q = p^2+w; z = csqrt(q);
   *       x = x+t; y = y+t; T(en,en) = x; T(na,na) = y;
   *       if real(p*conj(z))<0, z = p-z; else z = p+z; end;
   *       if w~=0
   *         % Quadratic divisor
   *         ev(na) = x+z; ev(en) = x-cdiv(w,z);
   *       else
   *         % Two linear divisors
   *         ev(na) = y; ev(en) = x;
   *       end
   *       % Annihilate T(en,en-1) entry by similarity transformation involving 
   *       % a rotation in the en,en-1 plane. 
   *       x = T(en,na); r = sqrt(x*x'+z*z'); p = x/r; q = z/r;
   *       % Row modification
   *       T(na:en,na:n) = [q',p';-p,q]*T(na:en,na:n);
   *       % Column modification
   *       T(1:en,na:en) = T(1:en,na:en)*[q,-p';p,q'];
   *       % Accumulate transformations: right-multiply S by the rotator.
   *       S(low:upp,na:en) = S(low:upp,na:en)*[q,-p';p,q'];
   *       en = en-2; na = en-1;
   *       its = 0;
   *     end
   *   end
   *   % All eigenvalues found, now perform backsubstitution to determine
   *   % eigenvectors of triangular form T.
   *   V = backsubst(T,n);
   *   % Left-multiply vectors by the transformation matrix to obtain eigenvectors
   *   % of original full matrix.
   *   for j=1:n
   *     m = min(j,upp);
   *     V(low:upp,j) = S(low:upp,low:m)*V(low:m,j);
   *   end
   */
#if 1
  xtcomplexfloat * restrict Ew = (xtcomplexfloat*)e;
  const xtcomplexfloat * restrict H0r;
  const xtcomplexfloat * restrict H1r;
  xtcomplexfloat * restrict H0w;
  xtcomplexfloat * restrict H1w;
  xtcomplexfloat * restrict Hrw = (xtcomplexfloat*)H;
  const xtcomplexfloat * restrict Pr;
  xtcomplexfloat * restrict Pw;

  int i,j,l,m,en,na,its,stp;
  xtcomplexfloat f,g,p,q,r,t,x,y,z,w;
  float32_t a,b,s;
  vbool1 ble,blt;

  const xtcomplexfloat c0f = BBE_CONSTCF32(0);

  /* Extract the isolated roots. */
  for ( i=0    ; i<low; i++ ) e[i] = H[HIDX(i,i)];
  for ( i=upp+1; i<N  ; i++ ) e[i] = H[HIDX(i,i)];

  /* Loop until all eogenvalues are found, or failed to converge. */
  t = c0f; its = 0;
  en = upp; na = en-1;
  while (en>=low) {
    /* Look for a negligible element on the subdiagonal, starting from
     * (en,en-1) upwards. If found at l, then the iteration will be performed
     * for submatrix T(l:en,l:en), otherwise l=low. Negligibility criterion
     * defined by [1] (23). */
    stp = N-low;
    z = *(xtcomplexfloat*)&H[HIDX(low,low)];
    H0r = (xtcomplexfloat*)&H[HIDX(low+1,low+1)];
    for ( l=low, i=low+1; i<=en; i++ ) {
      x = xtcomplexfloat_loadi(H0r, -(int)sz_f32c); /* H[HIDX(i,i-1)] */
      y = z; /* H[HIDX(i-1,i-1)] */
      xtcomplexfloat_loadxp(z, H0r, (stp--)*sz_f32c); /* H[HIDX(i,i)] */
      a = IT_CMAGCF32(x);
      b = XT_MUL_S(EPS, XT_ADD_S(IT_CMAGCF32(y), IT_CMAGCF32(z)));
      XT_MOVNEZ(l, i, BBE_MOVAB1(XT_OLE_S(a,b)));
    }
    /* Analyse the current state of the working matrix. */
    if (l<na) {
      /* No eigenvalue revealed. We should either:
       * A) perform a QR iteration with regular double shift, or
       * B) perform a QR iteration with exceptional shift, or
       * C) declare a failure! */
      if ( its>0 && 0==(its%10) ) {
        /* Try exceptional shift [1] (19) */
        r = Hrw[HIDX(en,en)];
        t = BBE_OPERATOR_ADDCF32(t,r);
        stp = N+1-low;
        H0r = H0w = (xtcomplexfloat*)&H[HIDX(low,low)];
        __Pragma("loop_count min=3");
        for ( i=low; i<=en; i++ ) {
          xtcomplexfloat_loadxp(x, H0r, stp*sz_f32c);
          x = BBE_OPERATOR_SUBCF32(x,r);
          xtcomplexfloat_storexp(x, H0w, (stp--)*sz_f32c);
        }
        s = XT_ADD_S(IT_CMAGCF32(Hrw[HIDX(en,na)]), IT_CMAGCF32(Hrw[HIDX(na,en-2)]));
        x = y = BBE_CMPLXF32(XT_CONST_S(0), 0.75f*s);
        w = BBE_CMPLXF32(XT_CONST_S(0), -0.4375f*s*s);
      } else if ( its<ITS_LIM ) {
        /* Apply the regular shift */
        x = Hrw[HIDX(en,en)];
        y = Hrw[HIDX(na,na)];
        w = BBE_MULCF32(Hrw[HIDX(na,en)], Hrw[HIDX(en,na)]); 
      } else {
        /* Failed to converge! */
        __Pragma("loop_count min=3");
        for ( i=0; i<N; i++ ) {
          Ew[i] = BBE_CMPLXF32(qNaNf.f, qNaNf.f);
        }
        return (0); 
      }
      /* Look for two small consecutive elements on the subdiagonal. If found,
       * the QR iteration should work on a smaller matrix. */
      {
        xtcomplexfloat _p,_q,_r;
        xtcomplexfloat h0,h1,h2,h3,h4,h5,h6,h7;

        m = l; 
        h0 = Hrw[HIDX(m-1,m-1)]; h1 = Hrw[HIDX(m  ,m  )];
        h2 = Hrw[HIDX(m+1,m+1)]; h3 = Hrw[HIDX(m+2,m+2)];
        h4 = Hrw[HIDX(m  ,m-1)]; h5 = Hrw[HIDX(m+1,m  )];
        h6 = Hrw[HIDX(m+2,m+1)]; h7 = Hrw[HIDX(m  ,m+1)];

        /* z = h1; r = csubf(x,z); f = csubf(y,z); */
        z = h1;
        r = BBE_OPERATOR_SUBCF32(x,z);
        f = BBE_OPERATOR_SUBCF32(y,z);
        /* g = csubf(cmulf(r,f),w); */
        g = BBE_NEGCF32(w); BBE_MULACF32(g,r,f);
        /* p = caddf(_cdivf(g,h5), h7); */
        p = BBE_OPERATOR_ADDCF32(IT_DIVCF32(g,h5,0), h7);
        /* q = csubf(csubf(csubf(h2,z),r),f); */
        q = BBE_OPERATOR_SUBCF32(h2,z);
        q = BBE_OPERATOR_SUBCF32(q,r);
        q = BBE_OPERATOR_SUBCF32(q,f);
        r = h6;
        /* s = cmagf(p) + cmagf(q) + cmagf(r); */
        s = XT_ADD_S(XT_ADD_S(IT_CMAGCF32(p), IT_CMAGCF32(q)), IT_CMAGCF32(r));
        /* p = crdivf(p,s); q = crdivf(q,s); r = crdivf(r,s); */
        crdiv_3xcf32(&p,&q,&r,s);

        H0r = (xtcomplexfloat*)&H[HIDX(l+3,l+2)];
        H1r = (xtcomplexfloat*)&H[HIDX(l+1,l+2)]; 
        stp = N-l;
        for ( i=l+1; i<na; i++ ) {
          h0 = h1; h1 = h2; h2 = h3; 
          h3 = xtcomplexfloat_loadi(H0r, sz_f32c);
          h4 = h5; h5 = h6; 
          xtcomplexfloat_loadxp(h6, H0r, (stp-2)*sz_f32c);
          xtcomplexfloat_loadxp(h7, H1r, (stp--)*sz_f32c);

          /* z = h1; r = csubf(x,z); f = csubf(y,z); */
          z = h1;
          _r = BBE_OPERATOR_SUBCF32(x,z);
          f = BBE_OPERATOR_SUBCF32(y,z);
          /* g = csubf(cmulf(r,f),w); */
          g = BBE_NEGCF32(w); BBE_MULACF32(g,_r,f);
          /* p = caddf(_cdivf(g,h5), h7); */
          _p = BBE_OPERATOR_ADDCF32(IT_DIVCF32(g,h5,0), h7);
          /* q = csubf(csubf(csubf(h2,z),r),f); */
          _q = BBE_OPERATOR_SUBCF32(h2,z);
          _q = BBE_OPERATOR_SUBCF32(_q,_r);
          _q = BBE_OPERATOR_SUBCF32(_q,f);
          _r = h6;
          /* s = cmagf(p) + cmagf(q) + cmagf(r); */
          s = XT_ADD_S(XT_ADD_S(IT_CMAGCF32(_p), IT_CMAGCF32(_q)), IT_CMAGCF32(_r));
          /* p = crdivf(p,s); q = crdivf(q,s); r = crdivf(r,s); */
          crdiv_3xcf32(&_p,&_q,&_r,s);

          /* a = cmagf(h4) * (cmagf(q) + cmagf(r)); */
          s = XT_ADD_S(IT_CMAGCF32(_q), IT_CMAGCF32(_r));
          a = XT_MUL_S(IT_CMAGCF32(h4), s);
          /* b = cmagf(p) * (cmagf(h0) + cmagf(z) + cmagf(h2)); */
          s = XT_ADD_S(XT_ADD_S(IT_CMAGCF32(h0), IT_CMAGCF32(z)), IT_CMAGCF32(h2));
          b = XT_MUL_S(IT_CMAGCF32(_p), s);
          /* if (a<=EPS*b) {
           *   m = i; p = p; q = q; r = r;
           * } */
          ble = XT_OLE_S(a, XT_MUL_S(EPS,b));
          XT_MOVNEZ(m,i,BBE_MOVAB1(ble));
          p = BBE_MOVCF32T(_p,p,ble);
          q = BBE_MOVCF32T(_q,q,ble);
          r = BBE_MOVCF32T(_r,r,ble);
        } /* i */
      }
#if 0
      /* !!!! */
      printf( "\nits: %2d  l:%2d  m:%2d  en:%2d  p: (%+.8e,%+.8e)  q: (%+.8e,%+.8e)  r: (%+.8e,%+.8e)",
              its, l, m, en, (float64_t)BBE_CREALCF32(p), (float64_t)BBE_CIMAGCF32(p), 
                             (float64_t)BBE_CREALCF32(q), (float64_t)BBE_CIMAGCF32(q), 
                             (float64_t)BBE_CREALCF32(r), (float64_t)BBE_CIMAGCF32(r) );
      /* !!!! */
#endif
      /* Perform a double QR step involving rows 1 to n and columns m to n. */
      francisQRStep(H, P, MOV_COMPLEX_FLOAT_FROM_CF32(p),
                          MOV_COMPLEX_FLOAT_FROM_CF32(q),
                          MOV_COMPLEX_FLOAT_FROM_CF32(r), 
                          l, m, en, low, upp, N);
      its++;
    } else if (l==na) {
      /* Revealed a pair of eigenvalues. */
      x = Hrw[HIDX(en,en)];
      y = Hrw[HIDX(na,na)];
      w = BBE_MULCF32(Hrw[HIDX(na,en)], Hrw[HIDX(en,na)]);
      f = BBE_OPERATOR_SUBCF32(y,x); p = IT_RCMULCF32(XT_CONST_S(3), f);
      f = w; BBE_MULACF32(f,p,p); z = IT_SQRTCF32(f);
      x = BBE_OPERATOR_ADDCF32(x,t); y = BBE_OPERATOR_ADDCF32(y,t);
      Hrw[HIDX(en,en)] = x; Hrw[HIDX(na,na)] = y;
      /* crealf(p)*crealf(z) + cimagf(p)*cimagf(z) */
      f = BBE_MULMCF32(p,z,0,3); BBE_MULMASCF32(f,p,z,0,12);
      blt = XT_OLT_S(BBE_CREALCF32(f), XT_CONST_S(0));
      BBE_ADDSUBCF32(f,g,p,z);
      z = BBE_MOVCF32T(f,g,blt);
      /* if (cmagf(w)>0) */
      if (BBE_MOVAB1(XT_OLT_S(XT_CONST_S(0), IT_CMAGCF32(w)))) {
        f = IT_DIVCF32(w,z,1);
        Ew[na] = BBE_OPERATOR_ADDCF32(x,z);
        Ew[en] = BBE_OPERATOR_SUBCF32(x,f);
      } else {
        Ew[na] = y; Ew[en] = x;
      }
      if (P) {
        /* Annihilate T(en,en-1) entry by similarity transformation involving 
         * a rotation in the en,en-1 plane. */ 
        x = Hrw[HIDX(en,na)];
        /* s = sqrtf(cabs2f(x)+cabs2f(z)); */
        f = BBE_MULMCF32(x,x,0,3); BBE_MULMASCF32(f,x,x,0,12);
        BBE_MULMASCF32(f,z,z,0,3); BBE_MULMASCF32(f,z,z,0,12); 
        s = XT_SQRT_S(BBE_CREALCF32(f));
        /* p = x/s; q = z/s; */
        p = x; q = z; crdiv_2xcf32(&p,&q,s);
        /* Row modification */
        H0r = H0w = (xtcomplexfloat*)&H[HIDX(na,na)];
        H1r = H1w = (xtcomplexfloat*)&H[HIDX(en,na)];
        __Pragma("loop_count min=2");
        for ( j=na; j<N; j++ ) {
          xtcomplexfloat_loadip(x, H0r, sz_f32c);
          xtcomplexfloat_loadip(y, H1r, sz_f32c);
          z = BBE_MULJCF32(x,q); BBE_MULJACF32(z,y,p);
          w = BBE_MULCF32(y,q); BBE_MULSCF32(w,x,p);
          xtcomplexfloat_storeip(z, H0w, sz_f32c);
          xtcomplexfloat_storeip(w, H1w, sz_f32c);
        }
        __Pragma("no_reorder");
        /* Column modification */
        H0r = H0w = (xtcomplexfloat*)&H[HIDX(0,en)];
        stp = N;
        __Pragma("loop_count min=1");
        for ( i=0; i<=en; i++ ) {
          x = xtcomplexfloat_loadi(H0r, -(int)sz_f32c);
          xtcomplexfloat_loadxp(y, H0r, stp*sz_f32c);
          z = BBE_MULCF32(x,q); BBE_MULACF32(z,y,p);
          w = BBE_MULJCF32(y,q); BBE_MULJSCF32(w,x,p);
          xtcomplexfloat_storei (z, H0w, -(int)sz_f32c);
          xtcomplexfloat_storexp(w, H0w, (stp--)*sz_f32c);
        }
        /* Accumulate transformations: right-multiply P by the rotator. */
        Pr = Pw = (xtcomplexfloat*)&P[low*N+en];
        __Pragma("loop_count min=2");
        for ( i=low; i<=upp; i++ ) {
          x = xtcomplexfloat_loadi(Pr, -(int)sz_f32c);
          xtcomplexfloat_loadxp(y, Pr, N*sz_f32c);
          z = BBE_MULCF32(x,q); BBE_MULACF32(z,y,p);
          w = BBE_MULJCF32(y,q); BBE_MULJSCF32(w,x,p);
          xtcomplexfloat_storei (z, Pw, -(int)sz_f32c);
          xtcomplexfloat_storexp(w, Pw, N*sz_f32c);
        }
      } /* P */
      en = na-1; na = en-1; its = 0;
    } else { /* l */
      /* Revealed a standalone eigenvalue. */
      Ew[en] = Hrw[HIDX(en,en)] = BBE_OPERATOR_ADDCF32(Hrw[HIDX(en,en)], t);
      en = na; na = en-1; its = 0;
    }
  } /* en */
#else
  int i,j,l,m,en,na,its;
  complex_float f,g,p,q,r,t,x,y,z,w;
  float32_t a,b,s;
  const complex_float c0f = _makecomplexf(0.f,0.f);

  /* Extract the isolated roots. */
  for ( i=0    ; i<low; i++ ) e[i] = H[HIDX(i,i)];
  for ( i=upp+1; i<N  ; i++ ) e[i] = H[HIDX(i,i)];

  /* Loop until all eogenvalues are found, or failed to converge. */
  t = c0f; its = 0;
  en = upp; na = en-1;
  while (en>=low) {
    /* Look for a negligible element on the subdiagonal, starting from
     * (en,en-1) upwards. If found at l, then the iteration will be performed
     * for submatrix T(l:en,l:en), otherwise l=low. Negligibility criterion
     * defined by [1] (23). */
    for ( l=en; l>low; l-- ) {
      x = H[HIDX(l,l-1)];
      y = H[HIDX(l-1,l-1)];
      z = H[HIDX(l,l)];
      if ( cmagf(x) <= EPS*(cmagf(y)+cmagf(z)) ) break;
    }
    /* Analyse the current state of the working matrix. */
    if (l<na) {
      /* No eigenvalue revealed. We should either:
       * A) perform a QR iteration with regular double shift, or
       * B) perform a QR iteration with exceptional shift, or
       * C) declare a failure! */
      if ( its>0 && 0==(its%10) ) {
        /* Try exceptional shift [1] (19) */
        r = H[HIDX(en,en)]; t = caddf(t,r);
        for ( i=low; i<=en; i++ ) H[HIDX(i,i)] = csubf(H[HIDX(i,i)], r);
        s = cmagf(H[HIDX(en,na)]) + cmagf(H[HIDX(na,en-2)]);
        x = y = _makecomplexf(0.75f*s, 0.f); 
        w = _makecomplexf(-0.4375f*s*s, 0.f);
      } else if ( its<ITS_LIM ) {
        /* Apply the regular shift */
        x = H[HIDX(en,en)];
        y = H[HIDX(na,na)];
        w = cmulf(H[HIDX(na,en)], H[HIDX(en,na)]); 
      } else {
        /* Failed to converge! */
        for ( i=0; i<N; i++ ) {
          e[i] = _makecomplexf( qNaNf.f, qNaNf.f );
        }
        return (0); 
      }
      /* Look for two small consecutive elements on the subdiagonal. If found,
       * the QR iteration should work on a smaller matrix. */
      m = na;
      do {
        m--; z = H[HIDX(m,m)]; r = csubf(x,z); f = csubf(y,z);
        g = csubf(cmulf(r,f),w);
        p = caddf(_cdivf(g,H[HIDX(m+1,m)]), H[HIDX(m,m+1)]);
        q = csubf(csubf(csubf(H[HIDX(m+1,m+1)],z),r),f);
        r = H[HIDX(m+2,m+1)];
        s = cmagf(p) + cmagf(q) + cmagf(r);
        p = crdivf(p,s);  q = crdivf(q,s); r = crdivf(r,s);
        if (m>l) {
          a = cmagf(H[HIDX(m,m-1)]) * (cmagf(q) + cmagf(r));
          b = cmagf(p) * (cmagf(H[HIDX(m-1,m-1)]) + cmagf(z) + 
                          cmagf(H[HIDX(m+1,m+1)]));
          if (a<=EPS*b) break;
        }
      } while (m>l);
      /* !!!! */
#if 0
      printf( "\nits: %2d  l:%2d  m:%2d  en:%2d  p: (%+.8e,%+.8e)  q: (%+.8e,%+.8e)  r: (%+.8e,%+.8e)",
              its, l, m, en, crealf(p), cimagf(p), crealf(q), cimagf(q), crealf(r), cimagf(r) );
#endif
      /* !!!! */
      /* Perform a double QR step involving rows 1 to n and columns m to n. */
      francisQRStep(H,P,p,q,r,l,m,en,low,upp,N); its++;
    } else if (l==na) {
      /* Revealed a pair of eigenvalues. */
      x = H[HIDX(en,en)];
      y = H[HIDX(na,na)];
      w = cmulf(H[HIDX(na,en)], H[HIDX(en,na)]);
      p = rcmulf(.5f, csubf(y,x));
      z = _csqrtf(caddf(cmulf(p,p),w));
      x = caddf(x,t); y = caddf(y,t);
      H[HIDX(en,en)] = x; H[HIDX(na,na)] = y;
      if ( crealf(p)*crealf(z) + cimagf(p)*cimagf(z) < 0 ) {
        z = csubf(p,z);
      } else {
        z = caddf(p,z);
      }
      if (cmagf(w)>0) {
        /* Quadratic divisor */
        e[na] = caddf(x,z);
        e[en] = csubf(x, _cdivf(w,z));
      } else {
        /* Two linear divisors */
        e[na] = y; e[en] = x;
      }
      if (P) {
        /* Annihilate T(en,en-1) entry by similarity transformation involving 
         * a rotation in the en,en-1 plane. */ 
        x = H[HIDX(en,na)]; s = sqrtf(cabs2f(x)+cabs2f(z));
        p = crdivf(x,s); q = crdivf(z,s);
        /* Row modification */
        for ( j=na; j<N; j++ ) {
          z = H[HIDX(na,j)];
          H[HIDX(na,j)] = caddf(cmuljf(z,q), cmuljf(H[HIDX(en,j)],p));
          H[HIDX(en,j)] = csubf(cmulf(q,H[HIDX(en,j)]), cmulf(p,z));
        }
        /* Column modification */
        for ( i=0; i<=en; i++ ) {
          z = H[HIDX(i,na)];
          H[HIDX(i,na)] = caddf(cmulf(q,z), cmulf(p,H[HIDX(i,en)]));
          H[HIDX(i,en)] = csubf(cmuljf(H[HIDX(i,en)],q), cmuljf(z,p));
        }
        /* Accumulate transformations: right-multiply S by the rotator. */
        for ( i=low; i<=upp; i++ ) {
          z = P[i*N+na]; 
          P[i*N+na] = caddf(cmulf(q,z), cmulf(p,P[i*N+en]));
          P[i*N+en] = csubf(cmuljf(P[i*N+en],q), cmuljf(z,p));
        }
      }
      en = na-1; na = en-1; its = 0;
    } else {
      /* Revealed a standalone eigenvalue. */
      H[HIDX(en,en)] = caddf(H[HIDX(en,en)], t);
      e[en] = H[HIDX(en,en)];
      en = na; na = en-1; its = 0;
    }
  } /* en */
#endif
  return (1);

} /* eigen_hqr_f() */

/* Perform a Francis QR step for complex upper-Hessenberg matrix, and
 * optionally update the transform matrix. */
void francisQRStep ( complex_float * restrict H,
                     complex_float * restrict P,
                     complex_float p, complex_float q, complex_float r,
                     int l, int m, int en, int low, int upp, int N )
{
  /*
   * MATLAB reference code:
   *
   *   function [H,P] = FrancisQRStep(H,P,l,m,en,na,n,low,upp,p,q,r)
   *   assert((l<=m)&&(m<=en-2)&&(low+2<=en)&&(en<=upp));
   *   t = 1; g = 0;
   *   for k = m:na
   *     notLast = k<na;
   *     if t==0, break; end;
   *     p = p/t; q = q/t; r = r/t;
   *     s = sqrt(p'*p+q'*q+r'*r);
   *     h = cnorm(p); % Note that cnorm(0) == 1
   *     if k>m
   *       H(k,k-1) = -h*(s*t);
   *     elseif l<m
   *       H(k,k-1) = -H(k,k-1);
   *     end
   *     t = cabs(p); x = t/s+1; y = h'*q/s; z = h'*r/s;
   *     q = h*q'/(t+s); r = h*r'/(t+s);
   *     % Accumulate transformations: right-multiply P by the reflector.
   *     for i=low:upp
   *       p = x*P(i,k)+y*P(i,k+1);
   *       if notLast
   *         p = p+z*P(i,k+2);
   *         P(i,k+2) = P(i,k+2)-p*r;
   *       end
   *       P(i,k+1) = P(i,k+1)-p*q;
   *       P(i,k) = P(i,k)-p;
   *     end
   *     % Row modification
   *     % k-th column
   *     p = H(k,k)+q*H(k+1,k)+r*g; 
   *     g = g-p*z;
   *     H(k+1,k) = H(k+1,k)-p*y;
   *     H(k,k) = H(k,k)-p*x;
   *     % Columns k+1..n
   *     for j=k+1:n
   *       p = H(k,j)+q*H(k+1,j);
   *       if notLast
   *         p = p+r*H(k+2,j);
   *         H(k+2,j) = H(k+2,j)-p*z;
   *       end
   *       H(k+1,j) = H(k+1,j)-p*y;
   *       H(k,j) = H(k,j)-p*x;
   *     end
   *     % Column modification
   *     % Rows l..k+1
   *     for i=1:k+1
   *       p = x*H(i,k)+y*H(i,k+1);
   *       if notLast
   *         p = p+z*H(i,k+2);
   *         H(i,k+2) = H(i,k+2)-p*r;
   *       end
   *       H(i,k+1) = H(i,k+1)-p*q;
   *       H(i,k) = H(i,k)-p;
   *     end
   *     % Row k+2
   *     if notLast
   *       p = x*g+y*H(k+2,k+1)+z*H(k+2,k+2);
   *       H(k+2,k+2) = H(k+2,k+2)-p*r;
   *       H(k+2,k+1) = H(k+2,k+1)-p*q;
   *       f = g-p;
   *     end
   *     % Row k+3
   *     if k<na-1
   *       p = z*H(k+3,k+2);
   *       H(k+3,k+2) = H(k+3,k+2)-p*r;
   *       g = -p*q;
   *       r = -p;
   *     else
   *       r = 0;
   *     end
   *     p = H(k+1,k);
   *     q = f; t = cmag(p)+cmag(q)+cmag(r);
   *   end
   */
#if 1
        xtcomplexfloat * restrict Hrw = (xtcomplexfloat*)H;
  const xtcomplexfloat * restrict H0r;
  const xtcomplexfloat * restrict H1r;
  const xtcomplexfloat * restrict H2r;
        xtcomplexfloat * restrict H0w;
        xtcomplexfloat * restrict H1w;
        xtcomplexfloat * restrict H2w;
  const xtcomplexfloat * restrict Pr;
        xtcomplexfloat * restrict Pw;

  int i,j,k;
  int na = en-1;
  int stp;
  xtcomplexfloat _p,_q,_r;
  xtcomplexfloat a,b,c,f,g,h,y,z;
  float32_t s,t,x,u;
  vbool1 ble;
  const xtcomplexfloat c1f = BBE_CMPLXF32(BBE_CONSTN_2XF32(0),BBE_CONSTN_2XF32(1));
  const xtcomplexfloat c0f = BBE_CONSTCF32(0);
  
  NASSERT( (l<=m) && (m<=en-2) && ((low+2)<=en) && (en<=upp) && (upp<N) );

  _p = MOV_CF32_FROM_COMPLEX_FLOAT(p);
  _q = MOV_CF32_FROM_COMPLEX_FLOAT(q);
  _r = MOV_CF32_FROM_COMPLEX_FLOAT(r);

  t = XT_CONST_S(1);
  g = c0f;
  for ( k=m; k<en && t>=FLT_MIN; k++ ) { 
    /* Construct a Householder reflector from the vector (p,q,r). */
    if (k>m) {
      crdiv_3xcf32(&_p,&_q,&_r,t);
    }
    /* s = sqrtf(cabs2f(p) + cabs2f(q) + cabs2f(r)); */
    s = XT_SQRT_S(XT_ADD_S(XT_ADD_S(IT_ABS2CF32(_p), IT_ABS2CF32(_q)), IT_ABS2CF32(_r)));
    u = IT_ABSCF32(_p); ble = XT_OLE_S(FLT_MIN, u);
    h = IT_CRDIVCF32(_p,u); h = BBE_MOVCF32T(h, c1f, ble);
    if (k>m) {
      Hrw[HIDX(k,k-1)] = IT_RCMULCF32(-s*t,h);
    } else if (l<m) {
      Hrw[HIDX(k,k-1)] = BBE_NEGCF32(Hrw[HIDX(k,k-1)]);
    }
    x = XT_ADD_S(IT_FDIVF32(u,s,0), XT_CONST_S(1));
    y = BBE_MULJCF32(_q,h); z = BBE_MULJCF32(_r,h);
    _q = BBE_CONJCF32(y); _r = BBE_CONJCF32(z);
    crdiv_2xcf32(&_q,&_r,s+u);
    crdiv_2xcf32(&y,&z,s);
    /*----------------------------------------------------------------*
     * Accumulate transformations: right-multiply P by the reflector. *
     *----------------------------------------------------------------*/
    if (P) {
      if (k<na) {
        Pr = Pw= (xtcomplexfloat*)&P[low*N+k+2];
        __Pragma("loop_count min=1");
        for ( i=low; i<=upp; i++ ) {
          a = xtcomplexfloat_loadi(Pr, -2*(int)sz_f32c); /* P[i*N+k] */
          b = xtcomplexfloat_loadi(Pr, -1*(int)sz_f32c); /* P[i*N+k+1] */
          xtcomplexfloat_loadxp(c, Pr, N*sz_f32c);       /* P[i*N+k+2] */
          _p = IT_RCMULCF32(x,a);
          BBE_MULACF32(_p,y,b);
          BBE_MULACF32(_p,z,c);
          a = BBE_OPERATOR_SUBCF32(a,_p);
          BBE_MULSCF32(b,_p,_q);
          BBE_MULSCF32(c,_p,_r);
          xtcomplexfloat_storei (a, Pw, -2*(int)sz_f32c);
          xtcomplexfloat_storei (b, Pw, -1*(int)sz_f32c);
          xtcomplexfloat_storexp(c, Pw, N*sz_f32c);
        }
      } else {  /* k */
        Pr = Pw= (xtcomplexfloat*)&P[low*N+k+1];
        __Pragma("loop_count min=1");
        for ( i=low; i<=upp; i++ ) {
          a = xtcomplexfloat_loadi(Pr, -(int)sz_f32c); /* P[i*N+k] */
          xtcomplexfloat_loadxp(b, Pr, N*sz_f32c);       /* P[i*N+k+1] */
          _p = IT_RCMULCF32(x,a);
          BBE_MULACF32(_p,y,b);
          a = BBE_OPERATOR_SUBCF32(a,_p);
          BBE_MULSCF32(b,_p,_q);
          xtcomplexfloat_storei (a, Pw, -(int)sz_f32c);
          xtcomplexfloat_storexp(b, Pw, N*sz_f32c);
        }
      } /* k */
    }
    /*----------------------------------------------------------------*
     * Row modification, involves columns k through N-1.              *
     *----------------------------------------------------------------*/
    /*========================================*
     * k-th column                            */
    /* p = H(k,k)+q*H(k+1,k)+r*g; */
    _p = Hrw[HIDX(k,k)]; BBE_MULACF32(_p,_q,Hrw[HIDX(k+1,k)]); BBE_MULACF32(_p,_r,g);
    /* g = g-p*z; */
    BBE_MULSCF32(g,_p,z);
    /* H(k+1,k) = H(k+1,k)-p*y; */
    a = Hrw[HIDX(k+1,k)]; BBE_MULSCF32(a,_p,y); Hrw[HIDX(k+1,k)] = a;
    /* H(k,k) = H(k,k)-p*x; */
    a = Hrw[HIDX(k,k)]; IT_RCMULSCF32(a,x,_p); Hrw[HIDX(k,k)] = a;
    /*========================================*
     * Columns k+1..N-1 or k+1..en            */
    i = ( P ? N-1 : en );
    if (k<na) {
      H0r = H0w = (xtcomplexfloat*)&H[HIDX(k  ,k+1)];
      H1r = H1w = (xtcomplexfloat*)&H[HIDX(k+1,k+1)];
      H2r = H2w = (xtcomplexfloat*)&H[HIDX(k+2,k+1)];
      __Pragma("loop_count min=1");
      for ( j=k+1; j<=i; j++ ) {
        xtcomplexfloat_loadip(a, H0r, sz_f32c);
        xtcomplexfloat_loadip(b, H1r, sz_f32c);
        xtcomplexfloat_loadip(c, H2r, sz_f32c);
        _p = a; BBE_MULACF32(_p,_q,b);
        BBE_MULACF32(_p,_r,c);
        IT_RCMULSCF32(a,x,_p);
        BBE_MULSCF32(b,y,_p);
        BBE_MULSCF32(c,z,_p);
        xtcomplexfloat_storeip(a, H0w, sz_f32c);
        xtcomplexfloat_storeip(b, H1w, sz_f32c);
        xtcomplexfloat_storeip(c, H2w, sz_f32c);
      }
    } else { /* k */
      H0r = H0w = (xtcomplexfloat*)&H[HIDX(k  ,k+1)];
      H1r = H1w = (xtcomplexfloat*)&H[HIDX(k+1,k+1)];
      __Pragma("loop_count min=1");
      for ( j=k+1; j<=i; j++ ) {
        xtcomplexfloat_loadip(a, H0r, sz_f32c);
        xtcomplexfloat_loadip(b, H1r, sz_f32c);
        _p = a; BBE_MULACF32(_p,_q,b);
        IT_RCMULSCF32(a,x,_p);
        BBE_MULSCF32(b,y,_p);
        xtcomplexfloat_storeip(a, H0w, sz_f32c);
        xtcomplexfloat_storeip(b, H1w, sz_f32c);
      }
    } /* k */
    __Pragma("no_reorder");
    /*----------------------------------------------------------------*
     * Column modification. Involves rows 0..k+3 (update of P is      *
     * requested) or l-4..k+3 (no need for transform accumulation)    *
     *----------------------------------------------------------------*/
    /*========================================*
     * Rows 0 (or l-4) through k+1            */
    j = ( P || (l<4) ? 0 : l-4 );
    if (k<na) {
      NASSERT(j<=k+1);
      H0r = H0w = (xtcomplexfloat*)&H[HIDX(j,k+2)];
      stp = N-j;
      __Pragma("loop_count min=1");
      for ( i=j; i<=k+1; i++ ) {
        a = xtcomplexfloat_loadi(H0r, -2*(int)sz_f32c);
        b = xtcomplexfloat_loadi(H0r, -1*(int)sz_f32c);
        xtcomplexfloat_loadxp(c, H0r, stp*sz_f32c);
        _p = IT_RCMULCF32(x,a);
        BBE_MULACF32(_p,y,b);
        BBE_MULACF32(_p,z,c);
        a = BBE_OPERATOR_SUBCF32(a,_p);
        BBE_MULSCF32(b,_p,_q);
        BBE_MULSCF32(c,_p,_r);
        xtcomplexfloat_storei (a, H0w, -2*(int)sz_f32c);
        xtcomplexfloat_storei (b, H0w, -1*(int)sz_f32c);
        xtcomplexfloat_storexp(c, H0w, (stp--)*sz_f32c);
      }
    } else { /* k */
      NASSERT(j<=k+1);
      H0r = H0w = (xtcomplexfloat*)&H[HIDX(j,k+1)];
      stp = N-j;
      __Pragma("loop_count min=1");
      for ( i=j; i<=k+1; i++ ) {
        a = xtcomplexfloat_loadi(H0r, -1*(int)sz_f32c);
        xtcomplexfloat_loadxp(b, H0r, stp*sz_f32c);
        _p = IT_RCMULCF32(x,a);
        BBE_MULACF32(_p,y,b);
        a = BBE_OPERATOR_SUBCF32(a,_p);
        BBE_MULSCF32(b,_p,_q);
        xtcomplexfloat_storei (a, H0w, -1*(int)sz_f32c);
        xtcomplexfloat_storexp(b, H0w, (stp--)*sz_f32c);
      }
    } /* k */
    /*========================================*
     * Row k+2                                */
    if (k<na) {
      /* p = x*g+y*H(k+2,k+1)+z*H(k+2,k+2); */
      _p = IT_RCMULCF32(x,g); 
      BBE_MULACF32(_p, y, Hrw[HIDX(k+2,k+1)]); 
      BBE_MULACF32(_p, z, Hrw[HIDX(k+2,k+2)]); 
      /* H(k+2,k+2) = H(k+2,k+2)-p*r; */
      a = Hrw[HIDX(k+2,k+2)]; BBE_MULSCF32(a,_p,_r); Hrw[HIDX(k+2,k+2)] = a;
      /* H(k+2,k+1) = H(k+2,k+1)-p*q; */
      a = Hrw[HIDX(k+2,k+1)]; BBE_MULSCF32(a,_p,_q); Hrw[HIDX(k+2,k+1)] = a;
      /* f = g-p; */
      f = BBE_OPERATOR_SUBCF32(g,_p);
    }
    /*========================================*
     * Row k+3                                */
    if (k<na-1) {
      /* p = z*H(k+3,k+2); */
      _p = BBE_MULCF32(z, Hrw[HIDX(k+3,k+2)]);
      /* H(k+3,k+2) = H(k+3,k+2)-p*r; */
      a = Hrw[HIDX(k+3,k+2)]; BBE_MULSCF32(a,_p,_r); Hrw[HIDX(k+3,k+2)] = a;
      /* g = -p*q; */
      g = BBE_NEGCF32(BBE_MULCF32(_p,_q));
      /* r = -p; */
      _r = BBE_NEGCF32(_p);
    } else {
      _r = c0f;
    }
    __Pragma("no_reorder");
    _p = Hrw[HIDX(k+1,k)];
    _q = f; t = XT_ADD_S(XT_ADD_S(IT_CMAGCF32(_p), IT_CMAGCF32(_q)), IT_CMAGCF32(_r));
  } /* k */
#else
  int i,j,k;
  int na = en-1;
  complex_float f,g,h,y,z;
  float32_t s,t,x;
  const complex_float c0f = _makecomplexf(0.f,0.f);
  
  NASSERT( (l<=m) && (m<=en-2) && ((low+2)<=en) && (en<=upp) && (upp<N) );

  t = 1.f; g = c0f;
  for ( k=m; k<en && t!=0.f; k++ ) {
    /* Construct a Householder reflector from the vector (p,q,r). */
    if (k>m) { p = crdivf(p,t); q = crdivf(q,t); r = crdivf(r,t); }
    s = sqrtf(cabs2f(p) + cabs2f(q) + cabs2f(r));
    h = cnormf(p); /* Note that cnormf(0) == 1 */
    if (k>m) {
      H[HIDX(k,k-1)] = rcmulf(-s*t,h);
    } else if (l<m) {
      H[HIDX(k,k-1)] = cnegf(H[HIDX(k,k-1)]);
    }
    t = _cabsf(p); x = t/s+1.f;
    y = cmuljf(q,h); z = cmuljf(r,h);
    q = crdivf(_conjf(y),s+t); r = crdivf(_conjf(z),s+t);
    y = crdivf(y,s); z = crdivf(z,s);
    /*----------------------------------------------------------------*
     * Accumulate transformations: right-multiply P by the reflector. *
     *----------------------------------------------------------------*/
    if (P) {
      for ( i=low; i<=upp; i++ ) {
        p = caddf(rcmulf(x, P[i*N+k]), cmulf(y, P[i*N+k+1]));
        if (k<na) {
          p = caddf(p, cmulf(z, P[i*N+k+2]));
          P[i*N+k+2] = csubf(P[i*N+k+2],  cmulf(p,r));
        }
        P[i*N+k+1] = csubf(P[i*N+k+1], cmulf(p,q));
        P[i*N+k] = csubf(P[i*N+k], p);
      }
    }
    /*----------------------------------------------------------------*
     * Row modification, involves columns k through N-1.              *
     *----------------------------------------------------------------*/
    /* k-th column */
    p = caddf(caddf(H[HIDX(k,k)], cmulf(q, H[HIDX(k+1,k)])), cmulf(r,g));
    g = csubf(g, cmulf(p,z));
    H[HIDX(k+1,k)] = csubf(H[HIDX(k+1,k)], cmulf(p,y));
    H[HIDX(k,k)] = csubf(H[HIDX(k,k)], rcmulf(x,p));
    /* Columns k+1..N-1 or k+1..en */
    i = ( P ? N-1 : en );
    for ( j=k+1; j<=i; j++ ) {
      p = caddf(H[HIDX(k,j)], cmulf(q, H[HIDX(k+1,j)]));
      if (k<na) {
        p = caddf(p, cmulf(r, H[HIDX(k+2,j)]));
        H[HIDX(k+2,j)] = csubf(H[HIDX(k+2,j)], cmulf(p,z));
      }
      H[HIDX(k+1,j)] = csubf(H[HIDX(k+1,j)], cmulf(p,y));
      H[HIDX(k,j)] = csubf(H[HIDX(k,j)], rcmulf(x,p));
    }
    /*----------------------------------------------------------------*
     * Column modification. Involves rows 0..k+3 (update of P is      *
     * requested) or l-3..k+3 (no need for transform accumulation)    *
     *----------------------------------------------------------------*/
    /* Rows 0 (or l-3) through k+1 */
    j = ( P || (l<3) ? 0 : l-3 );
    for ( i=j; i<=k+1; i++ ) {
      p = caddf(rcmulf(x, H[HIDX(i,k)]), cmulf(y, H[HIDX(i,k+1)]));
      if (k<na) {
        p = caddf(p, cmulf(z, H[HIDX(i,k+2)]));
        H[HIDX(i,k+2)] = csubf(H[HIDX(i,k+2)], cmulf(p,r));
      }
      H[HIDX(i,k+1)] = csubf(H[HIDX(i,k+1)], cmulf(p,q));
      H[HIDX(i,k)] = csubf(H[HIDX(i,k)], p);
    }
    /* Row k+2 */
    if (k<na) {
      p = caddf(caddf(rcmulf(x,g), cmulf(y, H[HIDX(k+2,k+1)])), cmulf(z, H[HIDX(k+2,k+2)]));
      H[HIDX(k+2,k+2)] = csubf(H[HIDX(k+2,k+2)], cmulf(p,r));
      H[HIDX(k+2,k+1)] = csubf(H[HIDX(k+2,k+1)], cmulf(p,q));
      f = csubf(g,p);
    }
    /* Row k+3 */
    if (k<na-1) {
      p = cmulf(z, H[HIDX(k+3,k+2)]);
      H[HIDX(k+3,k+2)] = csubf(H[HIDX(k+3,k+2)], cmulf(p,r));
      g = cnegf(cmulf(p,q)); r = cnegf(p);
    } else {
      r = c0f;
    }
    p = H[HIDX(k+1,k)];
    q = f; t = cmagf(p) + cmagf(q) + cmagf(r);
  } /* k */
#endif
} /* francisQRStep() */

#endif /* HAVE_VFPU */

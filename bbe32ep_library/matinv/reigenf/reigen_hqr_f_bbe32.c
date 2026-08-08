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
    Francis QR algorithm for real upper-Hessenberg matrix
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* #include <stdio.h> !!!! */

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

#if HAVE_VFPU

#define sz_f32   sizeof(float32_t)

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
#endif

/* Perform a Francis QR step for real upper-Hessenberg matrix, and
 * optionally update the transform matrix. */
static void francisQRStep ( float32_t * restrict H,
                            float32_t * restrict P,
                            float32_t p, float32_t q, float32_t r,
                            int l, int m, int en, int low, int upp, int N );

/*
 * Francis QR algorithm for real upper-Hessenberg matrix. 
 * For NxN matrix A, the function computes the full set of N eigenvalues e[N]
 * from its upper-Hessenberg form H. It also outputs a quasi-triangular form T 
 * and an updated orthogonal matrix P (optional) such that A == P*T*P'. 
 * Input:
 *   N               Matrix size
 *   low,upp         Lower/upper boundaries of block diagonal structure. 
 *                   Set to 0,N-1 if QR should be applied to full matrix H.
 * Input/Output:
 *   H[N*(N+3)/2-1]  Upper-Hessenberg form (in), quasi-triangular form with
 *                   2x2 blocks on the main diagonal (out). Zero elements
 *                   below the first subdiagonal are not stored.
 *   P[N*N]          Transformation matrix reducing original matrix to 
 *                   upper-Hessenberg form (in), or to quasi-triangular
 *                   form T (out). This argument is optional, set to zero
 *                   if transformation matrix is not required
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

int reigen_hqr_f( complex_float * restrict e, 
                  float32_t     * restrict H, 
                  float32_t     * restrict P, 
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
   *   function [ev,V,T,S,cnt] = rhqr2(H,P,low,upp)
   *   n = length(H);
   *   if ~exist('P','var'), P = eye(n); end;
   *   if ~exist('low','var'), low = 1; end;
   *   if ~exist('upp','var'), upp = n; end;
   *   % Validate the input matrix: it must be real, square, and upper-Hessenberg.
   *   msg = 'Argument H must be a real square upper-Hessenberg matrix!';
   *   assert(isreal(H),msg);
   *   assert(n^2==numel(H),msg);
   *   for i=3:n, assert(all(H(i,1:i-2)==0),msg); end;
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
   *       if abs(T(l,l-1))<=EPS*(abs(T(l-1,l-1))+abs(T(l,l))), break; end;
   *       l = l-1;
   *     end
   *     % Analyse the current state of the working matrix:
   *     %  A) no eigenvalue revealed. We should either:
   *     %     A1) perform a QR iteration with regular double shift, or
   *     %     A2) perform a QR iteration with exceptional shift, or
   *     %     A3) declare a failure!
   *     %  B) single eigenvalue revealed, or
   *     %  C) a pair of eigenvalues revealed
   *     if l<na
   *       % Select a double shift and perform a QR iteration.
   *       if its>0 && mod(its,10)==0
   *         % Try exceptional shift [2] (19)
   *         t = t+T(en,en);
   *         for i=low:en, T(i,i) = T(i,i)-T(en,en); end;
   *         s = abs(T(en,na))+abs(T(na,en-2));
   *         x = 0.75*s; y = 0.75*s; w = -0.4375*s^2;
   *       elseif its<45
   *         % Regular shift
   *         x = T(en,en); y = T(na,na); w = T(na,en)*T(en,na);
   *       else
   *         % Declare a failure
   *         fprintf('Failed to converge at n=%d\n',n);
   *         ev = []; V = []; T = []; cnt = [];
   *         return;
   *       end
   *       % Look for two small consecutive elements on the subdiagonal. If found,
   *       % the QR iteration should work on a smaller matrix.
   *       for m = en-2:-1:l
   *         z = T(m,m); r = x-z; s = y-z;
   *         p = (r*s-w)/T(m+1,m)+T(m,m+1);
   *         q = T(m+1,m+1)-z-r-s;
   *         r = T(m+2,m+1);
   *         s = abs(p)+abs(q)+abs(r);
   *         p = p/s; q = q/s; r = r/s;
   *         if m>l
   *           a = abs(T(m,m-1))*(abs(q)+abs(r));
   *           b = abs(p)*(abs(T(m-1,m-1))+abs(z)+abs(T(m+1,m+1)));
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
   *       T(en,en) = T(en,en)+t;
   *       ev(en) = T(en,en);
   *       en = na;
   *       na = en-1;
   *       its = 0;
   *     else
   *       % Revealed a pair of eigenvalues
   *       assert(l==na);
   *       cnt(en) = its; cnt(na) = -its;
   *       x = T(en,en); y = T(na,na); w = T(na,en)*T(en,na);
   *       p = (y-x)/2; q = p^2+w; z = sqrt(abs(q));
   *       x = x+t; T(en,en) = x; y = y+t; T(na,na) = y; 
   *       if q>=0
   *         % Pair of real eigenvalues
   *         if p<0, z = p-z; else z = p+z; end;
   *         if w~=0
   *           % Quadratic divisor
   *           ev(na) = x+z; ev(en) = x-w/z;
   *         else
   *           % Two linear divisors
   *           ev(na) = y; ev(en) = x;
   *         end
   *         % Annihilate T(en,en-1) entry by similarity transformation involving 
   *         % a rotation in the en,en-1 plane. 
   *         x = T(en,na); r = sqrt(x^2+z^2); p = x/r; q = z/r;
   *         % Row modification
   *         for j=na:n
   *           z = T(na,j); T(na,j) = q*z+p*T(en,j);
   *           T(en,j) = q*T(en,j)-p*z;
   *         end
   *         % Column modification
   *         for i=1:en
   *           z = T(i,na); T(i,na) = q*z+p*T(i,en);
   *           T(i,en) = q*T(i,en)-p*z;
   *         end
   *         % Accumulate transformations: right-multiply S by the rotator.
   *         for i=low:upp
   *           z = S(i,na); S(i,na) = q*z+p*S(i,en);
   *           S(i,en) = q*S(i,en)-p*z;
   *         end
   *       else
   *         % Complex eigenvalue and its conjugate
   *         ev(na) = x+p+1j*z; ev(en) = x+p-1j*z;
   *       end
   *       en = en-2;
   *       na = en-1;
   *       its = 0;
   *     end
   *   end
   *   % All eigenvalues found, now perform backsubstitution to determine
   *   % eigenvectors of quasitriangular form T. This function replicates
   *   % the procedure described in [3].
   *   V = backsubst(T,ev,n);
   *   % Left-multiply vectors by the transformation matrix to obtain eigenvectors
   *   % of original full matrix.
   *   for j=1:n
   *     if imag(ev(j))>0
   *       m = min(j+1,upp);
   *       % Transform re/im parts of j-th vector
   *       V(low:upp,[j,j+1]) = S(low:upp,low:m)*V(low:m,[j,j+1]);
   *       % Combine re/im parts into complex data.
   *       V(:,[j,j+1]) = V(:,[j,j+1])*[1,1;1j,-1j];
   *     elseif imag(ev(j))==0
   *       m = min(j,upp);
   *       V(low:upp,j) = S(low:upp,low:m)*V(low:m,j);
   *     end
   *   end
   */

#if 1
  const xtfloat * restrict H0r;
  const xtfloat * restrict H1r;
  const xtfloat * restrict Pr;
  xtfloat * restrict Hw;
  xtfloat * restrict Pw;
  xtcomplexfloat * restrict Ew = (xtcomplexfloat*)e;

  int i,j,l,m,en,na,its,stp;
  float32_t a,b,p,q,r,s,t;
  float32_t x,y,z,w;
  vbool1 ble,beqz,bltz;

  NASSERT(N>1);
  /* Extract the isolated roots. */
  for ( i=0    ; i<low; i++ ) Ew[i] = BBE_CMPLXF32(H[HIDX(i,i)], XT_CONST_S(0));
  for ( i=upp+1; i<N  ; i++ ) Ew[i] = BBE_CMPLXF32(H[HIDX(i,i)], XT_CONST_S(0));

  /* Loop until all eogenvalues are found, or failed to converge. */
  t = 0.f; its = 0;
  en = upp; na = en-1;
  while (en>=low) {
    /* Look for a negligible element on the subdiagonal, starting from
     * (en,en-1) upwards. If found at l, then the iteration will be performed
     * for submatrix T(l:en,l:en), othEwise l=low. Negligibility criterion
     * defined by [1] (23). */
#if 1
    stp = N-low;
    z = H[HIDX(low,low)];
    H0r = (xtfloat*)&H[HIDX(low+1,low+1)];
    for ( l=low, i=low+1; i<=en; i++ ) {
      x = xtfloat_loadi(H0r, -(int)sz_f32);   /* H[HIDX(i,i-1)]   */
      y = z;                                  /* H[HIDX(i-1,i-1)] */
      xtfloat_loadxp(z, H0r, (stp--)*sz_f32); /* H[HIDX(i,i)]     */
      r = XT_ABS_S(x); 
      s = XT_MUL_S(EPS, XT_ADD_S(XT_ABS_S(y), XT_ABS_S(z)));
      XT_MOVNEZ(l, i, BBE_MOVAB1(XT_OLE_S(r,s)));
    }
#else
    for ( l=en; l>low; l-- ) {
      x = H[HIDX(l,l-1)];
      y = H[HIDX(l-1,l-1)];
      z = H[HIDX(l,l)];
      if ( fabsf(x) <= EPS*(fabsf(y)+fabsf(z)) ) break;
    }
#endif
    /* Analyse the current state of the working matrix. */
    if (l<na) {
      /* No eigenvalue revealed. We should either:
       * A) perform a QR iteration with regular double shift, or
       * B) perform a QR iteration with exceptional shift, or
       * C) declare a failure! */
      if ( its>0 && 0==(its%10) ) {
        /* Try exceptional shift [1] (19) */
        r = H[HIDX(en,en)]; t += r;
        stp = N+1-low;
        H0r = Hw = (xtfloat*)&H[HIDX(low,low)];
        __Pragma("loop_count min=3");
        for ( i=low; i<=en; i++ ) {
          xtfloat_loadxp(s, H0r, stp*sz_f32);
          s = XT_SUB_S(s,r);
          xtfloat_storexp(s, Hw, (stp--)*sz_f32);
        }
        s = XT_ADD_S(XT_ABS_S(H[HIDX(en,na)]), XT_ABS_S(H[HIDX(na,en-2)]));
        x = y = 0.75f*s; w = -0.4375f*s*s;
      } else if ( its<ITS_LIM ) {
        /* Apply the regular shift */
        x = H[HIDX(en,en)];
        y = H[HIDX(na,na)];
        w = H[HIDX(na,en)] * H[HIDX(en,na)]; 
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
#if 1
      {
        float32_t _p,_q,_r;
        float32_t h0,h1,h2,h3,h4,h5,h6,h7;

        m = l; 
        h0 = H[HIDX(m-1,m-1)]; h1 = H[HIDX(m  ,m  )];
        h2 = H[HIDX(m+1,m+1)]; h3 = H[HIDX(m+2,m+2)];
        h4 = H[HIDX(m  ,m-1)]; h5 = H[HIDX(m+1,m  )];
        h6 = H[HIDX(m+2,m+1)]; h7 = H[HIDX(m  ,m+1)];

        z = h1; r = x-z; s = y-z;
        /* p = (r*s-w)/h5 + h7; */
        a = XT_NEG_S(w); XT_MADD_S(a,r,s);
        p = XT_ADD_S(IT_FDIVF32(a,h5,0),h7);
        /* q = h2-z-r-s; */
        q = XT_SUB_S(XT_SUB_S(XT_SUB_S(h2,z),r),s);
        r = h6;
        s = XT_ADD_S(XT_ABS_S(p), XT_ADD_S(XT_ABS_S(q), XT_ABS_S(r)));
        /* p = p/s; q = q/s; r = r/s; */
        a = XT_RECIP0_S(s);
        b = XT_CONST_S(1); XT_MSUB_S(b,a,s); XT_MADD_S(a,b,a);
        b = p; p = XT_MUL_S(a,p); XT_MSUB_S(b,s,p); XT_MADD_S(p,a,b);
        b = q; q = XT_MUL_S(a,q); XT_MSUB_S(b,s,q); XT_MADD_S(q,a,b);
        b = r; r = XT_MUL_S(a,r); XT_MSUB_S(b,s,r); XT_MADD_S(r,a,b);

        H0r = (xtfloat*)&H[HIDX(l+3,l+2)];
        H1r = (xtfloat*)&H[HIDX(l+1,l+2)]; 
        stp = N-l;
        for ( i=l+1; i<na; i++ ) {
          h0 = h1; h1 = h2; h2 = h3; 
          h3 = xtfloat_loadi(H0r, sz_f32);
          h4 = h5; h5 = h6; 
          xtfloat_loadxp(h6, H0r, (stp-2)*sz_f32);
          xtfloat_loadxp(h7, H1r, (stp--)*sz_f32);

          z = h1; _r = x-z; s = y-z;
          /* p = (r*s-w)/h5 + h7; */
          a = XT_NEG_S(w); XT_MADD_S(a,_r,s);
          _p = XT_ADD_S(IT_FDIVF32(a,h5,0),h7);
          /* q = h2-z-r-s; */
          _q = XT_SUB_S(XT_SUB_S(XT_SUB_S(h2,z),_r),s);
          _r = h6;
          s = XT_ADD_S(XT_ABS_S(_p), XT_ADD_S(XT_ABS_S(_q), XT_ABS_S(_r)));
          /* p = p/s; q = q/s; r = r/s; */
          a = XT_RECIP0_S(s);
          b = XT_CONST_S(1); XT_MSUB_S(b,a,s); XT_MADD_S(a,b,a);
          b = _p; _p = XT_MUL_S(a,_p); XT_MSUB_S(b,s,_p); XT_MADD_S(_p,a,b);
          b = _q; _q = XT_MUL_S(a,_q); XT_MSUB_S(b,s,_q); XT_MADD_S(_q,a,b);
          b = _r; _r = XT_MUL_S(a,_r); XT_MSUB_S(b,s,_r); XT_MADD_S(_r,a,b);

          s = XT_ADD_S(XT_ABS_S(_q), XT_ABS_S(_r)); 
          a = XT_MUL_S(XT_ABS_S(h4),s);
          s = XT_ADD_S(XT_ADD_S(XT_ABS_S(h0), XT_ABS_S(z)), XT_ABS_S(h2));
          b = XT_MUL_S(XT_ABS_S(_p),s);
          ble = XT_OLE_S(a, XT_MUL_S(EPS,b));
          XT_MOVNEZ(m,i,BBE_MOVAB1(ble));
          p = BBE_MOVF32T(_p,p,ble);
          q = BBE_MOVF32T(_q,q,ble);
          r = BBE_MOVF32T(_r,r,ble);
        } /* i */
      }
#elif 0
      {
        float32_t _p,_q,_r;
        float32_t h0,h1,h2,h3,h4,h5,h6,h7;
        m = l; 
        h0 = H[HIDX(m-1,m-1)]; h1 = H[HIDX(m  ,m  )];
        h2 = H[HIDX(m+1,m+1)]; h3 = H[HIDX(m+2,m+2)];
        h4 = H[HIDX(m  ,m-1)]; h5 = H[HIDX(m+1,m  )];
        h6 = H[HIDX(m+2,m+1)]; h7 = H[HIDX(m  ,m+1)];
        z = h1; r = x-z; s = y-z;
        p = (r*s-w)/h5 + h7;
        q = h2-z-r-s;
        r = h6;
        s = fabsf(p) + fabsf(q) + fabsf(r);
        p = p/s; q = q/s; r = r/s;
        for ( i=l+1; i<na; i++ ) {
          h0 = h1; h1 = h2; h2 = h3; h3 = H[HIDX(i+2,i+2)];
          h4 = h5; h5 = h6; h6 = H[HIDX(i+2,i+1)];
          h7 = H[HIDX(i,i+1)];
          z = h1; _r = x-z; s = y-z;
          _p = (_r*s-w)/h5 + h7;
          _q = h2-z-_r-s;
          _r = h6;
          s = fabsf(_p) + fabsf(_q) + fabsf(_r);
          _p = _p/s; _q = _q/s; _r = _r/s;
          a = fabsf(h4)*(fabsf(_q) + fabsf(_r));
          b = fabsf(_p)*(fabsf(h0) + fabsf(z) + fabsf(h2));
          if (a<=EPS*b) {
            m = i; p = _p; q = _q; r = _r;
          }
        }
      }
#else
      m = na;
      do {
        m--; z = H[HIDX(m,m)]; r = x-z; s = y-z;
        p = (r*s-w)/H[HIDX(m+1,m)] + H[HIDX(m,m+1)];
        q = H[HIDX(m+1,m+1)]-z-r-s;
        r = H[HIDX(m+2,m+1)];
        s = fabsf(p) + fabsf(q) + fabsf(r);
        p = p/s; q = q/s; r = r/s;
        if (m>l) {
          a = fabsf(H[HIDX(m,m-1)]) * (fabsf(q) + fabsf(r));
          b = fabsf(p) * (fabsf(H[HIDX(m-1,m-1)]) + fabsf(z) +
                          fabsf(H[HIDX(m+1,m+1)]));
          if (a<=EPS*b) break;
        }
      } while (m>l);
#endif
      /* Perform a double QR step involving rows 1 to n and columns m to n. */
      /* !!!! */
#if 0
      printf("its:%d l:%d m:%d en:%d p:%.7e q:%.7e r:%.7e\n",
             its,l,m,en,(float64_t)p,(float64_t)q,(float64_t)r);
#endif
      /* !!!! */
      francisQRStep(H,P,p,q,r,l,m,en,low,upp,N); its++;
    } else if (l==na) {
      /* Revealed a pair of eigenvalues. */
      x = H[HIDX(en,en)];
      y = H[HIDX(na,na)];
      w = H[HIDX(na,en)] * H[HIDX(en,na)];
      p = (y-x)/2; q = p*p+w; 
      z = XT_SQRT_S(XT_ABS_S(q));
      H[HIDX(en,en)] = x = x+t;
      H[HIDX(na,na)] = y = y+t;
      if (BBE_MOVAB1(XT_OLE_S(XT_CONST_S(0),q))) {
        /* Pair of real eigenvalues */
        XT_ADDSUB_S(a,b,p,z); 
        bltz = XT_OLT_S(p, XT_CONST_S(0));
        z = BBE_MOVF32T(a,b,bltz);
        beqz = XT_UEQ_S(w, XT_CONST_S(0));
        a = XT_ADD_S(x,z); b = XT_SUB_S(x, IT_FDIVF32(w,z,1)); 
        Ew[na] = BBE_CMPLXF32(XT_CONST_S(0), BBE_MOVF32T(y,a,beqz));
        Ew[en] = BBE_CMPLXF32(XT_CONST_S(0), BBE_MOVF32T(x,b,beqz));
        if (P) {
          /* Annihilate T(en,en-1) entry by similarity transformation involving 
           * a rotation in the en,en-1 plane. */
          x = H[HIDX(en,na)];
          /* r = sqrtf(x*x+z*z); */
          s = XT_MUL_S(x,x); XT_MADD_S(s,z,z); r = XT_SQRT_S(s);
          /* p = x/r; q = z/r; */
          a = XT_RECIP0_S(r);
          b = XT_CONST_S(1); XT_MSUB_S(b,a,r); XT_MADD_S(a,b,a);
          p = XT_MUL_S(x,a); XT_MSUB_S(x,p,r); XT_MADD_S(p,x,a);
          q = XT_MUL_S(z,a); XT_MSUB_S(z,q,r); XT_MADD_S(q,z,a);
          /* Row modification */
          H0r = Hw = (xtfloat*)&H[HIDX(en,na)];
          __Pragma("loop_count min=1");
          for ( j=na; j<N; j++ ) {
            /* x = H[HIDX(na,j)];
             * H[HIDX(na,j)] = q*x + p*H[HIDX(en,j)];
             * H[HIDX(en,j)] = q*H[HIDX(en,j)] - p*x; */
            x = xtfloat_loadx(H0r, (na-N)*sz_f32);
            xtfloat_loadip(y, H0r, sz_f32);
            z = XT_MUL_S(q,x); XT_MADD_S(z,p,y);
            w = XT_MUL_S(q,y); XT_MSUB_S(w,p,x);
            xtfloat_storex (z, Hw, (na-N)*sz_f32);
            xtfloat_storeip(w, Hw, sz_f32);
          }
          __Pragma("no_reorder");
          /* Column modification */
          H0r = Hw = (xtfloat*)&H[HIDX(0,en)];
          stp = N;
          __Pragma("loop_count min=1");
          for ( i=0; i<=en; i++ ) {
            /* x = H[HIDX(i,na)];
             * H[HIDX(i,na)] = q*x + p*H[HIDX(i,en)];
             * H[HIDX(i,en)] = q*H[HIDX(i,en)] - p*x; */
            x = xtfloat_loadi(H0r, -(int)sz_f32);
            xtfloat_loadxp(y, H0r, stp*sz_f32);
            z = XT_MUL_S(q,x); XT_MADD_S(z,p,y);
            w = XT_MUL_S(q,y); XT_MSUB_S(w,p,x);
            xtfloat_storei (z, Hw, -(int)sz_f32);
            xtfloat_storexp(w, Hw, (stp--)*sz_f32);
          }
          /* Accumulate transformations: right-multiply P by the rotator. */
          Pr = Pw = (xtfloat*)&P[low*N+en];
          __Pragma("loop_count min=1");
          for ( i=low; i<=upp; i++ ) {
            /* x = P[i*N+na]; 
             * P[i*N+na] = q*x + p*P[i*N+en];
             * P[i*N+en] = q*P[i*N+en] - p*x; */
            x = xtfloat_loadi(Pr, -(int)sz_f32);
            xtfloat_loadxp(y, Pr, N*sz_f32);
            z = XT_MUL_S(q,x); XT_MADD_S(z,p,y);
            w = XT_MUL_S(q,y); XT_MSUB_S(w,p,x);
            xtfloat_storei (z, Pw, -(int)sz_f32);
            xtfloat_storexp(w, Pw, N*sz_f32);
          }
        } /* P */
      } else {
        /* Complex eigenvalue and its conjugate */
        Ew[na] = BBE_CMPLXF32( z,x+p);
        Ew[en] = BBE_CMPLXF32(-z,x+p);
      } /* q */
      en = na-1; na = en-1; its = 0;
    } else {
      /* Revealed a standalone eigenvalue. */
      H[HIDX(en,en)] += t;
      Ew[en] = BBE_CMPLXF32(XT_CONST_S(0), H[HIDX(en,en)]);
      en = na; na = en-1; its = 0;
    }
  } /* en */

  return (1);
#else
  int i,j,l,m,en,na,its;
  float32_t a,b,p,q,r,s,t;
  float32_t x,y,z,w;

  /* Extract the isolated roots. */
  for ( i=0    ; i<low; i++ ) e[i] = _makecomplexf(H[HIDX(i,i)], 0.f);
  for ( i=upp+1; i<N  ; i++ ) e[i] = _makecomplexf(H[HIDX(i,i)], 0.f);

  /* Loop until all eogenvalues are found, or failed to converge. */
  t = 0.f; its = 0;
  en = upp; na = en-1;
  while (en>=low) {
    /* Look for a negligible element on the subdiagonal, starting from
     * (en,en-1) upwards. If found at l, then the iteration will be performed
     * for submatrix T(l:en,l:en), othEwise l=low. Negligibility criterion
     * defined by [1] (23). */
    for ( l=en; l>low; l-- ) {
      x = H[HIDX(l,l-1)];
      y = H[HIDX(l-1,l-1)];
      z = H[HIDX(l,l)];
      if ( fabsf(x) <= EPS*(fabsf(y)+fabsf(z)) ) break;
    }
    /* Analyse the current state of the working matrix. */
    if (l<na) {
      /* No eigenvalue revealed. We should either:
       * A) perform a QR iteration with regular double shift, or
       * B) perform a QR iteration with exceptional shift, or
       * C) declare a failure! */
      if ( its>0 && 0==(its%10) ) {
        /* Try exceptional shift [1] (19) */
        r = H[HIDX(en,en)]; t += r;
        for ( i=low; i<=en; i++ ) H[HIDX(i,i)] -= r;
        s = fabsf(H[HIDX(en,na)]) + fabsf(H[HIDX(na,en-2)]);
        x = y = 0.75f*s; w = -0.4375f*s*s;
      } else if ( its<ITS_LIM ) {
        /* Apply the regular shift */
        x = H[HIDX(en,en)];
        y = H[HIDX(na,na)];
        w = H[HIDX(na,en)] * H[HIDX(en,na)]; 
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
        m--; z = H[HIDX(m,m)]; r = x-z; s = y-z;
        p = (r*s-w)/H[HIDX(m+1,m)] + H[HIDX(m,m+1)];
        q = H[HIDX(m+1,m+1)]-z-r-s;
        r = H[HIDX(m+2,m+1)];
        s = fabsf(p) + fabsf(q) + fabsf(r);
        p = p/s; q = q/s; r = r/s;
        if (m>l) {
          a = fabsf(H[HIDX(m,m-1)]) * (fabsf(q) + fabsf(r));
          b = fabsf(p) * (fabsf(H[HIDX(m-1,m-1)]) + fabsf(z) +
                          fabsf(H[HIDX(m+1,m+1)]));
          if (a<=EPS*b) break;
        }
      } while (m>l);
      /* Perform a double QR step involving rows 1 to n and columns m to n. */
      francisQRStep(H,P,p,q,r,l,m,en,low,upp,N); its++;
    } else if (l==na) {
      /* Revealed a pair of eigenvalues. */
      x = H[HIDX(en,en)];
      y = H[HIDX(na,na)];
      w = H[HIDX(na,en)] * H[HIDX(en,na)];
      p = (y-x)/2; q = p*p+w; z = sqrtf(fabsf(q));
      H[HIDX(en,en)] = x = x+t;
      H[HIDX(na,na)] = y = y+t;
      if (q>=0.f) {
        /* Pair of real eigenvalues */
        z = ( p<0.f ? p-z : p+z );
        e[na] = _makecomplexf( w!=0.f ? x+z   : y, 0.f ); 
        e[en] = _makecomplexf( w!=0.f ? x-w/z : x, 0.f );
        if (P) {
          /* Annihilate T(en,en-1) entry by similarity transformation involving 
           * a rotation in the en,en-1 plane. */
          x = H[HIDX(en,na)];
          r = sqrtf(x*x+z*z); p = x/r; q = z/r;
          /* Row modification */
          for ( j=na; j<N; j++ ) {
            z = H[HIDX(na,j)];
            H[HIDX(na,j)] = q*z + p*H[HIDX(en,j)];
            H[HIDX(en,j)] = q*H[HIDX(en,j)] - p*z;
          }
          /* Column modification */
          for ( i=0; i<=en; i++ ) {
            z = H[HIDX(i,na)];
            H[HIDX(i,na)] = q*z + p*H[HIDX(i,en)];
            H[HIDX(i,en)] = q*H[HIDX(i,en)] - p*z;
          }
          /* Accumulate transformations: right-multiply P by the rotator. */
          for ( i=low; i<=upp; i++ ) {
            z = P[i*N+na]; 
            P[i*N+na] = q*z + p*P[i*N+en];
            P[i*N+en] = q*P[i*N+en] - p*z;
          }
        }
      } else {
        /* Complex eigenvalue and its conjugate */
        e[na] = _makecomplexf( x+p, z );
        e[en] = _makecomplexf( x+p, -z );
      }
      en = na-1; na = en-1; its = 0;
    } else {
      /* Revealed a standalone eigenvalue. */
      H[HIDX(en,en)] += t;
      e[en] = _makecomplexf( H[HIDX(en,en)], 0.f );
      en = na; na = en-1; its = 0;
    }
  } /* en */

  return (1);
#endif
} /* reigen_hqr_f() */

/* Perform a Francis QR step for real upper-Hessenberg matrix, and
 * optionally update the transform matrix. */
void francisQRStep ( float32_t * restrict H,
                     float32_t * restrict P,
                     float32_t p, float32_t q, float32_t r,
                     int l, int m, int en, int low, int upp, int N )
{
  /*
   * MATLAB reference code:
   *
   *   function [H,P] = FrancisQRStep(H,P,l,m,en,na,n,low,upp,p,q,r)
   *   assert((l<=m)&&(m<=en-2)&&(low+2<=en)&&(en<=upp));
   *   RMIN = realmin(class(H));
   *   t = 1; g = 0;
   *   for k = m:na
   *     notLast = k<na;
   *     if t<RMIN, break; end;
   *     p = p/t; q = q/t; r = r/t;
   *     s = sqrt(p^2+q^2+r^2);
   *     if p<0, s = -s; end;
   *     if k>m
   *       H(k,k-1) = -s*t;
   *     elseif l<m
   *       H(k,k-1) = -H(k,k-1);
   *     end
   *     x = p/s+1; y = q/s; z = r/s;
   *     q = q/(p+s); r = r/(p+s);
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
   *       t = g-p;
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
   *     q = t; t = abs(p)+abs(q)+abs(r);
   *   end
   */
#if 1
  const xtfloat * restrict Pr;
  const xtfloat * restrict H0r;
  const xtfloat * restrict H1r;
  const xtfloat * restrict H2r;
  xtfloat * restrict Pw;
  xtfloat * restrict H0w;
  xtfloat * restrict H1w;
  xtfloat * restrict H2w;
  float32_t a,b,c,e,g,h,s,t,x,y,z;
  int i,j,k;
  int stp;
  int na = en-1;

  NASSERT( (l<=m) && (m<=en-2) && ((low+2)<=en) && (en<=upp) && (upp<N) );

  t = XT_CONST_S(1); g = XT_CONST_S(0);
  for ( k=m; k<en && t>=FLT_MIN; k++ ) {
    /* Construct a Householder reflector from the vector (p,q,r). */
    if (k>m) { 
      /* p /= t; q /= t; r /= t; */
      h = XT_RECIP0_S(t);
      e = XT_CONST_S(1); XT_MSUB_S(e,t,h); XT_MADD_S(h,e,h);
      e = p; p = XT_MUL_S(h,p); XT_MSUB_S(e,t,p); XT_MADD_S(p,e,h);
      e = q; q = XT_MUL_S(h,q); XT_MSUB_S(e,t,q); XT_MADD_S(q,e,h);
      e = r; r = XT_MUL_S(h,r); XT_MSUB_S(e,t,r); XT_MADD_S(r,e,h);
    }
    /* s = sqrtf(p*p+q*q+r*r); */
    h = XT_MUL_S(p,p); XT_MADD_S(h,q,q); XT_MADD_S(h,r,r);
    s = XT_SQRT_S(h);
    /* if (p<0) s = -s; */
    BBE_NEGF32T(s, s, XT_OLT_S(p, XT_CONST_S(0)));
    if (k>m) {
      H[HIDX(k,k-1)] = -s*t;
    } else if (l<m) {
      H[HIDX(k,k-1)] = -H[HIDX(k,k-1)];
    }
    /* x = p/s+1.f; y = q/s;  z = r/s; */
    h = XT_RECIP0_S(s); e = XT_CONST_S(1); XT_MSUB_S(e,s,h); XT_MADD_S(h,e,h);
    x = XT_MUL_S(h,p); e = p; XT_MSUB_S(e,s,x); XT_MADD_S(x,e,h);
    y = XT_MUL_S(h,q); e = q; XT_MSUB_S(e,s,y); XT_MADD_S(y,e,h);
    z = XT_MUL_S(h,r); e = r; XT_MSUB_S(e,s,z); XT_MADD_S(z,e,h);
    x = XT_ADD_S(x, XT_CONST_S(1)); p = XT_ADD_S(p,s);
    /* q = q/(p+s); r = r/(p+s); */
    h = XT_RECIP0_S(p); e = XT_CONST_S(1); XT_MSUB_S(e,p,h); XT_MADD_S(h,e,h);
    e = q; q = XT_MUL_S(h,q); XT_MSUB_S(e,p,q); XT_MADD_S(q,e,h);
    e = r; r = XT_MUL_S(h,r); XT_MSUB_S(e,p,r); XT_MADD_S(r,e,h);
    /*----------------------------------------------------------------*
     * Accumulate transformations: right-multiply P by the reflector. *
     *----------------------------------------------------------------*/
    if (P) {
      if (k<na) {
        Pr = Pw = (xtfloat*)&P[low*N+k+2];
        __Pragma("loop_count min=1");
        __Pragma("concurrent");
        for ( i=low; i<=upp; i++ ) {
          a = xtfloat_loadi(Pr, -2*(int)sz_f32);
          b = xtfloat_loadi(Pr, -1*(int)sz_f32);
          xtfloat_loadxp(c, Pr, N*sz_f32);
          p = XT_MUL_S(a,x); XT_MADD_S(p,b,y); XT_MADD_S(p,c,z);
          b = xtfloat_loadx(Pr, -(N+1)*sz_f32);
          a = XT_SUB_S(a,p); XT_MSUB_S(b,p,q); XT_MSUB_S(c,p,r);
          xtfloat_storei(a, Pw, -2*(int)sz_f32);
          xtfloat_storei(b, Pw, -1*(int)sz_f32);
          xtfloat_storexp(c, Pw, N*sz_f32);
        }
      } else { /* k */
        Pr = Pw = (xtfloat*)&P[low*N+k+1];
        __Pragma("loop_count min=1");
        for ( i=low; i<=upp; i++ ) {
          a = xtfloat_loadi(Pr, -1*(int)sz_f32);
          xtfloat_loadxp(b, Pr, N*sz_f32);
          p = XT_MUL_S(a,x); XT_MADD_S(p,b,y);
          a = XT_SUB_S(a,p); XT_MSUB_S(b,p,q);
          xtfloat_storei(a, Pw, -1*(int)sz_f32);
          xtfloat_storexp(b, Pw, N*sz_f32);
        }
      } /* k */
    } /* P */
    /*----------------------------------------------------------------*
     * Row modification, involves columns k through N-1.              *
     *----------------------------------------------------------------*/
    /* k-th column */
    p = H[HIDX(k,k)] + q*H[HIDX(k+1,k)] + r*g;
    g -= p*z; H[HIDX(k+1,k)] -= p*y;
    H[HIDX(k,k)] -= p*x;
    /* Columns k+1..N-1 or k+1..en */
    i = ( P ? N-1 : en );
    if (k<na) {
      H0r = H0w = (xtfloat*)&H[HIDX(k  ,k+1)];
      H1r = H1w = (xtfloat*)&H[HIDX(k+1,k+1)];
      H2r = H2w = (xtfloat*)&H[HIDX(k+2,k+1)];
      __Pragma("loop_count min=1");
      for ( j=k+1; j<=i; j++ ) {
        xtfloat_loadip(a, H0r, sz_f32);
        xtfloat_loadip(b, H1r, sz_f32);
        xtfloat_loadip(c, H2r, sz_f32);
        p = a; XT_MADD_S(p,q,b); XT_MADD_S(p,r,c);
        XT_MSUB_S(a,x,p); XT_MSUB_S(b,y,p); XT_MSUB_S(c,z,p);
        xtfloat_storeip(a, H0w, sz_f32);
        xtfloat_storeip(b, H1w, sz_f32);
        xtfloat_storeip(c, H2w, sz_f32);
      }
    } else { /* k */
      H0r = H0w = (xtfloat*)&H[HIDX(k  ,k+1)];
      H1r = H1w = (xtfloat*)&H[HIDX(k+1,k+1)];
      __Pragma("loop_count min=1");
      for ( j=k+1; j<=i; j++ ) {
        xtfloat_loadip(a, H0r, sz_f32);
        xtfloat_loadip(b, H1r, sz_f32);
        p = a; XT_MADD_S(p,q,b);
        XT_MSUB_S(a,x,p); XT_MSUB_S(b,y,p);
        xtfloat_storeip(a, H0w, sz_f32);
        xtfloat_storeip(b, H1w, sz_f32);
      }
    } /* k */
    __Pragma("no_reorder");

    /*----------------------------------------------------------------*
     * Column modification. Involves rows 0..k+3 (update of P is      *
     * requested) or l-4..k+3 (no need for transform accumulation)    *
     *----------------------------------------------------------------*/
    /* Rows 0 (or l-4) through k+1 */
    j = ( P || (l<4) ? 0 : l-4 );
    if (k<na) {
      stp = (N-j);
      H0r = H0w = (xtfloat*)&H[HIDX(j,k+2)];
      __Pragma("loop_count min=1");
      for ( i=j; i<=k+1; i++ ) {
        a = xtfloat_loadi(H0r, -2*(int)sz_f32);
        b = xtfloat_loadi(H0r, -1*(int)sz_f32);
        xtfloat_loadxp(c, H0r, stp*sz_f32);
        p = XT_MUL_S(a,x); XT_MADD_S(p,b,y); XT_MADD_S(p,c,z);
        a = XT_SUB_S(a,p); XT_MSUB_S(b,p,q); XT_MSUB_S(c,p,r);
        xtfloat_storei(a, H0w, -2*(int)sz_f32);
        xtfloat_storei(b, H0w, -1*(int)sz_f32);
        xtfloat_storexp(c, H0w, (stp--)*sz_f32);
      }
    } else { /* k */
      stp = (N-j);
      H0r = H0w = (xtfloat*)&H[HIDX(j,k+1)];
      __Pragma("loop_count min=1");
      for ( i=j; i<=k+1; i++ ) {
        a = xtfloat_loadi(H0r, -1*(int)sz_f32);
        xtfloat_loadxp(b, H0r, stp*sz_f32);
        p = XT_MUL_S(a,x); XT_MADD_S(p,b,y);
        a = XT_SUB_S(a,p); XT_MSUB_S(b,p,q);
        xtfloat_storei(a, H0w, -1*(int)sz_f32);
        xtfloat_storexp(b, H0w, (stp--)*sz_f32);
      }
    } /* k */

    /* Row k+2 */
    if (k<na) {
      p = x*g + y*H[HIDX(k+2,k+1)] + z*H[HIDX(k+2,k+2)];
      H[HIDX(k+2,k+2)] -= p*r;
      H[HIDX(k+2,k+1)] -= p*q;
      t = g-p;
    }
    /* Row k+3 */
    if (k<na-1) {
      p = z*H[HIDX(k+3,k+2)];
      H[HIDX(k+3,k+2)] -= p*r;
      g = -p*q; r = -p;
    } else {
      r = XT_CONST_S(0);
    }

    __Pragma("no_reorder");
    p = H[HIDX(k+1,k)];
    q = t; t = XT_ADD_S(XT_ABS_S(p), XT_ADD_S(XT_ABS_S(q), XT_ABS_S(r)));
  } /* k */
#else
  float32_t x,y,z,f,g,s,t;
  int i,j,k;
  int na = en-1;

  NASSERT( (l<=m) && (m<=en-2) && ((low+2)<=en) && (en<=upp) && (upp<N) );

  t = 1.f; g = 0.f;
  for ( k=m; k<en && t!=0.f; k++ ) {
    /* Construct a Householder reflector from the vector (p,q,r). */
    if (k>m) { p /= t; q /= t; r /= t; };
    s = sqrtf(p*p+q*q+r*r);
    if (p<0) s = -s;
    if (k>m) {
      H[HIDX(k,k-1)] = -s*t;
    } else if (l<m) {
      H[HIDX(k,k-1)] = -H[HIDX(k,k-1)];
    }
    x = p/s+1.f; y = q/s; z = r/s;
    q = q/(p+s); r = r/(p+s);
    /*----------------------------------------------------------------*
     * Accumulate transformations: right-multiply P by the reflector. *
     *----------------------------------------------------------------*/
    if (P) {
      for ( i=low; i<=upp; i++ ) {
        p = x*P[i*N+k] + y*P[i*N+k+1];
        if (k<na) {
          p += z*P[i*N+k+2];
          P[i*N+k+2] -= p*r;
        }
        P[i*N+k+1] -= p*q; P[i*N+k] -= p;
      }
    }
    /*----------------------------------------------------------------*
     * Row modification, involves columns k through N-1.              *
     *----------------------------------------------------------------*/
    /* k-th column */
    p = H[HIDX(k,k)] + q*H[HIDX(k+1,k)] + r*g;
    g -= p*z; H[HIDX(k+1,k)] -= p*y;
    H[HIDX(k,k)] -= p*x;
    /* Columns k+1..N-1 or k+1..en */
    i = ( P ? N-1 : en );
    for ( j=k+1; j<=i; j++ ) {
      p = H[HIDX(k,j)] + q*H[HIDX(k+1,j)];
      if (k<na) {
        p += r*H[HIDX(k+2,j)];
        H[HIDX(k+2,j)] -= p*z;
      }
      H[HIDX(k+1,j)] -= p*y;
      H[HIDX(k,j)] -= p*x;
    }
    /*----------------------------------------------------------------*
     * Column modification. Involves rows 0..k+3 (update of P is      *
     * requested) or l-3..k+3 (no need for transform accumulation)    *
     *----------------------------------------------------------------*/
    /* Rows 0 (or l-3) through k+1 */
    j = ( P || (l<3) ? 0 : l-3 );
    for ( i=j; i<=k+1; i++ ) {
      p = x*H[HIDX(i,k)] + y*H[HIDX(i,k+1)];
      if (k<na) {
        p += z*H[HIDX(i,k+2)];
        H[HIDX(i,k+2)] -= p*r;
      }
      H[HIDX(i,k+1)] -= p*q;
      H[HIDX(i,k)] -= p;
    }
    /* Row k+2 */
    if (k<na) {
      p = x*g + y*H[HIDX(k+2,k+1)] + z*H[HIDX(k+2,k+2)];
      H[HIDX(k+2,k+2)] -= p*r;
      H[HIDX(k+2,k+1)] -= p*q;
      f = g-p;
    }
    /* Row k+3 */
    if (k<na-1) {
      p = z*H[HIDX(k+3,k+2)];
      H[HIDX(k+3,k+2)] -= p*r;
      g = -p*q; r = -p;
    } else {
      r = 0.f;
    }
    p = H[HIDX(k+1,k)];
    q = f; t = fabsf(p) + fabsf(q) + fabsf(r);
  } /* k */
#endif
} /* francisQRStep() */

#endif /* HAVE_VFPU */

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
    Hessenberg Form Of Real Matrices, Block and Stream Order
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
#include "vfpu_math.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if HAVE_VFPU

#define sz_f32    sizeof(float32_t)

/* Reduction of a single matrix to upper-Hessenberg form. Reduced matrix overwrites 
 * input data. Both the block and stream order are supported.
 * This function is based on orthes() algorithm from
 *   "Similarity Reduction of a General Matrix to Hessenberg Form" by R.S. Martin and
 *   J.H. Wilkinson, Handbook for Automatie Computation, Vol.II Linear Algebra, 
 *   Contribution II/13. */
#if 0
static void rhessf ( float32_t * restrict P, 
                     float32_t * restrict A, 
                     int N, int stride )
{
  /*
   * MATLAB reference code and a simple test:
   *
   *   n = 5; A = 1-2*rand(n,n); B = A;
   *   % Reduce the matrix to upper-Hessenberg form
   *   P = eye(n);
   *   tol = realmin(class(A));
   *   for m=2:n-1
   *     u = A(m:n,m-1);
   *     h = u'*u;
   *     if h>tol
   *       % Construct and apply a Householder reflector.
   *       g = sqrt(h);
   *       if u(1)>0, g = -g; end;
   *       h = h-g*u(1);
   *       u(1) = u(1)-g;
   *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
   *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
   *       % Accumulate transformations.
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

  float32_t f,g,h,s;
  int m,n,p;
  /* Initialize P with identity matrix. */
  if (P) {
    for ( m=0; m<N; m++ ) {
      for ( n=0; n<N; n++ ) {
        P[(m*N+n)*stride] = ( m==n ? 1.f : 0.f );
      }
    }
  }
  /* Perform N-2 Householder reflections. */
  for ( m=1; m<N-1; m++ ) {
    /* Squared norm of the pivot column */
    for ( h=0.f, n=m; n<N; n++ ) { f = A[(n*N+m-1)*stride]; h += f*f; }
    /* Skip the transformation if the norm is too small to 
     * be evaluated to working precision. */
    if (h<=FLT_MIN) { A[(m*(N+1)-1)*stride] = 0.f; continue; }
    /* In-place construction of the first column of Householder matrix. */
    f = A[(m*(N+1)-1)*stride];
    g = ( f>0 ? -sqrtf(h) : sqrtf(h) );
    h -= g*f;
    A[(m*(N+1)-1)*stride] = f-g;
    /* Left-hand transformation */
    for ( n=m; n<N; n++ ) {
      for ( s=0.f, p=m; p<N; p++ ) {
        s += A[(p*N+n)*stride] * A[(p*N+m-1)*stride];
      }
      for ( p=m; p<N; p++ ) {
        A[(p*N+n)*stride] -= (s/h)*A[(p*N+m-1)*stride];
      }
    }
    /* Right-hand transformation */
    for ( n=0; n<N; n++ ) {
      for ( s=0.f, p=m; p<N; p++ ) {
        s += A[(n*N+p)*stride] * A[(p*N+m-1)*stride];
      }
      for ( p=m; p<N; p++ ) {
        A[(n*N+p)*stride] -= (s/h)*A[(p*N+m-1)*stride];
      }
    }
    /* Accumulate transformations */
    if (P) {
      for ( n=1; n<N; n++ ) {
        for ( s=0.f, p=m; p<N; p++ ) {
          s += P[(n*N+p)*stride] * A[(p*N+m-1)*stride];
        }
        for ( p=m; p<N; p++ ) {
          P[(n*N+p)*stride] -= (s/h)*A[(p*N+m-1)*stride];
        }
      }
    }
    A[(m*(N+1)-1)*stride] = g;
  } /* m */

} /* rhessf() */
#elif 0
static void rhessfVec(float32_t * restrict P,
                      float32_t * restrict A,
                      int N, int stride)
{
    /*
    * MATLAB reference code and a simple test:
    *
    *   n = 5; A = 1-2*rand(n,n); B = A;
    *   % Reduce the matrix to upper-Hessenberg form
    *   P = eye(n);
    *   tol = realmin(class(A));
    *   for m=2:n-1
    *     u = A(m:n,m-1);
    *     h = u'*u;
    *     if h>tol
    *       % Construct and apply a Householder reflector.
    *       g = sqrt(h);
    *       if u(1)>0, g = -g; end;
    *       h = h-g*u(1);
    *       u(1) = u(1)-g;
    *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
    *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
    *       % Accumulate transformations.
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
    const xb_vecN_2xf32 * A_r;
    const xb_vecN_2xf32 * A_r0;
    const xb_vecN_2xf32 * restrict A_r1;
    const xb_vecN_2xf32 * restrict P_r;
          xb_vecN_2xf32 * A_rw;
          xb_vecN_2xf32 * restrict A_w;
          xb_vecN_2xf32 * restrict P_w;

    const xb_vecN_2xf32 NORM_MIN = FLT_MIN;
    const xb_vecN_2xf32 c0 = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 c1 = BBE_CONSTN_2XF32(1);

    xb_vecN_2xf32 f, h, g, s, x;
    xb_vecN_2xf32 a0, a1;
    vboolN_2 bskip, bmz, beq;

    int m, n, p;
    
    if (P){
        P_w = (xb_vecN_2xf32*)P;

        for (n = 0; n < N*N; n++) {
            BBE_SVN_2XF32_XP(c0, P_w, stride*sz_f32);
        }

        P_w = (xb_vecN_2xf32*)P;
        for (n = 0; n < N; n++) {
            BBE_SVN_2XF32_XP(c1, P_w, (N+1)*stride*sz_f32);
        }
    }

    /* Perform N-2 Householder reflections. */
    for (m = 1; m < N - 1; m++) {
        /* Squared norm of the pivot column */
        h = BBE_ZERON_2XF32();
        A_r = (xb_vecN_2xf32*)&A[(m*N + m - 1)*stride];
        for (n = m; n<N; n++) {
            BBE_LVN_2XF32_XP(f, A_r, N*stride*sz_f32);
            BBE_MULAN_2XF32(h, f, f);
        }

        /* Skip the transformation if the norm is too small to
        * be evaluated to working precision. */
        A_rw = (xb_vecN_2xf32*)&A[(m*N + m - 1)*stride];

        bskip = BBE_OLEN_2XF32(h, NORM_MIN);
        BBE_SVN_2XF32T_I(c0, A_rw, 0, bskip);

        /* In-place construction of the first column of Householder matrix. */
        f = BBE_LVN_2XF32_I(A_rw, 0);
        bmz = BBE_OLTN_2XF32(c0, f);
        g = BBE_SQRTN_2XF32(h);
        BBE_NEGN_2XF32T(g, g, bmz);
        BBE_MULSN_2XF32(h, g, f);
        f = BBE_SUBN_2XF32(f, g);
        BBE_SVN_2XF32F_I(f, A_rw, 0, bskip);

        /* Left-hand transformation */
        for (n = 0; n < (N - m); n++) {
            s = BBE_ZERON_2XF32();
            A_r0 = (xb_vecN_2xf32 *)&A[(m*N + m + n)*stride];
            A_r1 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
            for (p = m; p < N; p++) {
                BBE_LVN_2XF32_XP(a0, A_r0, N*stride*sz_f32);
                BBE_LVN_2XF32_XP(a1, A_r1, N*stride*sz_f32);
                BBE_MULAN_2XF32(s, a0, a1);
            }
            s = IT_DIVN_2XF32(s, h);
            A_r0 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
            A_r1 = A_w = (xb_vecN_2xf32 *)&A[(m*N + m + n)*stride];
            for (p = m; p < N; p++) {
                BBE_LVN_2XF32_XP(a0, A_r0, N*stride*sz_f32);
                BBE_LVN_2XF32_XP(a1, A_r1, N*stride*sz_f32);
                BBE_MULSN_2XF32(a1, a0, s);
                BBE_SVN_2XF32F_XP(a1, A_w, N*stride*sz_f32, bskip);
            }
        }
        /* Right-hand transformation */
        for (n = 0; n < N; n++) {
            s = BBE_ZERON_2XF32();
            A_r0 = (xb_vecN_2xf32 *)&A[(n*N + m)*stride];
            A_r1 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
            for (p = m; p < N; p++) {
                BBE_LVN_2XF32_XP(a0, A_r0,   stride*sz_f32);
                BBE_LVN_2XF32_XP(a1, A_r1, N*stride*sz_f32);
                BBE_MULAN_2XF32(s, a0, a1);
            }
            s = IT_DIVN_2XF32(s, h);
            A_r0  = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
            A_r1 = A_w = (xb_vecN_2xf32 *)&A[(n*N + m)*stride];
            for (p = m; p < N; p++) {
                BBE_LVN_2XF32_XP(a0, A_r0, N*stride*sz_f32);
                BBE_LVN_2XF32_XP(a1, A_r1,   stride*sz_f32);
                BBE_MULSN_2XF32(a1, a0, s);
                BBE_SVN_2XF32F_XP(a1, A_w, stride*sz_f32, bskip);
            }
        }

        /* Accumulate transformations */
        if (P) {
            for (n = 1; n < N; n++) {
                s = BBE_ZERON_2XF32();
                A_r = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
                P_r = (xb_vecN_2xf32 *)&P[(n*N + m)*stride];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r, N*stride*sz_f32);
                    BBE_LVN_2XF32_XP(a1, P_r,   stride*sz_f32);
                    BBE_MULAN_2XF32(s, a0, a1);
                }
                s = IT_DIVN_2XF32(s, h);
                A_r = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
                P_r = P_w = (xb_vecN_2xf32 *)&P[(n*N + m)*stride];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r, N*stride*sz_f32);
                    BBE_LVN_2XF32_XP(a1, P_r, stride*sz_f32);
                    BBE_MULSN_2XF32(a1, a0, s);
                    BBE_SVN_2XF32F_XP(a1, P_w, stride*sz_f32, bskip);
                }
            }
        }

        A_w = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*stride];
        BBE_SVN_2XF32F_I(g, A_w, 0, bskip);
    }
} /* rhessfVec() */
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

/* Real-valued, Stream Order
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_hess_nxnsf ( 
                float32_t * restrict P,     /* P[N*N][L] */
                float32_t * restrict A,     /* A[N*N][L] */
                int N, int L )
#if 0
{
  int k;
  const int STEP = BBE_SIMD_WIDTH / 2;

  NASSERT_ALIGN( P, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
#if 0
  for ( k=0; k<L; k++ ) {
    rhessf( ( P ? P+k : 0 ), A+k, N, L );
  }
#else
  for (k = 0; k<(L/STEP); k++) {
      rhessfVec(P, A, N, L);
      A += STEP;
      if(P) P+= STEP;
  }
#endif
} /* reigen_hessnxnsf() */
#else
{
    /*
    * MATLAB reference code and a simple test:
    *
    *   n = 5; A = 1-2*rand(n,n); B = A;
    *   % Reduce the matrix to upper-Hessenberg form
    *   P = eye(n);
    *   tol = realmin(class(A));
    *   for m=2:n-1
    *     u = A(m:n,m-1);
    *     h = u'*u;
    *     if h>tol
    *       % Construct and apply a Householder reflector.
    *       g = sqrt(h);
    *       if u(1)>0, g = -g; end;
    *       h = h-g*u(1);
    *       u(1) = u(1)-g;
    *       A(m:n,m:n) = A(m:n,m:n)-u*(u'*A(m:n,m:n)/h);
    *       A(1:n,m:n) = A(1:n,m:n)-(A(1:n,m:n)*u/h)*u';
    *       % Accumulate transformations.
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

    const xb_vecN_2xf32 * A_r;
    const xb_vecN_2xf32 * A_r0;
    const xb_vecN_2xf32 * restrict A_r1;
    const xb_vecN_2xf32 * restrict P_r;
    xb_vecN_2xf32 * A_rw;
    xb_vecN_2xf32 * restrict A_w;
    xb_vecN_2xf32 * restrict P_w;

    const xb_vecN_2xf32 NORM_MIN = FLT_MIN;
    const xb_vecN_2xf32 c0 = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 c1 = BBE_CONSTN_2XF32(1);
    const int STEP = BBE_SIMD_WIDTH / 2;

    xb_vecN_2xf32 f, h, g, s;
    xb_vecN_2xf32 a0, a1;
    vboolN_2 bskip, bmz;

    int k, m, n, p;

    NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));
    
    if (P){
        const int step_norm = (N + 1)*L*sz_f32;
        const int step_back = STEP * sz_f32 - (N - 1)*(N + 1)*L*sz_f32;
        int step;
        P_w = (xb_vecN_2xf32*)P;
        for (n = 0; n < (N*N*L / STEP); n++) {
            BBE_SVN_2XF32_IP(c0, P_w, (BBE_SIMD_WIDTH/2)*sz_f32);
        }
#if 0
        for (k = 0; k < (L / STEP); k++){
            P_w = (xb_vecN_2xf32*)&P[k*STEP];
            for (n = 0; n < N; n++) {
                BBE_SVN_2XF32_XP(c1, P_w, (N + 1)*L*sz_f32);
            }
        }
#else
        P_w = (xb_vecN_2xf32*)P;
        n = 0;
        for (k = 0; k < (L / STEP)*N; k++) {
            n = BBE_ADDMOD16U(n, (N<<16) | 1);
            step = step_norm;
            XT_MOVEQZ(step, step_back, n);
            BBE_SVN_2XF32_XP(c1, P_w, step);
        }
#endif
    }

    for (k = 0; k < (L / STEP); k++){
        /* Perform N-2 Householder reflections. */
        for (m = 1; m < N - 1; m++) {
            /* Squared norm of the pivot column */
            h = BBE_ZERON_2XF32();
            A_r = (xb_vecN_2xf32*)&A[(m*N + m - 1)*L];
            for (n = m; n<N; n++) {
                BBE_LVN_2XF32_XP(f, A_r, N*L*sz_f32);
                BBE_MULAN_2XF32(h, f, f);
            }

            /* Skip the transformation if the norm is too small to
            * be evaluated to working precision. */
            A_rw = (xb_vecN_2xf32*)&A[(m*N + m - 1)*L];

            bskip = BBE_OLEN_2XF32(h, NORM_MIN);
            BBE_SVN_2XF32T_I(c0, A_rw, 0, bskip);

            /* In-place construction of the first column of Householder matrix. */
            f = BBE_LVN_2XF32_I(A_rw, 0);
            bmz = BBE_OLTN_2XF32(c0, f);
            g = BBE_SQRTN_2XF32(h);
            BBE_NEGN_2XF32T(g, g, bmz);
            BBE_MULSN_2XF32(h, g, f);
            f = BBE_SUBN_2XF32(f, g);
            BBE_SVN_2XF32F_I(f, A_rw, 0, bskip);

            /* Left-hand transformation */
            for (n = 0; n < (N - m); n++) {
                s = BBE_ZERON_2XF32();
                A_r0 = (xb_vecN_2xf32 *)&A[(m*N + m + n)*L];
                A_r1 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r0, N*L*sz_f32);
                    BBE_LVN_2XF32_XP(a1, A_r1, N*L*sz_f32);
                    BBE_MULAN_2XF32(s, a0, a1);
                }
                s = IT_FDIVN_2XF32(s, h, 0);
                A_r0 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                A_r1 = A_w = (xb_vecN_2xf32 *)&A[(m*N + m + n)*L];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r0, N*L*sz_f32);
                    BBE_LVN_2XF32_XP(a1, A_r1, N*L*sz_f32);
                    BBE_MULSN_2XF32(a1, a0, s);
                    BBE_SVN_2XF32F_XP(a1, A_w, N*L*sz_f32, bskip);
                }
            }
            /* Right-hand transformation */
            for (n = 0; n < N; n++) {
                s = BBE_ZERON_2XF32();
                A_r0 = (xb_vecN_2xf32 *)&A[(n*N + m)*L];
                A_r1 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r0,   L*sz_f32);
                    BBE_LVN_2XF32_XP(a1, A_r1, N*L*sz_f32);
                    BBE_MULAN_2XF32(s, a0, a1);
                }
                s = IT_FDIVN_2XF32(s, h, 0);
                A_r0 = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                A_r1 = A_w = (xb_vecN_2xf32 *)&A[(n*N + m)*L];
                for (p = m; p < N; p++) {
                    BBE_LVN_2XF32_XP(a0, A_r0, N*L*sz_f32);
                    BBE_LVN_2XF32_XP(a1, A_r1, L*sz_f32);
                    BBE_MULSN_2XF32(a1, a0, s);
                    BBE_SVN_2XF32F_XP(a1, A_w, L*sz_f32, bskip);
                }
            }

            /* Accumulate transformations */
            if (P) {
                for (n = 1; n < N; n++) {
                    s = BBE_ZERON_2XF32();
                    A_r = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                    P_r = (xb_vecN_2xf32 *)&P[(n*N + m)*L];
                    for (p = m; p < N; p++) {
                        BBE_LVN_2XF32_XP(a0, A_r, N*L*sz_f32);
                        BBE_LVN_2XF32_XP(a1, P_r, L*sz_f32);
                        BBE_MULAN_2XF32(s, a0, a1);
                    }
                    s = IT_FDIVN_2XF32(s, h, 0);
                    A_r = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
                    P_r = P_w = (xb_vecN_2xf32 *)&P[(n*N + m)*L];
                    for (p = m; p < N; p++) {
                        BBE_LVN_2XF32_XP(a0, A_r, N*L*sz_f32);
                        BBE_LVN_2XF32_XP(a1, P_r, L*sz_f32);
                        BBE_MULSN_2XF32(a1, a0, s);
                        BBE_SVN_2XF32F_XP(a1, P_w, L*sz_f32, bskip);
                    }
                }
            }

            A_w = (xb_vecN_2xf32 *)&A[(m*N + m - 1)*L];
            BBE_SVN_2XF32F_I(g, A_w, 0, bskip);
        }

        A += STEP;
        if (P) P += STEP;
    }
} /* reigen_hessnxnsf() */
#endif

#endif /* HAVE_VFPU */

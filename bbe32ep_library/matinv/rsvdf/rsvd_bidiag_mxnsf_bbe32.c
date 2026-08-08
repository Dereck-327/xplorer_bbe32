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
    Real Data, Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <float.h>
#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* SVD common declarations */
#include "svd_common.h"

#if HAVE_VFPU

#define sz_f32   sizeof(float32_t)

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

/* Real-valued, Stream Order
 * Note that output matrices U replace input matrices A.
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void rsvd_bidiag_mxnsf (
                float32_t * restrict D, /* D[N][L]   */
                float32_t * restrict F, /* F[N-1][L] */
                float32_t * restrict V, /* V[N*N][L] */
                float32_t * restrict A, /* A[M*N][L] */
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
   *   % Simple test matrix
   *   M = 7; N = 5;
   *   A = single(1-2*rand(M,N)); R = A;
   *   % Perform the bidiagonalization
   *   tol = realmin(class(A));
   *   D = cast(zeros(1,N),'like',A);
   *   F = cast(zeros(1,N-1),'like',A);
   *   for i=1:N
   *     l = i+1;
   *     % Left-hand transformation
   *     s = sum(A(i:M,i).^2);
   *     if s>=tol
   *       f = A(i,i); g = sqrt(s); if f>0, g = -g; end;
   *       h = f*g-s; A(i,i) = f-g;
   *       for j=l:N
   *         s = sum(A(i:M,i).*A(i:M,j))/h;
   *         A(i:M,j) = A(i:M,j)+s*A(i:M,i);
   *       end
   *     else
   *       g = 0;
   *     end
   *     D(i) = g;
   *     % Right-hand transformation
   *     if i<N-1
   *       s = sum(A(i,l:N).^2);
   *       if s>=tol
   *         f = A(i,l); g = sqrt(s); if f>0, g = -g; end;
   *         h = f*g-s; A(i,l) = f-g;
   *         for j=l:M
   *           s = sum(A(i,l:N).*A(j,l:N))/h;
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
   *   V = cast(rand(N),'like',A);
   *   V(N-1,N-1) = 1; V(N,N) = 1;
   *   for i=N-2:-1:1
   *     l = i+1; g = F(i);
   *     if g~=0
   *       h = A(i,l)*g;
   *       for j=l+1:N
   *         s = sum(V(l+1:N,j).*A(i,l+1:N)')/h;
   *         V(l,j) = s*A(i,l);
   *         V(l+1:N,j) = V(l+1:N,j)+s*A(i,l+1:N)';
   *       end
   *       V(l:N,l) = A(i,l:N)'/g;
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
   *       h = A(i,i)*g;
   *       for j=l:N
   *         s = sum(A(l:M,i).*A(l:M,j))/h;
   *         A(i,j) = s*A(i,i);
   *         A(l:M,j) = A(l:M,j)+s*A(l:M,i);
   *       end
   *       A(i:M,i) = A(i:M,i)/g;
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

        xb_vecN_2xf32 * restrict D_rw = (xb_vecN_2xf32*)D;
        xb_vecN_2xf32 * restrict F_rw = (xb_vecN_2xf32*)F;
        xb_vecN_2xf32 * restrict V_rw = (xb_vecN_2xf32*)V;
  const xb_vecN_2xf32 *          V_r;
        xb_vecN_2xf32 * restrict V_w0;
        xb_vecN_2xf32 * restrict V_w1;
        xb_vecN_2xf32 * restrict A_rw = (xb_vecN_2xf32*)A;
  const xb_vecN_2xf32 *          A_r0;
  const xb_vecN_2xf32 *          A_r1;
        xb_vecN_2xf32 * restrict A_w;

  const xb_vecN_2xf32 tol = FLT_MIN;
  const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
  const xb_vecN_2xf32 c1f = BBE_CONSTN_2XF32(1);

  xb_vecN_2xf32 f,g,h,r,s,t;
  vboolN_2 btol,bgtz,bnz;
  int i,j,k,l,p;

  #define SIDX(i,j,k)  (((i)*N+(j))*L+(k))

  NASSERT_ALIGN( D, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( F, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A, 2*BBE_SIMD_WIDTH );
  NASSERT( M>=N && N>1 );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  for ( k=0; k<L; k+=BBE_SIMD_WIDTH/2 ) {
    /*
     * Perform the bidiagonalization.
     */
#if 1
    for ( i=0; i<N; i++ ) {
      l = i+1;
      /*============= Left-hand transformation =============*/
      A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,i,k)];
      __Pragma("loop_count min=1");
      for ( s=c0f, j=i; j<M; j++ ) {
        BBE_LVN_2XF32_XP(f, A_r0, N*L*sz_f32);
        BBE_MULAN_2XF32(s,f,f);
      }
      /* Mask out columns whose L2 norm is finite, but so small that
       * it cannot be evaluated to working precision. */
      btol = BBE_ULEN_2XF32(tol,s);
      /* Construct a Householder vector */
      f = BBE_LVN_2XF32_X(A_rw, SIDX(i,i,k)*sz_f32);
      g = c0f; BBE_SQRTN_2XF32T(g,s,btol);
      bgtz = BBE_OLTN_2XF32(c0f,f);
      BBE_NEGN_2XF32T(g,g,bgtz);
      BBE_SVN_2XF32_X(g, D_rw, SIDX(0,i,k)*sz_f32);
      h = BBE_NEGN_2XF32(s);
      BBE_MULAN_2XF32(h,f,g);
      t = BBE_SUBN_2XF32(f,g);
      BBE_SVN_2XF32T_X(t, A_rw, SIDX(i,i,k)*sz_f32,btol);
      /* Apply the Householder reflector from the left to rows i..M-1,
       * columns i+1..N-1 */
      for ( j=l; j<N; j++ ) {
        A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,i,k)];
        A_r1 = (xb_vecN_2xf32*)&A[SIDX(i,j,k)];
        __Pragma("loop_count min=1");
        for ( s=c0f, p=i; p<M; p++ ) {
          BBE_LVN_2XF32_XP(r, A_r0, N*L*sz_f32);
          BBE_LVN_2XF32_XP(t, A_r1, N*L*sz_f32);
          BBE_MULAN_2XF32(s,r,t);
        }
        t = BBE_RECIP0N_2XF32(h);
        r = BBE_CONSTN_2XF32(1);
        BBE_MULSN_2XF32(r,t,h);
        BBE_MULAN_2XF32(t,r,t);
        f = BBE_MULN_2XF32(s,t);
        BBE_MULSN_2XF32(s,f,h);
        BBE_MULAN_2XF32(f,s,t);
        A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,i,k)];
        A_r1 = (xb_vecN_2xf32*)&A[SIDX(i,j,k)];
        A_w  = (xb_vecN_2xf32*)&A[SIDX(i,j,k)];
        __Pragma("loop_count min=1");
        for ( p=i; p<M; p++ ) {
          BBE_LVN_2XF32_XP(r, A_r0, N*L*sz_f32);
          BBE_LVN_2XF32_XP(t, A_r1, N*L*sz_f32);
          BBE_MULAN_2XF32(t,r,f);
          BBE_SVN_2XF32T_XP(t, A_w, N*L*sz_f32, btol);
        }
      } /* j */
      __Pragma("no_reorder");
      /*============= Right-hand transformation =============*/
      if (l<N-1) {
        A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l,k)];
        __Pragma("loop_count min=1");
        for ( s=c0f, j=l; j<N; j++ ) {
          BBE_LVN_2XF32_XP(f, A_r0, L*sz_f32);
          BBE_MULAN_2XF32(s,f,f);
        }
        /* Mask out rows whose L2 norm is finite, but so small that
         * it cannot be evaluated to working precision. */
        btol = BBE_ULEN_2XF32(tol,s);
        /* Construct a Householder vector. */
        f = BBE_LVN_2XF32_X(A_rw, SIDX(i,l,k)*sz_f32);
        g = c0f; BBE_SQRTN_2XF32T(g,s,btol);
        bgtz = BBE_OLTN_2XF32(c0f,f);
        BBE_NEGN_2XF32T(g,g,bgtz);
        BBE_SVN_2XF32_X(g, F_rw, SIDX(0,i,k)*sz_f32);
        h = BBE_NEGN_2XF32(s);
        BBE_MULAN_2XF32(h,f,g);
        t = BBE_SUBN_2XF32(f,g);
        BBE_SVN_2XF32T_X(t, A_rw, SIDX(i,l,k)*sz_f32,btol);
        /* Apply the Householder reflector from the right to rows i+1..M-1,
         * columns i+1..N-1 */
        for ( j=l; j<M; j++ ) {
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l,k)];
          A_r1 = (xb_vecN_2xf32*)&A[SIDX(j,l,k)];
          __Pragma("loop_count min=1");
          for ( s=c0f, p=l; p<N; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, L*sz_f32);
            BBE_LVN_2XF32_XP(t, A_r1, L*sz_f32);
            BBE_MULAN_2XF32(s,r,t);
          }
          t = BBE_RECIP0N_2XF32(h);
          r = BBE_CONSTN_2XF32(1);
          BBE_MULSN_2XF32(r,t,h);
          BBE_MULAN_2XF32(t,r,t);
          f = BBE_MULN_2XF32(s,t);
          BBE_MULSN_2XF32(s,f,h);
          BBE_MULAN_2XF32(f,s,t);
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l,k)];
          A_r1 = (xb_vecN_2xf32*)&A[SIDX(j,l,k)];
          A_w  = (xb_vecN_2xf32*)&A[SIDX(j,l,k)];
          __Pragma("loop_count min=1");
          for ( p=l; p<N; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, L*sz_f32);
            BBE_LVN_2XF32_XP(t, A_r1, L*sz_f32);
            BBE_MULAN_2XF32(t,r,f);
            BBE_SVN_2XF32T_XP(t, A_w, L*sz_f32, btol);
          }
        } /* j */
      } else if (l<N) {
        g = BBE_LVN_2XF32_X(A_rw, SIDX(N-2,N-1,k)*sz_f32);
        BBE_SVN_2XF32_X(g, F_rw, SIDX(0,i,k)*sz_f32);
      } /* l */
      __Pragma("no_reorder");
    } /* i */
#else
    {
      int m;
      for ( m=k; m<k+BBE_SIMD_WIDTH/2; m++ ) {
        const float32_t tol = FLT_MIN;
        float32_t f,g,h,s;
        for ( i=0; i<N; i++ ) {
          l = i+1;
          /* Left-hand transformation */
          for ( s=0.f, j=i; j<M; j++ ) { 
            f = A[SIDX(j,i,m)]; s += f*f; 
          }
          /* Use inverted condition for robust NaN propagation. */
          if (!(s<tol)) {
            f = A[SIDX(i,i,m)];
            g = sqrtf(s); if (f>0) g = -g;
            h = f*g-s; A[SIDX(i,i,m)] = f-g;
            for ( j=l; j<N; j++ ) {
              for ( s=0.f, p=i; p<M; p++ ) {
                s += A[SIDX(p,i,m)]*A[SIDX(p,j,m)];
              }
              s /= h;
              for ( p=i; p<M; p++ ) {
                A[SIDX(p,j,m)] += s*A[SIDX(p,i,m)];
              }
            }
          } else {
            g = 0.f;
          }
          D[SIDX(0,i,m)] = g;
          /* Right-hand transformation */
          if (l<N-1) {
            for ( s=0.f, j=l; j<N; j++ ) {
              f = A[SIDX(i,j,m)]; s += f*f;
            }
            if (!(s<tol)) {
              f = A[SIDX(i,l,m)];
              g = sqrtf(s); if (f>0) g = -g;
              h = f*g-s; A[SIDX(i,l,m)] = f-g;
              for ( j=l; j<M; j++ ) {
                for ( s=0.f, p=l; p<N; p++ ) {
                  s += A[SIDX(i,p,m)]*A[SIDX(j,p,m)];
                }
                s /= h;
                for ( p=l; p<N; p++ ) {
                  A[SIDX(j,p,m)] += s*A[SIDX(i,p,m)];
                }
              }
            } else {
              g = 0.f;
            }
            F[SIDX(0,i,m)] = g;
          } else if (l<N) {
            F[SIDX(0,i,m)] = A[SIDX(N-2,N-1,m)];
          } /* l */
        } /* i */
      } /* m */
    }
#endif
    /*
     * If needed, accumulate right-hand transformations from reflection axes 
     * stored in  the work matrix A.
     */
#if 1
    if (V) {
      /* V[IDX_K(N-2,N-2,k)] = V[IDX_K(N-1,N-1,k)] = 1.f; */
      BBE_SVN_2XF32_X(c1f, V_rw, SIDX(N-2,N-2,k)*sz_f32);
      BBE_SVN_2XF32_X(c1f, V_rw, SIDX(N-1,N-1,k)*sz_f32);
      __Pragma("no_reorder");
      for ( i=N-3; i>=0; i-- ) {
        /* l = i+1; g = F[IDX_K(0,i,k)]; */
        l = i+1; g = BBE_LVN_2XF32_X(F_rw, SIDX(0,i,k)*sz_f32);
        bnz = BBE_NOTBN_2(BBE_UEQN_2XF32(c0f,g));
        /* h = A[IDX_K(i,l,k)]*g; */
        t = BBE_LVN_2XF32_X(A_rw, SIDX(i,l,k)*sz_f32);
        h = BBE_MULN_2XF32(t,g);
        for ( j=l+1; j<N; j++ ) {
          /* for ( s=0.f, p=l+1; p<N; p++ ) {
           *   s += V[IDX_K(p,j,k)]*A[IDX_K(i,p,k)];
           * } */
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l+1,k)];
          V_r  = (xb_vecN_2xf32*)&V[SIDX(l+1,j,k)];
          for ( s=c0f, p=l+1; p<N; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, L*sz_f32);
            BBE_LVN_2XF32_XP(t, V_r, N*L*sz_f32);
            BBE_MULAN_2XF32(s,r,t);
          }
          /* f = s/h; */
          t = BBE_RECIP0N_2XF32(h);
          r = BBE_CONSTN_2XF32(1);
          BBE_MULSN_2XF32(r,t,h);
          BBE_MULAN_2XF32(t,r,t);
          f = BBE_MULN_2XF32(s,t);
          BBE_MULSN_2XF32(s,f,h);
          BBE_MULAN_2XF32(f,s,t);
          /* V[IDX_K(l,j,k)] = f*A[IDX_K(i,l,k)]; */
          r = BBE_LVN_2XF32_X(A_rw, SIDX(i,l,k)*sz_f32);
          t = c0f; BBE_MULN_2XF32T(t,r,f,bnz);
          BBE_SVN_2XF32_X(t, V_rw, SIDX(l,j,k)*sz_f32);
          /* for ( p=l+1; p<N; p++ ) {
           *   V[IDX_K(p,j,k)] += s*A[IDX_K(i,p,k)];
           * } */
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l+1,k)];
          V_r  = (xb_vecN_2xf32*)&V[SIDX(l+1,j,k)];
          V_w0 = (xb_vecN_2xf32*)&V[SIDX(l+1,j,k)];
          for ( p=l+1; p<N; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, L*sz_f32);
            BBE_LVN_2XF32_XP(t, V_r, N*L*sz_f32);
            BBE_MULAN_2XF32T(t,r,f,bnz);
            BBE_SVN_2XF32_XP(t, V_w0, N*L*sz_f32);
          }
        } /* j */
        /* for ( p=l; p<N; p++ ) {
         *   V[IDX_K(p,l,k)] = A[IDX_K(i,p,k)]/g;
         * } */
        h = c0f; BBE_RECIP0N_2XF32T(h,g,bnz);
        r = BBE_CONSTN_2XF32(1);
        BBE_MULSN_2XF32(r,g,h);
        BBE_MULAN_2XF32T(h,r,h,bnz);
        A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,l,k)];
        V_w0 = (xb_vecN_2xf32*)&V[SIDX(l,l,k)];
        for ( p=l; p<N; p++ ) {
          BBE_LVN_2XF32_XP(r, A_r0, L*sz_f32);
          t = BBE_MULN_2XF32(r,h);
          BBE_MULSN_2XF32(r,t,g);
          BBE_MULAN_2XF32(t,r,h);
          BBE_SVN_2XF32_XP(t, V_w0, N*L*sz_f32);
        }
        __Pragma("no_reorder");
        /* V[IDX_K(l,l,k)] += 1.f; */
        r = BBE_LVN_2XF32_X(V_rw, SIDX(l,l,k)*sz_f32);
        t = BBE_ADDN_2XF32(r,c1f);
        BBE_SVN_2XF32_X(t, V_rw, SIDX(l,l,k)*sz_f32);
        __Pragma("no_reorder");
      } /* i */
      /* V[IDX_K(0,0,k)] = 1.f;
       * for ( i=1; i<N; i++ ) {
       *   V[IDX_K(i,0,k)] = V[IDX_K(0,i,k)] = 0.f;
       * } */
      BBE_SVN_2XF32_X(c1f, V_rw, SIDX(0,0,k)*sz_f32);
      V_w0 = (xb_vecN_2xf32*)&V[SIDX(1,0,k)];
      V_w1 = (xb_vecN_2xf32*)&V[SIDX(0,1,k)];
      for ( i=1; i<N; i++ ) {
        BBE_SVN_2XF32_XP(c0f, V_w0, N*L*sz_f32);
        BBE_SVN_2XF32_XP(c0f, V_w1, L*sz_f32);
      }
    } /* V */
#else
    if (V) {
      float32_t g,h,s;
      int m;
      for ( m=k; m<k+BBE_SIMD_WIDTH/2; m++ ) {
        V[SIDX(N-2,N-2,m)] = V[SIDX(N-1,N-1,m)] = 1.f;
        for ( i=N-3; i>=0; i-- ) {
          l = i+1; g = F[SIDX(0,i,m)];
          if (g!=0) {
            h = A[SIDX(i,l,m)]*g;
            for ( j=l+1; j<N; j++ ) {
              for ( s=0.f, p=l+1; p<N; p++ ) {
                s += V[SIDX(p,j,m)]*A[SIDX(i,p,m)];
              }
              s /= h;
              V[SIDX(l,j,m)] = s*A[SIDX(i,l,m)];
              for ( p=l+1; p<N; p++ ) {
                V[SIDX(p,j,m)] += s*A[SIDX(i,p,m)];
              }
            }
            for ( p=l; p<N; p++ ) {
              V[SIDX(p,l,m)] = A[SIDX(i,p,m)]/g;
            }
          } else {
            V[SIDX(l,l,m)] = 0.f;
            for ( j=l+1; j<N; j++ ) {
              V[SIDX(j,l,m)] = V[SIDX(l,j,m)] = 0.f;
            }
          } /* g */
          V[SIDX(l,l,m)] += 1.f;
        } /* i */
        V[SIDX(0,0,m)] = 1.f;
        for ( i=1; i<N; i++ ) {
          V[SIDX(i,0,m)] = V[SIDX(0,i,m)] = 0.f;
        }
      } /* m */
    } /* V */
#endif
    /*
     * If needed, accumulate left-hand transformations. Matrix U overrides 
     * intermediate data stored in A. 
     */
#if 1
    if (needU) {
      for ( i=N-1; i>=0; i-- ) {
        /* l = i+1; g = D[IDX_K(0,i,k)]; */
        l = i+1; g = BBE_LVN_2XF32_X(D_rw, SIDX(0,i,k)*sz_f32);
        bnz = BBE_NOTBN_2(BBE_UEQN_2XF32(c0f,g));
        /* h = A[IDX_K(i,i,k)]*g; */
        t = BBE_LVN_2XF32_X(A_rw, SIDX(i,i,k)*sz_f32);
        h = BBE_MULN_2XF32(t,g);
        for ( j=l; j<N; j++ ) {
          /* for ( s=0.f, p=l; p<M; p++ ) {
           *   s += A[IDX_K(p,i,k)]*A[IDX_K(p,j,k)];
           * } */
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(l,i,k)];
          A_r1 = (xb_vecN_2xf32*)&A[SIDX(l,j,k)];
          for ( s=c0f, p=l; p<M; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, N*L*sz_f32);
            BBE_LVN_2XF32_XP(t, A_r1, N*L*sz_f32);
            BBE_MULAN_2XF32(s,r,t);
          }
          /* f = s/h; */
          t = BBE_RECIP0N_2XF32(h);
          r = BBE_CONSTN_2XF32(1);
          BBE_MULSN_2XF32(r,t,h);
          BBE_MULAN_2XF32(t,r,t);
          f = BBE_MULN_2XF32(s,t);
          BBE_MULSN_2XF32(s,f,h);
          BBE_MULAN_2XF32(f,s,t);
          /* A[IDX_K(i,j,k)] = s*A[IDX_K(i,i,k)]; */
          r = BBE_LVN_2XF32_X(A_rw, SIDX(i,i,k)*sz_f32);
          t = c0f; BBE_MULN_2XF32T(t,r,f,bnz);
          BBE_SVN_2XF32_X(t, A_rw, SIDX(i,j,k)*sz_f32);
          /* for ( p=l; p<M; p++ ) {
           *   A[IDX_K(p,j,k)] += s*A[IDX_K(p,i,k)];
           * } */
          A_r0 = (xb_vecN_2xf32*)&A[SIDX(l,i,k)];
          A_r1 = (xb_vecN_2xf32*)&A[SIDX(l,j,k)];
          A_w  = (xb_vecN_2xf32*)&A[SIDX(l,j,k)];
          for ( p=l; p<M; p++ ) {
            BBE_LVN_2XF32_XP(r, A_r0, N*L*sz_f32);
            BBE_LVN_2XF32_XP(t, A_r1, N*L*sz_f32);
            BBE_MULAN_2XF32T(t,r,f,bnz);
            BBE_SVN_2XF32_XP(t, A_w, N*L*sz_f32);
          }
        } /* j */
        /* for ( j=i; j<M; j++ ) {
         *   A[IDX_K(j,i,k)] /= g;
         * } */
        h = c0f; BBE_RECIP0N_2XF32T(h,g,bnz);
        r = BBE_CONSTN_2XF32(1);
        BBE_MULSN_2XF32(r,g,h);
        BBE_MULAN_2XF32T(h,r,h,bnz);
        A_r0 = (xb_vecN_2xf32*)&A[SIDX(i,i,k)];
        A_w  = (xb_vecN_2xf32*)&A[SIDX(i,i,k)];
        for ( j=i; j<M; j++ ) {
          BBE_LVN_2XF32_XP(r, A_r0, N*L*sz_f32);
          t = BBE_MULN_2XF32(r,h);
          BBE_MULSN_2XF32(r,t,g);
          BBE_MULAN_2XF32(t,r,h);
          BBE_SVN_2XF32_XP(t, A_w, N*L*sz_f32);
        }
        __Pragma("no_reorder");
        /* A[IDX_K(i,i,k)] += 1.f; */
        r = BBE_LVN_2XF32_X(A_rw, SIDX(i,i,k)*sz_f32);
        t = BBE_ADDN_2XF32(r,c1f);
        BBE_SVN_2XF32_X(t, A_rw, SIDX(i,i,k)*sz_f32);
        __Pragma("no_reorder");
      } /* i */
    } /* needU */
#else
    if (needU) {
      float32_t g,h,s;
      int i,j,l,m,p;
      for ( m=k; m<k+BBE_SIMD_WIDTH/2; m++ ) {
        for ( i=N-1; i>=0; i-- ) {
          l = i+1; g = D[SIDX(0,i,m)];
          if (g!=0) {
            h = A[SIDX(i,i,m)]*g;
            for ( j=l; j<N; j++ ) {
              for ( s=0.f, p=l; p<M; p++ ) {
                s += A[SIDX(p,i,m)]*A[SIDX(p,j,m)];
              }
              s /= h;
              A[SIDX(i,j,m)] = s*A[SIDX(i,i,m)];
              for ( p=l; p<M; p++ ) {
                A[SIDX(p,j,m)] += s*A[SIDX(p,i,m)];
              }
            }
            for ( j=i; j<M; j++ ) {
              A[SIDX(j,i,m)] /= g;
            }
          } else {
            for ( j=i; j<M; j++ ) {
              A[SIDX(j,i,m)] = 0.f;
            }
            for ( j=l; j<N; j++ ) {
              A[SIDX(i,j,m)] = 0.f;
            }
          }
          A[SIDX(i,i,m)] += 1.f;
        } /* i */
      } /* m */
    } /* needU */
#endif
  } /* k */

  #undef SIDX

} /* rsvd_bidiag_mxnsf() */

#endif /* HAVE_VFPU */

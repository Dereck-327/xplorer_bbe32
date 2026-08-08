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
    Golub-Reinsch SVD algorithm for real upper bidiagonal matrices stored in
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
/* NaN values for single precision routines */
#include "nanf_tbl.h"
/* SVD common declarations */
#include "svd_common.h"

#if HAVE_VFPU

#define EPS          FLT_EPSILON
#define ITS_LIM      75 /* Iterations count limit */

#define MAX(a,b)     ((a)>(b)?(a):(b))
#define SIDX(i,k)    ((i)*L+(k))

#define sz_i16       sizeof(int16_t)
#define sz_f32       sizeof(float32_t)
#define sz_vbn2      sizeof(vboolN_2)

/* Internal function of Golub-Kahan SVD step implementation.
 * Summarize index data and fetch data for Wilkinson's shift computation. */
typedef void gks_fetch_fxn_t( 
                    int16_t * a_its, int16_t * rng_l, int16_t * rng_k, vboolN_2  * a_m,
                    float32_t * a_a, float32_t * a_b, float32_t * a_x, float32_t * a_y,
                    float32_t * a_c, float32_t * a_s,
              const int16_t   * a_l, const int16_t   * a_k,
              const float32_t * D  , const float32_t * F,
              int itsLim, int N, int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Compute Wilkinson's shift. */
typedef void gks_wilkShift_fxn_t(
                    float32_t * a_a,       float32_t * a_b,
                    float32_t * a_c,       float32_t * a_s,
              const float32_t * a_x, const float32_t * a_y,
              int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Compute real Givens's rotator matrix: 
 *   G(a,b) <- [c,s;-s,c]: [a,b]*[c,s;-s,c] == [*,0], c^2+s^2 == 1 */
typedef void gks_givens_fxn_t( 
                    float32_t * a_c,        float32_t * a_s,
              const float32_t * a_a,  const float32_t * a_b, 
              int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * 1st updating step of SVD step iteration: B <- B*G(a,b). */
typedef void gks_step1_fxn_t(
                    float32_t * D  ,       float32_t * F,
                    float32_t * a_a,       float32_t * a_b,
              const float32_t * a_c, const float32_t * a_s,
              const int16_t   * a_l, const vboolN_2  * a_m,
              int n, int L );

/*Internal function of Golub-Kahan SVD step implementation.
 * 2nd updating step of SVD step iteration: B <- G(a,b)'*B. */
typedef void gks_step2_fxn_t(
                    float32_t * D,         float32_t * F,
                    float32_t * a_a,       float32_t * a_b,
              const float32_t * a_c, const float32_t * a_s,
              const int16_t   * a_k, const vboolN_2  * a_m,
              int n, int L );

/* Internal function of Golub-Kahan SVD step implementation.
 * Accumulate transformation: U <- U*G(a,b) or V <- V*G(a,b). */
typedef void gks_accum_fxn_t(
                    float32_t * W,   const float32_t * a_c,
              const float32_t * a_s, const vboolN_2  * a_m,
              int n, int M, int N, int L );

/* Golub-Reinsch SVD scratch area */
typedef struct tag_rgrsvdsf_scratch
{
  float32_t *tol,*a,*b,*c,*s,*x,*y;
  int16_t *its,*k,*l,*n;
  vboolN_2 *m;

} rgrsvdsf_scratch_t;

/* Golub-Kahan SVD Step internal functions */
typedef struct tag_rgrsvdsf_gksApi
{
  gks_fetch_fxn_t     * fetch;
  gks_wilkShift_fxn_t * wilkShift;
  gks_givens_fxn_t    * givens;
  gks_step1_fxn_t     * step1;
  gks_step2_fxn_t     * step2;
  gks_accum_fxn_t     * accum;

} rgrsvdsf_gksApi_t;

/* Golub-Kahan SVD Step internal functions set for generic matrix sizes. */
static const rgrsvdsf_gksApi_t rgrsvdsf_gksApi_mxn = {
  rgrsvdsf_gks_fetch, rgrsvdsf_gks_wilkShift, rgrsvdsf_gks_givens,
  rgrsvdsf_gks_step1, rgrsvdsf_gks_step2    , rgrsvdsf_gks_accum
};

/* Golub-Reinsch SVD for real upper bidiagonal matrices. Function processes L 
 * matrices stored in stream order.
 * Based on [1] Algorithm 8.6.2. */
static void golubReinschSVD( 
                 const rgrsvdsf_scratch_t * scr,
                 const rgrsvdsf_gksApi_t  * gksApi,
                       float32_t * restrict D, /* [N][L]   (in/out) */
                       float32_t * restrict F, /* [N-1][L] (in/tmp) */
                       float32_t * restrict U, /* [M*N][L] (in/out) */
                       float32_t * restrict V, /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Zero conditionally accessed scratch arrays. */
static void wipeScratch(
                 const rgrsvdsf_scratch_t * scr,
                 int L );

/* Derive tolerance level from maximum absolute value over a matrix. */
static void setTolLvl( float32_t * restrict a_tol, /* [L]      (out) */
                 const float32_t * restrict D,     /* [N][L]   (in)  */
                 const float32_t * restrict F,     /* [N-1][L] (in(  */
                 int N, int L );

/* Look for a diagonal submatrix at the bottom and deflate the SVD. */
static void tryDeflate(
                       int16_t   * restrict a_k,   /* [L]      (in/out) */
                       int16_t   * restrict a_its, /* [L]      (in/out) */
                       float32_t * restrict D,     /* [N][L]   (in/out) */
                       float32_t * restrict F,     /* [N-1][L] (in/out) */
                       float32_t * restrict V,     /* [N*N][L] (in/out) */
                 int N, int L );

/* Search for a zero on the superdiagonal or on the main diagonal. If found, 
 * split the matrix. */
static void trySplit(  int16_t   * restrict a_n,   /* [L]      (tmp)    */
                       int16_t   * restrict a_l,   /* [L]      (out)    */
                 const int16_t   * restrict a_k,   /* [L]      (in)     */
                 const float32_t * restrict a_tol, /* [L]      (in)     */
                       float32_t * restrict D,     /* [N][L]   (in/out) */
                       float32_t * restrict F,     /* [N-1][L] (in/out) */
                       float32_t * restrict U,     /* [M*N][L] (in/out) */
                 int k_up, int M, int N, int L );

/* Golub-Kahan SVD step for principal submatrices in rows/cols a_l[p]..a_k[p].
 * p=0..L-1. Function processes L matrices stored in stream order. 
 * Returns the highest index k among all matrices with a_its[]<ITS_LIM.
 * Based on [1] Algorithm 8.6.1. */
static int16_t golubKahanSVDStep( 
                 const rgrsvdsf_scratch_t * scr,
                 const rgrsvdsf_gksApi_t  * gksApi,
                       float32_t * restrict D, /* [N][L]   (in/out) */
                       float32_t * restrict F, /* [N-1][L] (in/out) */
                       float32_t * restrict U, /* [M*N][L] (in/out) */
                       float32_t * restrict V, /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Sort singular values in descending order, and permute left- and right-singular
 * vectors (columns of matrices U and V) accordingly. Data are reordered in-place. */
static void sortSingVal(
                       float32_t * restrict D,   /* [N][L]   (in/out) */
                       float32_t * restrict U,   /* [M*N][L] (in/out) */
                       float32_t * restrict V,   /* [N*N][L] (in/out) */
                 int M, int N, int L );

/* Compute Givens's rotation matrix.
 * Given scalars a and b, compute scalars s and c such that (in MATLAB notation):
 *   A) [a,b]*[c,s;-s,c] == [*,0] and 
 *   B) c^2+s^2 == 1 */
static void givens( float32_t * restrict c, 
                    float32_t * restrict s, 
                    float32_t a, float32_t b );

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

void rgrsvdsf( void * pScr,
               float32_t * restrict D,
               float32_t * restrict F,
               float32_t * restrict U,
               float32_t * restrict V,
               int M, int N, int L )
{
  const rgrsvdsf_gksApi_t * gksApi = &rgrsvdsf_gksApi_mxn;
  rgrsvdsf_scratch_t scr;

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
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
    int LW2 = L/(BBE_SIMD_WIDTH/2);
    void * p = pScr;
    memset(&scr, 0, sizeof(scr));
    scr.a   = (float32_t*)p; p = scr.a + L;
    scr.b   = (float32_t*)p; p = scr.b + L;
    scr.x   = (float32_t*)p; p = scr.x + L;
    scr.y   = (float32_t*)p; p = scr.y + L;
    scr.c   = (float32_t*)p; p = scr.c + L;
    scr.s   = (float32_t*)p; p = scr.s + L;
    scr.tol = (float32_t*)p; p = scr.tol + L;
    scr.its = (int16_t  *)p; p = scr.its + Lix;
    scr.k   = (int16_t  *)p; p = scr.k + Lix;
    scr.l   = (int16_t  *)p; p = scr.l + Lix;
    scr.n   = (int16_t  *)p; p = scr.n + Lix;
    scr.m   = (vboolN_2 *)p; p = scr.m + (N-1)*LW2;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)rgrsvdsf_getScratchSize(M,N,L) );
  }

  golubReinschSVD(&scr,gksApi,D,F,U,V,M,N,L);

} /* rgrsvdsf() */

size_t rgrsvdsf_getScratchSize( int M, int N, int L )
{
  int Lix = (L+BBE_SIMD_WIDTH-1)/BBE_SIMD_WIDTH*BBE_SIMD_WIDTH;
  int LW2 = L/(BBE_SIMD_WIDTH/2);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));
  NASSERT(M>=N);
  NASSERT(N>1);
  return ( 7*L  *sz_f32 +   /* Floating-point data */
           4*Lix*sz_i16 +   /* Indexing data       */
       (N-1)*LW2*sz_vbn2 ); /* Vector boolean data */
}

/* Golub-Reinsch SVD for real upper bidiagonal matrices. Function processes L 
 * matrices stored in stream order.
 * Based on [1] Algorithm 8.6.2. */
void golubReinschSVD( const rgrsvdsf_scratch_t * scr,
                      const rgrsvdsf_gksApi_t  * gksApi,
                      float32_t * restrict D, /* [N][L]   (in/out) */
                      float32_t * restrict F, /* [N-1][L] (in/tmp) */
                      float32_t * restrict U, /* [M*N][L] (in/out) */
                      float32_t * restrict V, /* [N*N][L] (in/out) */
                      int M, int N, int L )
{
  /*
   * MATLAB reference code:
   *
   *   function varargout = rsvd(A)
   *   assert(isreal(A),'Input argument A must be real');
   *   nout = max(1,nargout);
   *   needU = nout>1; needV = nout>2;
   *   n = size(A,2);
   *   EPS = eps(cast(1,'like',A));
   *   % Convert input matrix A to bidiagonal form: A = U*B*V'.
   *   if nout>1
   *     if needV, [U,B,V] = rbidiag(A); else [U,B] = rbidiag(A); end;
   *     d = B(1:n+1:n*n); f = B(n+1:n+1:n*n);
   *   else
   *     b = rbidiag(A);
   *     d = b(1,:); f = b(2,2:n);
   *   end
   *   % Preserve the data type for empty matrices because it impacts other
   *   % arguments when passed to a MATLAB function.
   *   if ~needU, U = cast([],'like',A); end;
   *   if ~needV, V = cast([],'like',A); end;
   *   tol = EPS*max(abs([d,f]));
   *   iterCnt = zeros(1,n);
   *   k = n;
   *   % This loop implements the Golub-Reinsch SVD, [1] Algorithm 8.6.2.
   *   while k>1
   *     % Look for a diagonal submatrix at the bottom.
   *     while k==1 || (k>1 && abs(f(k-1))<=EPS*(abs(d(k-1))+abs(d(k))))
   *       % Tiny elements on the superdiagonal are flushed to zero.
   *       if k>1, f(k-1) = 0; end;
   *       % If necessary, negate the result and adjust the matrix V accordingly.
   *       if d(k)<0
   *         d(k) = -d(k); 
   *         if ~isempty(V), V(:,k) = -V(:,k); end;
   *       end
   *       % Deflate the problem
   *       k = k-1;
   *     end
   *     % Search for a zero on the superdiagonal. If found, split the matrix.
   *     l = k;
   *     while l>1
   *       if abs(f(l-1))<=EPS*(abs(d(l-1))+abs(d(l)))
   *         % Set a tiny superdiagonal element to zero and decouple the problem
   *         f(l-1) = 0; break;
   *       elseif abs(d(l-1))<=tol
   *         % If d(l-1) is zero, set f(l-1) to zero through plane rotations
   *         % (l-1,l),(l-1,l+1),...,(l-1,k).
   *         c = 1; s = 0;
   *         for p=l:k
   *           a = f(p-1); f(p-1) = -s*a;
   *           a = c*a; b = d(p);
   *           % G(a,b) <- [c,s;-s,c]: [c,s;-s,c]'*[a;b] == [*;0]
   *           [c,s] = givens(a,b);
   *           d(p) = c*a-s*b;
   *           if ~isempty(U)
   *             u = U(:,l-1); v = U(:,p);
   *             U(:,l-1) = s*u+c*v;
   *             U(:,p) = c*u-s*v;
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
   *   % Sort singular values in descending order. We use a kind of the bubble 
   *   % sort due to simplicity. Its miserable running time of O(n^2) is 
   *   % completely buried under the overall complexity of O(n^3).
   *   if ~isempty(d)
   *     for k=1:n-1
   *       l = k;
   *       for p=k+1:n, 
   *         if d(p)>d(l), l = p; end;
   *       end;
   *       if k~=l
   *         t = d(l); d(l) = d(k); d(k) = t;
   *         if ~isempty(U), u = U(:,l); U(:,l) = U(:,k); U(:,k) = u; end;
   *         if ~isempty(V), v = V(:,l); V(:,l) = V(:,k); V(:,k) = v; end;
   *       end
   *     end
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
    tryDeflate(a_k,a_its,D,F,V,N,L);
    /* Search for a zero on the superdiagonal or on the main diagonal. 
     * If found, split the matrix. */
    trySplit(a_n,a_l,a_k,a_tol,D,F,U,k_up,M,N,L);
#if 0
    for ( p=0; p<L; p++ ) {
      fprintf(f_a_l, "%d ", (int)scr->l[p]);
      fprintf(f_a_k, "%d ", (int)scr->k[p]);
      fprintf(f_a_its, "%d ", (int)scr->its[p]);
    }
    fprintf(f_a_l, "\n");
    fprintf(f_a_k, "\n");
    fprintf(f_a_its, "\n");
#endif
    /* Perform a single Golub-Kahan SVD step for each matrix */
    k_up = golubKahanSVDStep(scr,gksApi,D,F,U,V,M,N,L);
  } /* k_up */
  /* Check if the algorithm failed to converge for any input matrix. */
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>0) {
      NASSERT(a_its[p]>=ITS_LIM);
      for ( i=0; i<N; i++ ) D[SIDX(i,p)] = qNaNf.f;
    }
  }
  /* Sort singular values in descending order. */
  sortSingVal(D,U,V,M,N,L);

} /* golubReinschSVD() */

/* Zero conditionally accessed scratch arrays. */
void wipeScratch(
          const rgrsvdsf_scratch_t * scr,
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
  for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
    BBE_SVN_2XF32_IP(c0f, A_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, B_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, C_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, S_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, X_w, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(c0f, Y_w, 2*BBE_SIMD_WIDTH);
  }
#else
  int p;
  for ( p=0; p<L; p++ ) {
    scr->a[p] = scr->b[p] = 
    scr->c[p] = scr->s[p] =
    scr->x[p] = scr->y[p] = 0.f;
  }
#endif
} /* wipeScratch() */

/* Derive tolerance level from maximum absolute value over a matrix. */
void setTolLvl( float32_t * restrict a_tol, /* [L]      (out) */
          const float32_t * restrict D,     /* [N][L]   (in)  */
          const float32_t * restrict F,     /* [N-1][L] (in(  */
          int N, int L )
{
  /*
   * MATLAB outline:
   *   tol = EPS*max(abs([d,f]));
   */
#if 1
  const xb_vecN_2xf32 * restrict D_r;
  const xb_vecN_2xf32 * restrict F_r;
  const xb_vecN_2xf32 * restrict TOL_r;
        xb_vecN_2xf32 * restrict TOL_w;

  int i,p;

  D_r = (xb_vecN_2xf32*)&D[SIDX(N-1,0)];
  TOL_w = (xb_vecN_2xf32*)a_tol;
  for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
    xb_vecN_2xf32 a;
    BBE_LVN_2XF32_IP(a, D_r, 2*BBE_SIMD_WIDTH);
    a = BBE_MULN_2XF32(EPS, BBE_ABSN_2XF32(a));
    BBE_SVN_2XF32_IP(a, TOL_w, 2*BBE_SIMD_WIDTH);
  } /* p */
  __Pragma("no_reorder");
  D_r = (xb_vecN_2xf32*)&D[SIDX(0,0)];
  F_r = (xb_vecN_2xf32*)&F[SIDX(0,0)];
  for ( i=0; i<N-1; i++ ) {
    TOL_r = TOL_w = (xb_vecN_2xf32*)a_tol;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      xb_vecN_2xf32 a,b,c;
      BBE_LVN_2XF32_IP(a, D_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(b, F_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(c, TOL_r, 2*BBE_SIMD_WIDTH);
      a = BBE_MAXN_2XF32(BBE_ABSN_2XF32(a), BBE_ABSN_2XF32(b));
      c = BBE_MAXN_2XF32(c, BBE_MULN_2XF32(EPS, a));
      BBE_SVN_2XF32_IP(c, TOL_w, 2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");
  } /* i */
#else
  int i,p;
  for ( p=0; p<L; p++ ) {
    a_tol[p] = fabsf(D[SIDX(N-1,p)])*EPS;
  }
  for ( i=N-2; i>=0; i-- ) {
    for ( p=0; p<L; p++ ) {
      a_tol[p] = MAX(a_tol[p], MAX(fabsf(D[SIDX(i,p)]), fabsf(F[SIDX(i,p)]))*EPS);
    }
  } /* i */
#endif
#if 0
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
                int16_t   * restrict a_k,   /* [L]      (in/out) */
                int16_t   * restrict a_its, /* [L]      (in/out) */
                float32_t * restrict D,     /* [N][L]   (in/out) */
                float32_t * restrict F,     /* [N-1][L] (in/out) */
                float32_t * restrict V,     /* [N*N][L] (in/out) */
          int N, int L )
{
  /*
   * MATLAB outline:
   *   % Look for a diagonal submatrix at the bottom.
   *   while k==1 || (k>1 && abs(f(k-1))<=EPS*(abs(d(k-1))+abs(d(k))))
   *     % Tiny elements on the superdiagonal are flushed to zero.
   *     if k>1, f(k-1) = 0; end;
   *     % If necessary, negate the result and adjust the matrix V accordingly.
   *     if d(k)<0
   *       d(k) = -d(k); 
   *       if ~isempty(V), V(:,k) = -V(:,k); end;
   *     end
   *     % Deflate the problem
   *     k = k-1;
   *   end
   */
#if 1
  const float32_t     * restrict V_r;
        float32_t     * restrict V_w;

  float32_t a,b,c;
  int i,p;

  for ( p=0; p<L; p++ ) {
    while (a_k[p]>=0) {
      if (a_k[p]>0) {
        a = XT_ABS_S(F[SIDX(a_k[p]-1,p)]); b = XT_ABS_S(D[SIDX(a_k[p]-1,p)]); c = XT_ABS_S(D[SIDX(a_k[p],p)]);
        if (a>EPS*(b+c)) break;
        /* Tiny elements on the superdiagonal are flushed to zero. */
        F[SIDX(a_k[p]-1,p)] = 0;
      }
      /* If necessary, negate the result and adjust the matrix V accordingly. */
      if (D[SIDX(a_k[p],p)]<0) {
        D[SIDX(a_k[p],p)] = -D[SIDX(a_k[p],p)]; 
        if (V) {
          V_r = V_w = &V[SIDX(a_k[p],p)];
          __Pragma("no_unroll");
          for ( i=0; i<N; i++ ) {
            float32_t a;
            xtfloat_loadxp(a, V_r, N*L*sz_f32);
            xtfloat_storexp(-a, V_w, N*L*sz_f32);
          }
        }
      }
      /* Deflate the problem */
      a_its[p] = 0; a_k[p]--;
    } /* a_k[p] */
  } /* p */
#else
  float32_t a,b,c;
  int i,p;

  for ( p=0; p<L; p++ ) {
    while (a_k[p]>=0) {
      if (a_k[p]>0) {
        a = fabsf(F[SIDX(a_k[p]-1,p)]); b = fabsf(D[SIDX(a_k[p]-1,p)]); c = fabsf(D[SIDX(a_k[p],p)]);
        if (a>EPS*(b+c)) break;
        /* Tiny elements on the superdiagonal are flushed to zero. */
        F[SIDX(a_k[p]-1,p)] = 0;
      }
      /* If necessary, negate the result and adjust the matrix V accordingly. */
      if (D[SIDX(a_k[p],p)]<0) {
        D[SIDX(a_k[p],p)] = -D[SIDX(a_k[p],p)]; 
        if (V) for ( i=0; i<N; i++ ) V[SIDX(i*N+a_k[p],p)] = -V[SIDX(i*N+a_k[p],p)];
      }
      /* Deflate the problem */
      a_its[p] = 0; a_k[p]--;
    } /* a_k[p] */
  } /* p */
#endif
} /* tryDeflate() */

/* Search for a zero on the superdiagonal or on the main diagonal. If found, 
 * split the matrix. */
void trySplit(  int16_t   * restrict a_n,   /* [L]      (tmp)    */
                int16_t   * restrict a_l,   /* [L]      (out)    */
          const int16_t   * restrict a_k,   /* [L]      (in)     */
          const float32_t * restrict a_tol, /* [L]      (in)     */
                float32_t * restrict D,     /* [N][L]   (in/out) */
                float32_t * restrict F,     /* [N-1][L] (in/out) */
                float32_t * restrict U,     /* [M*N][L] (in/out) */
          int k_up, int M, int N, int L )
{
  /*
   * MATLAB outline:
   *   l = k;
   *   while l>1
   *     if abs(f(l-1))<=EPS*(abs(d(l-1))+abs(d(l)))
   *       % Set a tiny superdiagonal element to zero and decouple the problem
   *       f(l-1) = 0; break;
   *     elseif abs(d(l-1))<=tol
   *       % If d(l-1) is zero, set f(l-1) to zero through plane rotations
   *       % (l-1,l),(l-1,l+1),...,(l-1,k).
   *       c = 1; s = 0;
   *       for p=l:k
   *         a = f(p-1); f(p-1) = -s*a;
   *         a = c*a; b = d(p);
   *         % G(a,b) <- [c,s;-s,c]: [c,s;-s,c]'*[a;b] == [*;0]
   *         [c,s] = givens(a,b);
   *         d(p) = c*a-s*b;
   *         if ~isempty(U)
   *           u = U(:,l-1); v = U(:,p);
   *           U(:,l-1) = s*u+c*v;
   *           U(:,p) = c*u-s*v;
   *         end
   *       end
   *       % Restart the search.
   *       l = k; continue;
   *     end
   *     l = l-1;
   *   end
   */
#if 1
  const xb_vecN_2xf32 * restrict D_r;
  const xb_vecN_2xf32 * restrict F_r;
        xb_vecN_2xf32 * restrict F_w;
  const xb_vecN_2xf32 * restrict TOL_r;
  const float32_t     * restrict U_r;
        float32_t     * restrict U_w;
  const xb_vecNx16    * restrict K_r;
        xb_vecNx16    * restrict L_w;
        xb_vecNx16    * restrict N_w;

  float32_t a,b,c,s,u,v;
  int i,j,p,n_lo;

  for ( p=0; p<L; p++ ) {
    a_l[p] = a_n[p] = 0;
  } /* p */
  n_lo = 0;
  do {
#if 1
    D_r = (xb_vecN_2xf32*)&D[SIDX(n_lo,0)];
    F_r = F_w = (xb_vecN_2xf32*)&F[SIDX(n_lo,0)];
    for ( i=n_lo+1; i<=k_up; i++ ) {
      xb_vecN_2xf32 a0,a1,b0,b1,c0,c1,t0,t1;
      xb_vecNx16 vi,vk;
      vboolN_2 g0,g1,h0,h1;
      vboolN f,g,h;
      K_r = (xb_vecNx16*)a_k;
      L_w = (xb_vecNx16*)a_l;
      N_w = (xb_vecNx16*)a_n;
      TOL_r = (xb_vecN_2xf32*)a_tol;
      vi = BBE_MOVVA16(i);
      for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
        c0 = BBE_LVN_2XF32_X(D_r, L*sz_f32); /* D[SIDX(i,p) */
        c1 = BBE_LVN_2XF32_X(D_r, L*sz_f32+2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XF32_IP(b0, D_r, 2*BBE_SIMD_WIDTH); /* D[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(b1, D_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a0, F_r, 2*BBE_SIMD_WIDTH); /* F[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(a1, F_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t0, TOL_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t1, TOL_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        a0 = BBE_ABSN_2XF32(a0); a1 = BBE_ABSN_2XF32(a1); 
        b0 = BBE_ABSN_2XF32(b0); b1 = BBE_ABSN_2XF32(b1); 
        c0 = BBE_ABSN_2XF32(c0); c1 = BBE_ABSN_2XF32(c1); 
        /* Zero on the superdiagonal? */
        g0 = BBE_OLEN_2XF32(a0, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(b0,c0)));
        g1 = BBE_OLEN_2XF32(a1, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(b1,c1)));
        /* Zero on the main diagonal? */
        h0 = BBE_OLEN_2XF32(b0,t0); h1 = BBE_OLEN_2XF32(b1,t1);
        /* Tiny superdiagonal elements are flushed to zero. 
         * Condition of i<=k does not matter here! */
        BBE_SVN_2XF32T_IP(BBE_CONSTN_2XF32(0), F_w, 2*BBE_SIMD_WIDTH, g0);
        BBE_SVN_2XF32T_IP(BBE_CONSTN_2XF32(0), F_w, 2*BBE_SIMD_WIDTH, g1);
        f = BBE_LENX16(vi,vk);
        g = BBE_ANDBN(f, BBE_JOINBN_2(g1,g0));
        h = BBE_ANDBN(f, BBE_JOINBN_2(h1,h0));
        BBE_SVNX16T_IP(vi, L_w, 2*BBE_SIMD_WIDTH, g);
        BBE_SVNX16T_IP(vi, N_w, 2*BBE_SIMD_WIDTH, h);
      } /* p */
      if ((L&(BBE_SIMD_WIDTH/2))!=0) {
        c0 = BBE_LVN_2XF32_X(D_r, L*sz_f32); /* D[SIDX(i,p) */
        BBE_LVN_2XF32_IP(b0, D_r, 2*BBE_SIMD_WIDTH); /* D[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(a0, F_r, 2*BBE_SIMD_WIDTH); /* F[SIDX(i-1,p)] */
        BBE_LVN_2XF32_IP(t0, TOL_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        a0 = BBE_ABSN_2XF32(a0);
        b0 = BBE_ABSN_2XF32(b0);
        c0 = BBE_ABSN_2XF32(c0);
        g0 = BBE_OLEN_2XF32(a0, BBE_MULN_2XF32(EPS, BBE_ADDN_2XF32(b0,c0)));
        h0 = BBE_OLEN_2XF32(b0, t0);
        BBE_SVN_2XF32T_IP(BBE_CONSTN_2XF32(0), F_w, 2*BBE_SIMD_WIDTH, g0);
        f = BBE_LENX16(vi,vk);
        g = BBE_ANDBN(f, BBE_JOINBN_2(BBE_MOVN_2_FROMN(BBE_NEQNX16(vi,vi)),g0));
        h = BBE_ANDBN(f, BBE_JOINBN_2(BBE_MOVN_2_FROMN(BBE_NEQNX16(vi,vi)),h0));
        BBE_SVNX16T_IP(vi, L_w, 2*BBE_SIMD_WIDTH, g);
        BBE_SVNX16T_IP(vi, N_w, 2*BBE_SIMD_WIDTH, h);
      } /* L */
      __Pragma("no_reorder");
    } /* i */
#else
    for ( i=n_lo+1; i<=k_up; i++ ) {
      for ( p=0; p<L; p++ ) {
        if (i<=a_k[p]) {
          vbool1 g,h;
          a = XT_ABS_S(F[SIDX(i-1,p)]); b = XT_ABS_S(D[SIDX(i-1,p)]); c = XT_ABS_S(D[SIDX(i,p)]);
          g = XT_OLE_S(a, EPS*(b+c)); h = XT_OLE_S(b, a_tol[p]);
          XT_MOVT_S(F[SIDX(i-1,p)], XT_CONST_S(0), g);
          XT_MOVNEZ(a_l[p], i, BBE_MOVAB1(g)); /* Zero on the superdiagonal? */
          XT_MOVNEZ(a_n[p], i, BBE_MOVAB1(h)); /* Zero on the main diagonal? */
        }
      } /* p */
    } /* i */
#endif
    n_lo = k_up+1;
    for ( p=0; p<L; p++ ) {
      if (a_n[p]>a_l[p]) {
        /* If D(n-1) is zero, set F(n-1) to zero through plane rotations
         * (n-1,n),(n-1,n+1),...,(n-1,k). */
        c = XT_CONST_S(1); s = XT_CONST_S(0);
        for ( j=a_n[p]; j<=a_k[p]; j++ ) {
          a = F[SIDX(j-1,p)]; F[SIDX(j-1,p)] = -s*a;
          a *= c; b = D[SIDX(j,p)];
          /* G(a,b) <- [c,s;-s,c]: [c,s;-s,c]'*[a;b] == [*;0] */
          givens(&c,&s,a,b);
          D[SIDX(j,p)] = c*a-s*b;
          if (U) {
            U_r = U_w = &U[SIDX(a_n[p]-1,p)];
            for ( i=0; i<M; i++ ) {
              xtfloat_loadxp(u, U_r, (j-a_n[p]+1)*L*sz_f32);
              xtfloat_loadxp(v, U_r, (N+a_n[p]-1-j)*L*sz_f32);
              a = s*u+c*v; b = c*u-s*v;
              xtfloat_storexp(a, U_w, (j-a_n[p]+1)*L*sz_f32);
              xtfloat_storexp(b, U_w, (N+a_n[p]-1-j)*L*sz_f32);
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
#else
  float32_t a,b,c,s,u,v;
  int i,j,p,n_lo;

  for ( p=0; p<L; p++ ) {
    a_l[p] = a_n[p] = 0;
  } /* p */
  n_lo = 0;
  do {
    for ( i=n_lo+1; i<=k_up; i++ ) {
      for ( p=0; p<L; p++ ) {
        if (i<=a_k[p]) {
          a = fabsf(F[SIDX(i-1,p)]); b = fabsf(D[SIDX(i-1,p)]); c = fabsf(D[SIDX(i,p)]);
          /* Zero on the superdiagonal? */
          if (a<=EPS*(b+c)) {
            F[SIDX(i-1,p)] = 0.f; a_l[p] = i;
          }
          /* Zero on the main diagonal? */
          if (b<=a_tol[p]) {
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
        c = 1.f; s = 0.f;
        for ( j=a_n[p]; j<=a_k[p]; j++ ) {
          a = F[SIDX(j-1,p)]; F[SIDX(j-1,p)] = -s*a;
          a *= c; b = D[SIDX(j,p)];
          /* G(a,b) <- [c,s;-s,c]: [c,s;-s,c]'*[a;b] == [*;0] */
          givens(&c,&s,a,b);
          D[SIDX(j,p)] = c*a-s*b;
          if (U) {
            for ( i=0; i<M; i++ ) {
              u = U[SIDX(i*N+a_n[p]-1,p)]; v = U[SIDX(i*N+j,p)];
              U[SIDX(i*N+a_n[p]-1,p)] = s*u+c*v;
              U[SIDX(i*N+j,p)] = c*u-s*v;
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
                      const rgrsvdsf_scratch_t * scr,
                      const rgrsvdsf_gksApi_t  * gksApi,
                      float32_t * restrict D, /* [N][L]   (in/out) */
                      float32_t * restrict F, /* [N-1][L] (in/out) */
                      float32_t * restrict U, /* [M*N][L] (in/out) */
                      float32_t * restrict V, /* [N*N][L] (in/out) */
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
   *   t = ((x-y)*(x+y)+(w-v)*(w+v))/(2*x*v);
   *   r = sqrt(t*t+1); if t<0, r = -r; end;
   *   % a <- t11-mu; b <- t12
   *   a = (z-y)*(z+y)+v*(x/(t+r)-v); b = z*f(l);
   *   % QR algorithm iteration
   *   x = d(l); y = f(l); 
   *   for p=l:k-1
   *     q = p+1;
   *     % G(a,b) <- [c,s;-s,c]: [a,b]*[c,s;-s,c] == [*,0]
   *     [c,s] = givens(a,b);
   *     % B <- B*G(a,b)
   *     if p>l, f(p-1) = c*a-s*b; end;
   *     z = d(q);
   *     a = c*x-s*y; b = -s*z;
   *     y = s*x+c*y; z = c*z;
   *     % V <- V*G(a,b)
   *     if ~isempty(V);
   *       u = V(:,p); v = V(:,q);
   *       V(:,p) = c*u-s*v;
   *       V(:,q) = s*u+c*v;
   *     end
   *     % G(a,b) <- [c,s;-s,c]: [c,s;-s,c]'*[a;b] == [*;0]
   *     [c,s] = givens(a,b);
   *     % B <- G(a,b)'*B 
   *     d(p) = c*a-s*b;
   *     a = c*y-s*z; x = s*y+c*z;
   *     if p<k-1, y = f(q); b = -s*y; y = c*y; end;
   *     % U <- U*G(a,b)
   *     if ~isempty(U);
   *       u = U(:,p); v = U(:,q);
   *       U(:,p) = c*u-s*v;
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
    for ( n=l_lo; n<k_up; n++ ) {
      gksApi->givens       (scr->c, scr->s, scr->a, scr->b                         , L);
      gksApi->step1        (D, F, scr->a, scr->b, scr->c, scr->s, scr->l, scr->m, n, L);
      if (V) gksApi->accum (V, scr->c, scr->s, scr->m, n, N, N                     , L);
      gksApi->givens       (scr->c, scr->s, scr->a, scr->b                         , L);
      gksApi->step2        (D, F, scr->a, scr->b, scr->c, scr->s, scr->k, scr->m, n, L);
      if (U) gksApi->accum (U, scr->c, scr->s, scr->m, n, M, N                     , L);
    } /* n */
  } /* l_lo, k_up */
  return (k_up);

} /* golubKahanSVDStep() */

/* Sort singular values in descending order, and permute left- and right-singular
 * vectors (columns of matrices U and V) accordingly. Data are reordered in-place. */
static void sortSingVal(
                       float32_t * restrict D,   /* [N][L]   (in/out) */
                       float32_t * restrict U,   /* [M*N][L] (in/out) */
                       float32_t * restrict V,   /* [N*N][L] (in/out) */
                       int M, int N, int L )
{
  float32_t s;
  int k,l,n,p;
  for ( p=0; p<L; p++ ) {
    /* We use a kind of the bubble sort due to its simplicity. Normally, the number
     * of singular N is too small for an advanced sorting algorithm to demonstrate
     * its power. */
    for ( k=0; k<N-1; k++ ) {
      s = D[SIDX(k,p)];
      for ( l=k, n=k+1; n<N; n++ ) {
        if (D[SIDX(n,p)]>s) { l = n; s = D[SIDX(n,p)]; }
      }
      if (k!=l) {
        s = D[SIDX(l,p)]; D[SIDX(l,p)] = D[SIDX(k,p)]; D[SIDX(k,p)] = s;
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

/* Compute Givens's rotation matrix: given scalars a and b, compute scalars 
 * s and c such that (in MATLAB notation):
 *   A) [a,b]*[c,s;-s,c] == [*,0] and 
 *   B) c^2+s^2 == 1
 * Sanity checks for zero inputs and NaN propagation:
 *   {
 *     float32_t c,s;
 *     givens(&c,&s,qNaNf.f,qNaNf.f); NASSERT(isnan(c) && isnan(s));
 *     givens(&c,&s,  123.f,qNaNf.f); NASSERT(isnan(c) && isnan(s));
 *     givens(&c,&s,qNaNf.f,  456.f); NASSERT(isnan(c) && isnan(s));
 *     givens(&c,&s,    0.f,qNaNf.f); NASSERT(isnan(c) && isnan(s));
 *     givens(&c,&s,qNaNf.f,    0.f); NASSERT(isnan(c) && isnan(s));
 *     givens(&c,&s,  123.f,    0.f); NASSERT(c==1.f && s==0.f);
 *     givens(&c,&s,    0.f, -456.f); NASSERT(c==0.f && s==1.f);
 *     givens(&c,&s,    0.f,    0.f); NASSERT(c==1.f && s==0.f);
 *   }
 * Note that behavior for infinite inputs is not specified.
 */
void givens( float32_t * c, float32_t * s, 
             float32_t a, float32_t b )
{
  /*
   * function [c,s] = givens(a,b)
   * if b~=0
   *   if abs(b)>abs(a)
   *     t = a/b; s = 1/sqrt(1+t^2); c = -t*s;
   *   else
   *     t = b/a; c = 1/sqrt(1+t^2); s = -t*c;
   *   end
   * else
   *   c = 1; s = 0;
   * end;
   */
#if 1
  float32_t f,g,h,t;
  vbool1 bnz,blt;
  /* if (fabsf(a)<fabsf(b) || isnan(a) || isnan(b)) {
   *   f = b; g = a;
   * } else {
   *   f = a; g = b;
   * } */
  blt = XT_ULT_S(XT_ABS_S(a), XT_ABS_S(b));
  f = a; XT_MOVT_S(f,b,blt);
  g = b; XT_MOVT_S(g,a,blt);
  /* if (b==0) f = 1.f; */
  bnz = XT_UNEQ_S(b, XT_CONST_S(0));
  XT_MOVF_S(f, XT_CONST_S(1), bnz);
  /* h <- ~1/f */
  h = XT_RECIP0_S(f);
  /* Newton-Raphson refinement iteration for a reciprocal */
  t = XT_CONST_S(1); XT_MSUBN_S(t,h,f); XT_MADDN_S(h,t,h);
  /* t <- ~g/f */
  t = XT_MUL_S(h,g);
  /* Modified Newton-Raphson iteration for a quotient */
  XT_MSUBN_S(g,t,f); XT_MADDN_S(t,g,h);
  /* f <- 1/sqrtf(1+t*t) */
  h = XT_CONST_S(1); XT_MADDN_S(h,t,t); f = XT_RSQRT_S(h);
  /* g <- -f*t */
  g = XT_MUL_S(XT_NEG_S(t),f);
  /* if (fabsf(a)<fabsf(b) || isnan(a) || isnan(b)) {
   *   *s = f; *c = g;
   * } else {
   *   *s = g; *c = f;
   * } */
  t = f; XT_MOVT_S(f,g,blt); XT_MOVT_S(g,t,blt);
  *s = g; *c = f;
#else
  float32_t t,cc,ss;
  if (!(b==0) || !(a==a)) {
    if (fabsf(b)>fabsf(a)) {
      t = a/b; ss = 1.f/sqrtf(1+t*t); cc = -t*ss;
    } else {
      t = b/a; cc = 1.f/sqrtf(1+t*t); ss = -t*cc;
    }
  } else {
    cc = 1.f; ss = 0.f;
  }
  *s = ss; *c = cc;
#endif
} /* givens() */

#endif /* HAVE_VFPU */

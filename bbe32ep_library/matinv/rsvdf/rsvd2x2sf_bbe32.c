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
    Real 2x2 stream ordered matrices
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <float.h>
#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
/* SVD common declarations */
#include "svd_common.h"
#include "vfpu_math.h"

#define sz_f32   sizeof(float32_t)

#if HAVE_VFPU

#define IT_DOUBLE_FDIVN_2XF32(r1,r0,a,b,c,scaleDivsr)  __IT_DOUBLE_FDIVN_2XF32(&(r1),&(r0),&(a),&(b),&(c),(scaleDivsr))
ATTRIBUTE_ALWAYS_INLINE
inline_ void __IT_DOUBLE_FDIVN_2XF32(void * pr1, void * pr0, const void * pn1, const void * pn0, const void * pd, int scaleDivsr)
#if 0
{
    xb_vecN_2xf32 n0, n1, d;
    xb_vecN_2xf32 r0, r1;
    n0 = *(xb_vecN_2xf32*)pn0;
    n1 = *(xb_vecN_2xf32*)pn1;
    d = *(xb_vecN_2xf32*)pd;

    r0 = IT_FDIVN_2XF32(n0, d, scaleDivsr);
    r1 = IT_FDIVN_2XF32(n1, d, scaleDivsr);

    *(xb_vecN_2xf32*)pr0 = r0;
    *(xb_vecN_2xf32*)pr1 = r1;
} /* __IT_DOUBLE_FDIVN_2XF32() */
#else
{
    xb_vecN_2xf32 n0, n1, d;
    xb_vecN_2xf32 r0, r1;
    n0 = *(xb_vecN_2xf32*)pn0;
    n1 = *(xb_vecN_2xf32*)pn1;
    d  = *(xb_vecN_2xf32*)pd;

    xb_vecN_2xf32 r, t;

    /* t <- ~1/d */
    if (scaleDivsr) {
        /* This variant correctly handles subnormals in the divisor. */
        vboolN_2 blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(d), FLT_MIN);
        BBE_MULN_2XF32T(d, d, 8388608.f, blt);
        t = BBE_RECIP0N_2XF32(d);
        /* Newton-Raphson refinement iteration for a reciprocal */
        r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(r, t, d); BBE_MULAN_2XF32(t, r, t);
        /* r0 <- ~n0/d */
        r = BBE_MULN_2XF32(n0, t);
        /* Modified Newton-Raphson iteration for a quotient */
        BBE_MULSN_2XF32(n0, d, r); BBE_MULAN_2XF32(r, t, n0);
        BBE_MULN_2XF32T(r, r, 8388608.f, blt);
        r0 = r;
        /* r1 <- ~n1/d */
        r = BBE_MULN_2XF32(n1, t);
        /* Modified Newton-Raphson iteration for a quotient */
        BBE_MULSN_2XF32(n1, d, r); BBE_MULAN_2XF32(r, t, n1);
        BBE_MULN_2XF32T(r, r, 8388608.f, blt);
        r1 = r;
    } else {
        /* This variant overflows when the divisor is subnormal. */
        t = BBE_RECIP0N_2XF32(d);
        /* Newton-Raphson refinement iteration for a reciprocal */
        r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(r, t, d); BBE_MULAN_2XF32(t, r, t);
        /* r0 <- ~n0/d */
        r = BBE_MULN_2XF32(n0, t);
        /* Modified Newton-Raphson iteration for a quotient */
        BBE_MULSN_2XF32(n0, d, r); BBE_MULAN_2XF32(r, t, n0);
        r0 = r;
        /* r1 <- ~n1/d */
        r = BBE_MULN_2XF32(n1, t);
        /* Modified Newton-Raphson iteration for a quotient */
        BBE_MULSN_2XF32(n1, d, r); BBE_MULAN_2XF32(r, t, n1);
        r1 = r;
    }
    *(xb_vecN_2xf32*)pr0 = r0;
    *(xb_vecN_2xf32*)pr1 = r1;
} /* __IT_DOUBLE_FDIVN_2XF32() */
#endif

/*-------------------------------------------------------------------------
Thin SVD For Real/Complex Stream Ordered Matrices

Description: compute the Thin Singular Value Decomposition of L complex 
(real) MxN matrices, with the number of rows greater than or equal to the
number of columns: M>=N. Input and output matrices are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Notes:
1. SVD implementation may perform in-place transformations of input matrices,
   so INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once U or V matrix is not required, set the corresponding output pointer
   parameter to zero to allow for a lower complexity implementation of the
   SVD algorithm.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.

Temporary:
  pScr          Scratch area. Required size (in bytes) is defined by 
                functions [r]svd<size>nf_getScratchSize(N,L)
Input:
  M,N           Matrix dimensions
  L             Number of matrices
  A[M*N][L]     MxN input matrices
Output:
  U[M*N][L]     MxN matrices comprised of M left-singular column 
                vectors If NULL, this agrument is assumed to be optional 
                and will not be computed
  s[N][L]       Nx1 vectors of singular values in descending order. In an 
                exceptional case when the iterative algorithm fails to 
                converge for a particular matrix, all elements of the
                respective vector are set to NaN.
  V[N*N][L]     NxN matrices comprised of N right-singular column
                vectors. If NULL, this agrument is assumed to be optional 
                and will not be computed
Restrictions:
  pScr,U,s,V,A  Must not overlap and must be aligned on 32-byte boundary 
  M,N           M>=N>1
  L             Must be a multiple of 8
---------------------------------------------------------------------------*/
void rsvd2x2sf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L )
{
  /*
   * MATLAB reference code:
   *
   *   function [U,S,V] = svd2x2(A)
   *   [M,N] = size(A); 
   *   assert((M==2)&&(N==2));
   *   rmin = realmin(class(A));
   *   tol = eps(cast(1,'like',A));
   *   V = eye(2);
   *   % Perform two Jacobi rotations of column vectors. The first rotation does
   *   % the vast of orthogonalization work, and the second rotation cleans small
   *   % errors due to round-off in sine/cosine computation. 
   *   for n=1:3
   *     % Compute a sine/consine pair. The rotation angle is selected in such a
   *     % way that the L2 norm of the second column does not exceed the norm of 
   *     % the first one.
   *     x = A(:,1); y = A(:,2); 
   *     p = x'*y; q = x'*x; r = y'*y;
   *     if q>=r
   *       if cmag(p)>tol*r
   *         f = p/q; g = 1-r/q;
   *         h = sqrt(4*f*conj(f)+g^2);
   *         cs = sqrt((1+g/h)/2);
   *         sn = -f/(h*cs);
   *       else
   *         cs = 1; sn = 0;
   *       end
   *     else
   *       f = p/r; g = 1-q/r;
   *       h = sqrt(4*f*conj(f)+g^2);
   *       sn = sqrt((1+g/h)/2);
   *       cs = -f/(h*sn);
   *     end
   *     % Apply the rotator to working matrix
   *     A = A*[cs,sn;-sn',cs'];
   *     % Accumulate transformations
   *     V = V*[cs,sn;-sn',cs'];
   *   end
   *   % Compute the first singular value and the corresponding left-singular
   *   % vector.
   *   s = A(:,1)'*A(:,1);
   *   if s<rmin
   *     % The first singular value is either zero or too small to be evaluated to
   *     % working precision. Note that the second singular value cannot exceed
   *     % the first one.
   *     s1 = 0; u1 = [1;0];
   *   else
   *     s1 = sqrt(s); u1 = A(:,1)/s1;
   *   end
   *   % Compute the first singular value and the corresponding left-singular
   *   % vector.
   *   s = A(:,2)'*A(:,2);
   *   if s<rmin
   *     % The second singular value is either zero or tiny. Select the orthogonal
   *     % complement of the first left-singular vector for the second singular
   *     % vector.
   *     s2 = 0; u2 = [-u1(2);u1(1)];
   *   else
   *     s2 = sqrt(s); u2 = A(:,2)/s2;
   *   end
   *   U = [u1,u2];
   *   S = diag([s1,s2]);
   */

  const xb_vecN_2xf32 * restrict A_r;
  const xb_vecN_2xf32 * P_r;
  const xb_vecN_2xf32 * Q_r;
  const xb_vecN_2xf32 * R_r;
  const xb_vecN_2xf32 * F_r;
  const xb_vecN_2xf32 * G_r;
  const xb_vecN_2xf32 * H_r;
  const xb_vecN_2xf32 * SN_r;
  const xb_vecN_2xf32 * restrict CS_r;
  const xb_vecN_2xf32 * restrict V_r;
  const xb_vecN_2xf32 * restrict U_r;
        xb_vecN_2xf32 * restrict A_w;
        xb_vecN_2xf32 * P_w;
        xb_vecN_2xf32 * Q_w;
        xb_vecN_2xf32 * R_w;
        xb_vecN_2xf32 * S_w;
        xb_vecN_2xf32 * F_w;
        xb_vecN_2xf32 * G_w;
        xb_vecN_2xf32 * H_w;
        xb_vecN_2xf32 * SN_w;
        xb_vecN_2xf32 * restrict CS_w;
        xb_vecN_2xf32 * restrict V_w;
        xb_vecN_2xf32 * restrict U_w;

  float32_t * restrict P;
  float32_t * restrict Q;
  float32_t * restrict R;
  float32_t * restrict F;
  float32_t * restrict G;
  float32_t * restrict H;
  float32_t * restrict SN;
  float32_t * restrict CS;
#if 0
  float32_t f,g,h,p,q,r;
  float32_t cs,sn;
#endif
  int k,n;

  const float32_t tol     = FLT_EPSILON;
  const float32_t rmin    = FLT_MIN;
  const int       iterNum = 3;

  xb_vecN_2xf32 a00, a01, a10, a11;
  xb_vecN_2xf32 v00, v01, v10, v11;
  xb_vecN_2xf32 u00, u01, u10, u11;
  xb_vecN_2xf32 _p, _q, _r, _f, _g, _h;
  xb_vecN_2xf32 _cs, _sn;
  xb_vecN_2xf32 _t;
  vboolN_2 bflrmin;
  vboolN_2 bqlr, bpler;
  vboolN_2 bqmer;

  #define IDX(i,j,k)    (((i)*2+(j))*L+(k))

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( U   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( s   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  {
    void * p = pScr;
    /* Partition the scratch area. */
    P  = (float32_t*)p; p = P + L;
    Q  = (float32_t*)p; p = Q + L;
    R  = (float32_t*)p; p = R + L;
    F  = (float32_t*)p; p = F + L;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)rsvd2x2sf_getScratchSize(2,2,L) );
    /* Partially reuse the scratch area. */
    G = &s[IDX(0,0,0)]; H = &s[IDX(0,1,0)]; SN = P; CS = F;
  }

  /* Initialize the transformation matrix. */
  if (V) {
    for ( k=0; k<L; k++ ) {
      V[IDX(0,0,k)] = 1.f; V[IDX(0,1,k)] = 0.f;
      V[IDX(1,0,k)] = 0.f; V[IDX(1,1,k)] = 1.f;
    }
  }
  /* Perform a number of Jacobi rotations of column vectors. The first rotation
   * does the vast of orthogonalization work, and the second rotation cleans
   * small errors due to round-off in sine/cosine computation. */
  for ( n=0; n<iterNum; n++ ) {
    /* Compute a sine/consine pair, part 1. */
#if 0
    for ( k=0; k<L; k++ ) {
      P[k] = A[IDX(0,0,k)]*A[IDX(0,1,k)] + A[IDX(1,0,k)]*A[IDX(1,1,k)];
      Q[k] = A[IDX(0,0,k)]*A[IDX(0,0,k)] + A[IDX(1,0,k)]*A[IDX(1,0,k)];
      R[k] = A[IDX(0,1,k)]*A[IDX(0,1,k)] + A[IDX(1,1,k)]*A[IDX(1,1,k)];
    } /* k */
#else
    A_r = (xb_vecN_2xf32 *)A;
    P_w = (xb_vecN_2xf32 *)P;
    Q_w = (xb_vecN_2xf32 *)Q;
    R_w = (xb_vecN_2xf32 *)R;
    for (k = 0; k < L/(BBE_SIMD_WIDTH/2); k++) {
        BBE_LVN_2XF32_XP(a00, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a01, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a10, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 3*L)*sz_f32);
        _p = BBE_MULN_2XF32(a00, a01); BBE_MULAN_2XF32(_p, a10,a11);
        _q = BBE_MULN_2XF32(a00, a00); BBE_MULAN_2XF32(_q, a10, a10);
        _r = BBE_MULN_2XF32(a01, a01); BBE_MULAN_2XF32(_r, a11, a11);
        BBE_SVN_2XF32_IP(_p, P_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_q, Q_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_r, R_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
#endif
    /* Compute a sine/consine pair, part 2. The rotation angle is selected in
     * such a way that the L2 norm of the second column does not exceed the
     * norm of the first one. */
#if 0
    for ( k=0; k<L; k++ ) {
      p = P[k]; q = Q[k]; r = R[k];
      if (q>=r) {
        f = p/q; g = 1-r/q;
      } else {
        f = p/r; g = 1-q/r;
      }
      F[k] = f; G[k] = g;
      H[k] = sqrtf(4*f*f+g*g);
    } /* k */
#elif 0
    P_r = (xb_vecN_2xf32 *)P;
    Q_r = (xb_vecN_2xf32 *)Q;
    R_r = (xb_vecN_2xf32 *)R;
    F_w = (xb_vecN_2xf32 *)F;
    G_w = (xb_vecN_2xf32 *)G;
    H_w = (xb_vecN_2xf32 *)H;

    for (k = 0; k < L/(BBE_SIMD_WIDTH/2); k++) {
        BBE_LVN_2XF32_IP(_p, P_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_r, R_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bqmer = BBE_OLEN_2XF32(_r,_q);
        _t = _r;
        _r = BBE_MOVN_2XF32T(_r, _q, bqmer);
        _q = BBE_MOVN_2XF32T(_q, _t, bqmer);

        IT_DOUBLE_FDIVN_2XF32(_g, _f, _r, _p, _q, 1);
        _g = BBE_SUBN_2XF32(BBE_CONSTN_2XF32(1), _g);

        _h = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), _f);
        _h = BBE_MULN_2XF32(_h, _h); BBE_MULAN_2XF32(_h, _g, _g);
        _h = BBE_SQRTN_2XF32(_h);

        BBE_SVN_2XF32_IP(_f, F_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_g, G_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_h, H_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
#elif 0
    P_r = (xb_vecN_2xf32 *)P;
    Q_r = (xb_vecN_2xf32 *)Q;
    R_r = (xb_vecN_2xf32 *)R;
    F_w = (xb_vecN_2xf32 *)F;
    G_w = (xb_vecN_2xf32 *)G;
    H_w = (xb_vecN_2xf32 *)H;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        BBE_LVN_2XF32_IP(_p, P_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_r, R_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bqmer = BBE_OLEN_2XF32(_r, _q);
        _t = _r;
        _r = BBE_MOVN_2XF32T(_r, _q, bqmer);
        _q = BBE_MOVN_2XF32T(_q, _t, bqmer);

        IT_DOUBLE_FDIVN_2XF32(_g, _f, _r, _p, _q, 0);
        _g = BBE_SUBN_2XF32(BBE_CONSTN_2XF32(1), _g);

        _h = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), _f);
        _h = BBE_MULN_2XF32(_h, _h); BBE_MULAN_2XF32(_h, _g, _g);
        _h = BBE_FSQRTN_2XF32(_h);

        BBE_SVN_2XF32_IP(_f, F_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_g, G_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_h, H_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
#else
    P_r = (xb_vecN_2xf32 *)P;
    Q_r = (xb_vecN_2xf32 *)Q;
    R_r = (xb_vecN_2xf32 *)R;
    F_w = (xb_vecN_2xf32 *)F;
    G_w = (xb_vecN_2xf32 *)G;
    H_w = (xb_vecN_2xf32 *)H;
    for (k = 0; k < L/(BBE_SIMD_WIDTH/2); k++) {
        BBE_LVN_2XF32_IP(_p, P_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_r, R_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bqmer = BBE_OLEN_2XF32(_r,_q);
        _t = _r;
        _r = BBE_MOVN_2XF32T(_r, _q, bqmer);
        _q = BBE_MOVN_2XF32T(_q, _t, bqmer);

        IT_DOUBLE_FDIVN_2XF32(_g, _f, _r, _p, _q, 0);
        _g = BBE_SUBN_2XF32(BBE_CONSTN_2XF32(1), _g);

        BBE_SVN_2XF32_IP(_f, F_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_g, G_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
    F_w = (xb_vecN_2xf32 *)F;
    G_w = (xb_vecN_2xf32 *)G;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        BBE_LVN_2XF32_IP(_f, F_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_g, G_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        _h = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), _f);
        _h = BBE_MULN_2XF32(_h, _h); BBE_MULAN_2XF32(_h, _g, _g);
        _h = BBE_FSQRTN_2XF32(_h);
        BBE_SVN_2XF32_IP(_h, H_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
#endif
    /* Compute a sine/consine pair, part 3. */
#if 0
    for ( k=0; k<L; k++ ) {
      p = P[k]; q = Q[k]; r = R[k];
      f = F[k]; g = G[k]; h = H[k];
      cs = sqrtf((1+g/h)/2);
      sn = -f/(h*cs);
      if (q<r) {
        f = cs; cs = sn; sn = f;
      } else if (fabs(p)<=tol*r) {
        cs = 1.f; sn = 0.f;
      }
      SN[k] = sn; CS[k] = cs;
    } /* k */
#elif 0
    P_r = (xb_vecN_2xf32 *)P;
    Q_r = (xb_vecN_2xf32 *)Q;
    R_r = (xb_vecN_2xf32 *)R;
    F_r = (xb_vecN_2xf32 *)F;
    G_r = (xb_vecN_2xf32 *)G;
    H_r = (xb_vecN_2xf32 *)H;
    SN_w = (xb_vecN_2xf32 *)SN;
    CS_w = (xb_vecN_2xf32 *)CS;

    for (k = 0; k<L /(BBE_SIMD_WIDTH / 2); k++) {
        BBE_LVN_2XF32_IP(_p, P_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_r, R_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_f, F_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_g, G_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_h, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bqlr = BBE_OLTN_2XF32(_q,_r);

        _t = IT_FDIVN_2XF32(_g, _h, 1); _cs = BBE_CONSTN_2XF32(3);
        BBE_MULAN_2XF32(_cs, BBE_CONSTN_2XF32(3), _t); _cs = BBE_SQRTN_2XF32(_cs);
        _sn = BBE_MULN_2XF32(_h, _cs); _sn = BBE_NEGN_2XF32(IT_FDIVN_2XF32(_f, _sn, 1));

        _f = _cs;
        _cs = BBE_MOVN_2XF32T(_sn, _cs, bqlr); _sn = BBE_MOVN_2XF32T(_f, _sn, bqlr);
        bpler = BBE_OLEN_2XF32(BBE_ABSN_2XF32(_p),BBE_MULN_2XF32(tol,_r));
        _cs = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), _cs, bpler); _sn = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _sn, bpler);

        BBE_SVN_2XF32_IP(_sn, SN_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_cs, CS_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
#else
    P_r = (xb_vecN_2xf32 *)P;
    Q_r = (xb_vecN_2xf32 *)Q;
    R_r = (xb_vecN_2xf32 *)R;
    F_r = (xb_vecN_2xf32 *)F;
    G_r = (xb_vecN_2xf32 *)G;
    H_r = (xb_vecN_2xf32 *)H;
    SN_w = (xb_vecN_2xf32 *)SN;
    CS_w = (xb_vecN_2xf32 *)G;
    for (k = 0; k<L /(BBE_SIMD_WIDTH / 2); k++) {
        BBE_LVN_2XF32_IP(_g, G_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_h, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
#if 1
        _t = IT_FDIVN_2XF32(_g, _h, 0); _cs = BBE_CONSTN_2XF32(3);
        BBE_MULAN_2XF32(_cs, BBE_CONSTN_2XF32(3), _t); _cs = BBE_FSQRTN_2XF32(_cs);
#else
        _t = IT_FDIVN_2XF32(_g, _h, 1); _cs = BBE_CONSTN_2XF32(3);
        BBE_MULAN_2XF32(_cs, BBE_CONSTN_2XF32(3), _t); _cs = BBE_SQRTN_2XF32(_cs);
#endif
        BBE_SVN_2XF32_IP(_cs, CS_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    } /* k */
    CS_r = (xb_vecN_2xf32 *)G;
    H_r  = (xb_vecN_2xf32 *)H;
    CS_w = (xb_vecN_2xf32 *)CS;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        BBE_LVN_2XF32_IP(_p, P_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_r, R_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_f, F_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_h, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_cs, CS_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bqlr = BBE_OLTN_2XF32(_q, _r);
        _sn = BBE_MULN_2XF32(_h, _cs); _sn = BBE_NEGN_2XF32(IT_FDIVN_2XF32(_f, _sn, 1));
        _f = _cs;
        _cs = BBE_MOVN_2XF32T(_sn, _cs, bqlr); _sn = BBE_MOVN_2XF32T(_f, _sn, bqlr);
        bpler = BBE_OLEN_2XF32(BBE_ABSN_2XF32(_p), BBE_MULN_2XF32(tol, _r));
        _cs = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), _cs, bpler); _sn = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _sn, bpler);

        BBE_SVN_2XF32_IP(_sn, SN_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_cs, CS_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
    }
    __Pragma("no_reorder");
#endif
    /* Apply the rotator to working matrix: A <- A*[cs,sn;-sn',cs'] */
#if 0
    for ( k=0; k<L; k++ ) {
      sn = SN[k]; cs = CS[k];
      f = A[IDX(0,0,k)]; g = A[IDX(0,1,k)];
      p = A[IDX(1,0,k)]; q = A[IDX(1,1,k)];
      A[IDX(0,0,k)] = f*cs-g*sn; A[IDX(0,1,k)] = f*sn+g*cs;
      A[IDX(1,0,k)] = p*cs-q*sn; A[IDX(1,1,k)] = p*sn+q*cs;
    } /* k */
#else
    A_r = A_w = (xb_vecN_2xf32 *)A;
    SN_r = (xb_vecN_2xf32 *)SN;
    CS_r = (xb_vecN_2xf32 *)CS;
    for (k = 0; k<L/(BBE_SIMD_WIDTH/2); k++) {
        BBE_LVN_2XF32_IP(_sn, SN_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_cs, CS_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        BBE_LVN_2XF32_XP(a00, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a01, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a10, A_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);

        _f = a00; _g = a01; _p = a10; _q = a11;

        a00 = BBE_MULN_2XF32(_f, _cs); BBE_MULSN_2XF32(a00, _g, _sn);
        a01 = BBE_MULN_2XF32(_f, _sn); BBE_MULAN_2XF32(a01, _g, _cs);
        a10 = BBE_MULN_2XF32(_p, _cs); BBE_MULSN_2XF32(a10, _q, _sn);
        a11 = BBE_MULN_2XF32(_p, _sn); BBE_MULAN_2XF32(a11, _q, _cs);

        BBE_SVN_2XF32_XP(a00, A_w, 1 * L*sz_f32);
        BBE_SVN_2XF32_XP(a01, A_w, 1 * L*sz_f32);
        BBE_SVN_2XF32_XP(a10, A_w, 1 * L*sz_f32);
        BBE_SVN_2XF32_XP(a11, A_w, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
    } /* k */
    __Pragma("no_reorder");
#endif
    /* Accumulate transformations: V <- V*[cs,sn;-sn',cs'] */
    if (V) {
#if 0
      for ( k=0; k<L; k++ ) {
        sn = SN[k]; cs = CS[k];
        f = V[IDX(0,0,k)]; g = V[IDX(0,1,k)];
        p = V[IDX(1,0,k)]; q = V[IDX(1,1,k)];
        V[IDX(0,0,k)] = f*cs-g*sn;
        V[IDX(0,1,k)] = f*sn+g*cs;
        V[IDX(1,0,k)] = p*cs-q*sn;
        V[IDX(1,1,k)] = p*sn+q*cs;
      } /* k */
#else
      V_r = V_w = (xb_vecN_2xf32 *)V;
      SN_r = (xb_vecN_2xf32 *)SN;
      CS_r = (xb_vecN_2xf32 *)CS;
      for (k = 0; k<L/(BBE_SIMD_WIDTH/2); k++) {
          BBE_LVN_2XF32_IP(_sn, SN_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
          BBE_LVN_2XF32_IP(_cs, CS_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

          BBE_LVN_2XF32_XP(v00, V_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(v01, V_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(v10, V_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(v11, V_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);

          _f = v00; _g = v01; _p = v10; _q = v11;
          v00 = BBE_MULN_2XF32(_f, _cs); BBE_MULSN_2XF32(v00, _g, _sn);
          v01 = BBE_MULN_2XF32(_f, _sn); BBE_MULAN_2XF32(v01, _g, _cs);
          v10 = BBE_MULN_2XF32(_p, _cs); BBE_MULSN_2XF32(v10, _q, _sn);
          v11 = BBE_MULN_2XF32(_p, _sn); BBE_MULAN_2XF32(v11, _q, _cs);

          BBE_SVN_2XF32_XP(v00, V_w, 1 * L*sz_f32);
          BBE_SVN_2XF32_XP(v01, V_w, 1 * L*sz_f32);
          BBE_SVN_2XF32_XP(v10, V_w, 1 * L*sz_f32);
          BBE_SVN_2XF32_XP(v11, V_w, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
      } /* k */
      __Pragma("no_reorder");
#endif
    } /* V */
  } /* n */
  /* Compute the first singular value and the corresponding left-singular vector. */
#if 0
  for ( k=0; k<L; k++ ) {
    f = A[IDX(0,0,k)]*A[IDX(0,0,k)] + A[IDX(1,0,k)]*A[IDX(1,0,k)];
    /* A singular value may be zero or too small to be evaluated to working precision */
    g = ( f<rmin ? 0.f : sqrtf(f) ); s[IDX(0,0,k)] = g;
    if (U) {
      U[IDX(0,0,k)] = ( f<rmin ? 1.f : A[IDX(0,0,k)]/g );
      U[IDX(1,0,k)] = ( f<rmin ? 0.f : A[IDX(1,0,k)]/g );
    }
  } /* k */
#elif 0
  if (U) {
      U_r = U_w = (xb_vecN_2xf32 *)U;
      for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a00, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a10, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a00, a00); BBE_MULAN_2XF32(_f, a10, a10);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);
          _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
          BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
         
          BBE_LVN_2XF32_XP(u00, U_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(u10, U_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          IT_DOUBLE_FDIVN_2XF32(u10,u00,a10,a00,_g,1);
          u00 = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), u00, bflrmin);
          u10 = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), u10, bflrmin);

          BBE_SVN_2XF32_XP(u00, U_w, 2 * L*sz_f32);
          BBE_SVN_2XF32_XP(u10, U_w, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);
      } /* k */
      __Pragma("no_reorder");
  } else {
      for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a00, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a10, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a00, a00); BBE_MULAN_2XF32(_f, a10, a10);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);
          _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
          BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
      } /* k */
  }
#else
  A_r = (xb_vecN_2xf32 *)A;
  S_w = (xb_vecN_2xf32 *)s;
  for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
      BBE_LVN_2XF32_XP(a00, A_r, 2 * L*sz_f32);
      BBE_LVN_2XF32_XP(a10, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

      _f = BBE_MULN_2XF32(a00, a00); BBE_MULAN_2XF32(_f, a10, a10);
      bflrmin = BBE_OLTN_2XF32(_f, rmin);
#if 1
      _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
#else
      _g = BBE_FSQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
#endif
      BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
  } /* k */

  if(U){
      A_r = (xb_vecN_2xf32 *)A;
      S_w = (xb_vecN_2xf32 *)s;
      U_r = U_w = (xb_vecN_2xf32 *)U;
      for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a00, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a10, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a00, a00); BBE_MULAN_2XF32(_f, a10, a10);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);

          BBE_LVN_2XF32_XP(u00, U_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(u10, U_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          BBE_LVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
#if 1
          IT_DOUBLE_FDIVN_2XF32(u10, u00, a10, a00, _g, 0);
#else
          IT_DOUBLE_FDIVN_2XF32(u10, u00, a10, a00, _g, 1);
#endif
          u00 = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), u00, bflrmin);
          u10 = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), u10, bflrmin);

          BBE_SVN_2XF32_XP(u00, U_w, 2 * L*sz_f32);
          BBE_SVN_2XF32_XP(u10, U_w, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);
      }
      __Pragma("no_reorder");
  }
#endif
  /* Compute the first singular value and the corresponding left-singular vector. */
#if 0
  for ( k=0; k<L; k++ ) {
    f = A[IDX(0,1,k)]*A[IDX(0,1,k)] + A[IDX(1,1,k)]*A[IDX(1,1,k)];
    /* Note that the second singular value cannot exceed the first one. */
    g = ( f<rmin ? 0.f : sqrtf(f) ); s[IDX(0,1,k)] = g;
    if (U) {
      /* If the second singular value is zero or too small, select the orthogonal
       * complement of the first left-singular vector for the second left-singular
       * vector. */
      U[IDX(0,1,k)] = ( f<rmin ? -U[IDX(1,0,k)] : A[IDX(0,1,k)]/g );
      U[IDX(1,1,k)] = ( f<rmin ?  U[IDX(0,0,k)] : A[IDX(1,1,k)]/g );
    }
  } /* k */
#elif 0
  A_r = (xb_vecN_2xf32 *)&A[L];
  S_w = (xb_vecN_2xf32 *)&s[L];
  if (U){
      U_r = (xb_vecN_2xf32 *)U;
      U_w = (xb_vecN_2xf32 *)&U[L];
      for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a01, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a01, a01); BBE_MULAN_2XF32(_f, a11, a11);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);
          _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
          BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
          
          BBE_LVN_2XF32_XP(u00, U_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(u01, U_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(u10, U_r, 1 * L*sz_f32);
          BBE_LVN_2XF32_XP(u11, U_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
          
          u01 = IT_FDIVN_2XF32(a01, _g, 1); u01 = BBE_MOVN_2XF32T(BBE_NEGN_2XF32(u10), u01, bflrmin);
          u11 = IT_FDIVN_2XF32(a11, _g, 1); u11 = BBE_MOVN_2XF32T(u00, u11, bflrmin);
          
          BBE_SVN_2XF32_XP(u01, U_w, 2 * L*sz_f32);
          BBE_SVN_2XF32_XP(u11, U_w, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);
      } /* k */
      __Pragma("no_reorder");
  } else {
      for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a01, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a01, a01); BBE_MULAN_2XF32(_f, a11, a11);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);
          _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
          BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
      } /* k */
  }
#else
  A_r = (xb_vecN_2xf32 *)&A[L];
  S_w = (xb_vecN_2xf32 *)&s[L];
  for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
      BBE_LVN_2XF32_XP(a01, A_r, 2 * L*sz_f32);
      BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

      _f = BBE_MULN_2XF32(a01, a01); BBE_MULAN_2XF32(_f, a11, a11);
      bflrmin = BBE_OLTN_2XF32(_f, rmin);
#if 1
      _g = BBE_SQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
#else
      _g = BBE_FSQRTN_2XF32(_f); _g = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), _g, bflrmin);
#endif
      BBE_SVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
  } /* k */
  if(U){
      A_r = (xb_vecN_2xf32 *)&A[L];
      S_w = (xb_vecN_2xf32 *)&s[L];
      U_r = (xb_vecN_2xf32 *)U;
      U_w = (xb_vecN_2xf32 *)&U[L];
      for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
          BBE_LVN_2XF32_XP(a01, A_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_XP(a11, A_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);

          _f = BBE_MULN_2XF32(a01, a01); BBE_MULAN_2XF32(_f, a11, a11);
          bflrmin = BBE_OLTN_2XF32(_f, rmin);
          BBE_LVN_2XF32_IP(_g, S_w, (BBE_SIMD_WIDTH / 2)*sz_f32);

          u01 = BBE_LVN_2XF32_X(U_r, 1 * L*sz_f32);
          u11 = BBE_LVN_2XF32_X(U_r, 3 * L*sz_f32);
#if 1
          IT_DOUBLE_FDIVN_2XF32(u11, u01, a11, a01, _g, 0);
#else
          IT_DOUBLE_FDIVN_2XF32(u11, u01, a11, a01, _g, 1);
#endif
          u10 = BBE_LVN_2XF32_X(U_r, 2 * L*sz_f32);
          BBE_LVN_2XF32_IP(u00, U_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

          u01 = BBE_MOVN_2XF32T(BBE_NEGN_2XF32(u10), u01, bflrmin);
          u11 = BBE_MOVN_2XF32T(u00, u11, bflrmin);

          BBE_SVN_2XF32_XP(u01, U_w, 2 * L*sz_f32);
          BBE_SVN_2XF32_XP(u11, U_w, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);
      } /* k */
      __Pragma("no_reorder");
  }
#endif

  #undef IDX

} /* rsvd2x2sf() */

size_t rsvd2x2sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 2==M && 2==N ); 
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return ( L*sz_f32 +  /* P: dot products               */
           L*sz_f32 +  /* Q: squared norm of 1st column */
           L*sz_f32 +  /* R: squared norm of 2nd column */
           L*sz_f32 ); /* F: normalized dot products    */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, rsvd2x2sf, ( void * pScr,
                       float32_t * restrict U,
                       float32_t * restrict s,
                       float32_t * restrict V,
                       float32_t * restrict A,
                       int L ) )

size_t rsvd2x2sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 2==M && 2==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (0);
}

#endif /* HAVE_VFPU */

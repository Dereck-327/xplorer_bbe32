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
    Complex 2x2 stream ordered matrices
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
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
/* SVD common declarations */
#include "svd_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

#if HAVE_VFPU
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
void svd2x2sf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
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
   *     s2 = 0; u2 = [-conj(u1(2));conj(u1(1))];
   *   else
   *     s2 = sqrt(s); u2 = A(:,2)/s2;
   *   end
   *   U = [u1,u2];
   *   S = diag([s1,s2]);
   */

        xb_vecN_2xf32 * restrict U_w;
  const xb_vecN_2xf32 *          S_r;
        xb_vecN_2xf32 * restrict S_w;
        xb_vecN_2xf32 * restrict V_w0;
        xb_vecN_2xf32 * restrict V_w1;
        xb_vecN_2xf32 * restrict V_w2;
        xb_vecN_2xf32 * restrict V_w3;
  const xb_vecN_2xf32 *          V_r;
        xb_vecN_2xf32 * restrict A_w0;
        xb_vecN_2xf32 * restrict A_w1;
        xb_vecN_2xf32 * restrict A_w2;
        xb_vecN_2xf32 * restrict A_w3;
  const xb_vecN_2xf32 *          A_r;
  const xb_vecN_2xf32 *          A_r0;
  const xb_vecN_2xf32 *          A_r1;
        xb_vecN_2xf32 * restrict P_w;
  const xb_vecN_2xf32 *          P_r;
        xb_vecN_2xf32 * restrict Q_w;
  const xb_vecN_2xf32 *          Q_r;
        xb_vecN_2xf32 * restrict R_w;
  const xb_vecN_2xf32 *          R_r;
        xb_vecN_2xf32 * restrict F_w;
  const xb_vecN_2xf32 *          F_r;
        xb_vecN_2xf32 * restrict G_w;
  const xb_vecN_2xf32 *          G_r;
        xb_vecN_2xf32 * restrict H_w;
  const xb_vecN_2xf32 *          H_r;
        xb_vecN_2xf32 * restrict SN_w;
  const xb_vecN_2xf32 *          SN_r;
        xb_vecN_2xf32 * restrict CS_w;
  const xb_vecN_2xf32 *          CS_r;

  complex_float * restrict P;
  float32_t     * restrict Q;
  float32_t     * restrict R;
  complex_float * restrict F;
  float32_t     * restrict G;
  float32_t     * restrict H;
  complex_float * restrict SN;
  complex_float * restrict CS;
#if 0
  const complex_float c0f = _makecomplexf(0.f,0.f);
  const complex_float c1f = _makecomplexf(1.f,0.f);

  const float32_t tol  = FLT_EPSILON;
  const float32_t rmin = FLT_MIN;

  float32_t g,h,q,r,t;
  complex_float a,b,c,d,f,p;
  complex_float cs,sn;
#endif
  const xb_vecN_2xf32 _c0f = BBE_CONSTN_2XF32(0);
  const xb_vecN_2xf32 _c1f = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0), BBE_CONSTN_2XF32(1), BBE_SELI_INTERLEAVE_2_LO);

  const xb_vecN_2xf32 _tol = FLT_EPSILON;
  const xb_vecN_2xf32 _rmin = FLT_MIN;
  const int            iterNum = 3;

  xb_vecN_2xf32 _g,_h,_h0,_h1,_q,_r,_t;
  xb_vecN_2xf32 _a,_b,_c,_d,_f,_f0,_f1,_p;
  xb_vecN_2xf32 _cs,_sn;
  vboolN_2 bge,blt,ble,beq,bneq;
  
  int k,n;

  #define IDX(i,j,k)    (((i)*2+(j))*L+(k))

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( U   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( s   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  {
    void * p = pScr;
    /* Partition the scratch area. All real data are duplicated to simplify mixed 
     * real/complex computations */
    P  = (complex_float*)p; p = P + L;
    Q  = (float32_t    *)p; p = Q + 2*L;
    R  = (float32_t    *)p; p = R + 2*L;
    F  = (complex_float*)p; p = F + L;
    G  = (float32_t    *)p; p = G + 2*L;
    /* Make sure that scratch arrays fit into the reserved space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)svd2x2sf_getScratchSize(2,2,L) );
    /* Partially reuse the scratch area. */
    H = &s[IDX(0,0,0)]; SN = P; CS = F;
  }

  /* Initialize the transformation matrix. */
#if 1
  if (V) {
    V_w0 = (xb_vecN_2xf32*)&V[IDX(0,0,0)]; V_w1 = (xb_vecN_2xf32*)&V[IDX(0,1,0)]; 
    V_w2 = (xb_vecN_2xf32*)&V[IDX(1,0,0)]; V_w3 = (xb_vecN_2xf32*)&V[IDX(1,1,0)]; 
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_SVN_2XF32_IP(_c1f, V_w0, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_c0f, V_w1, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_c0f, V_w2, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_c1f, V_w3, 2*BBE_SIMD_WIDTH);
    }
  }
  __Pragma("no_reorder");
#else
  if (V) {
    for ( k=0; k<L; k++ ) {
      V[IDX(0,0,k)] = c1f; V[IDX(0,1,k)] = c0f;
      V[IDX(1,0,k)] = c0f; V[IDX(1,1,k)] = c1f;
    }
  }
#endif

  /* Perform a number of Jacobi rotations of column vectors. The first rotation
   * does the vast of orthogonalization work, and the second rotation cleans
   * small errors due to round-off in sine/cosine computation. */
  for ( n=0; n<iterNum; n++ ) {
    /* Compute a sine/consine pair, part 1. */
#if 1
    A_r = (xb_vecN_2xf32*)A;
    P_w = (xb_vecN_2xf32*)P;
    Q_w = (xb_vecN_2xf32*)Q;
    R_w = (xb_vecN_2xf32*)R;
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_LVN_2XF32_XP(_a, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_c, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_d, A_r, 2*BBE_SIMD_WIDTH-3*L*sz_f32c);
      /* p = caddf(cmuljf(A[IDX(0,1,k)], A[IDX(0,0,k)]),
       *           cmuljf(A[IDX(1,1,k)], A[IDX(1,0,k)])); */
      _p = BBE_MULMN_2XF32(_a,_b,2,8);
      BBE_MULMASN_2XF32(_p,_a,_b,0,7);
      BBE_MULMASN_2XF32(_p,_c,_d,2,8);
      BBE_MULMASN_2XF32(_p,_c,_d,0,7);
      /* q = cabs2f(A[IDX(0,0,k)]) + cabs2f(A[IDX(1,0,k)]); */
      _q = BBE_MULMN_2XF32(_a,_a,0,12);
      BBE_MULMASN_2XF32(_q,_c,_c,0,12);
      _q = BBE_ADDN_2XF32(_q, BBE_SHFLN_2XF32I(_q, BBE_SHFLI_SWAP_2));
      /* r = cabs2f(A[IDX(0,1,k)]) + cabs2f(A[IDX(1,1,k)]); */
      _r = BBE_MULMN_2XF32(_b,_b,0,12);
      BBE_MULMASN_2XF32(_r,_d,_d,0,12);
      _r = BBE_ADDN_2XF32(_r, BBE_SHFLN_2XF32I(_r, BBE_SHFLI_SWAP_2));
      BBE_SVN_2XF32_IP(_p, P_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_q, Q_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_r, R_w, 2*BBE_SIMD_WIDTH);
    } /* k */
    __Pragma("no_reorder");
#else
    for ( k=0; k<L; k++ ) {
      p = caddf(cmuljf(A[IDX(0,1,k)], A[IDX(0,0,k)]),
                cmuljf(A[IDX(1,1,k)], A[IDX(1,0,k)]));
      q = cabs2f(A[IDX(0,0,k)]) + cabs2f(A[IDX(1,0,k)]);
      r = cabs2f(A[IDX(0,1,k)]) + cabs2f(A[IDX(1,1,k)]);
      P[k] = p;
      Q[2*k+0] = Q[2*k+1] = q;
      R[2*k+0] = R[2*k+1] = r;
    } /* k */
#endif
    /* Compute a sine/consine pair, part 2. The rotation angle is selected in
     * such a way that the L2 norm of the second column does not exceed the
     * norm of the first one. */
#if 1
    P_r = (xb_vecN_2xf32*)P;
    Q_r = (xb_vecN_2xf32*)Q;
    R_r = (xb_vecN_2xf32*)R;
    F_w = (xb_vecN_2xf32*)F;
    G_w = (xb_vecN_2xf32*)G;
    H_w = (xb_vecN_2xf32*)H;
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_LVN_2XF32_IP(_p, P_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_q, Q_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_r, R_r, 2*BBE_SIMD_WIDTH);
      /* if (q>=r) {
       *   f = crdivf(p,q); g = 1-r/q;
       * } else {
       *   f = crdivf(p,r); g = 1-q/r;
       * } */
      bge = BBE_OGEN_2XF32(_q,_r);
      _a = BBE_MOVN_2XF32T(_r,_q,bge);
      _b = BBE_MOVN_2XF32T(_q,_r,bge);
      _c = BBE_RECIP0N_2XF32(_b);
      _d = BBE_CONSTN_2XF32(1);
      BBE_MULSN_2XF32(_d,_b,_c);
      BBE_MULAN_2XF32(_c,_d,_c);
      _f = BBE_MULN_2XF32(_c,_p);
      BBE_MULSN_2XF32(_p,_b,_f);
      BBE_MULAN_2XF32(_f,_c,_p);
      _g = BBE_MULN_2XF32(_a,_c);
      BBE_MULSN_2XF32(_a,_b,_g);
      BBE_MULAN_2XF32(_g,_a,_c);
      _g = BBE_SUBN_2XF32(BBE_CONSTN_2XF32(1), _g);
      /* h = sqrtf(4*cabs2f(f)+g*g); */
      _t = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), _f);
      _t = BBE_MULN_2XF32(_t,_t);
      _t = BBE_ADDN_2XF32(_t, BBE_SHFLN_2XF32I(_t, BBE_SHFLI_SWAP_2));
      BBE_MULAN_2XF32(_t,_g,_g);
      _h = BBE_FSQRTN_2XF32(_t);
      BBE_SVN_2XF32_IP(_f, F_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_g, G_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_h, H_w, 2*BBE_SIMD_WIDTH);
    } /* k */
    __Pragma("no_reorder");
#else
    for ( k=0; k<L; k++ ) {
      p = P[k]; q = Q[2*k+0]; r = R[2*k+0];
      if (q>=r) {
        f = crdivf(p,q); g = 1-r/q;
      } else {
        f = crdivf(p,r); g = 1-q/r;
      }
      h = sqrtf(4*cabs2f(f)+g*g);
      F[k] = f; 
      G[2*k+0] = G[2*k+1] = g;
      H[2*k+0] = H[2*k+1] = h;
    } /* k */
#endif
    /* Compute a sine/consine pair, part 3. */
#if 1
    G_w = (xb_vecN_2xf32*)G;
    G_r = (xb_vecN_2xf32*)G;
    H_r = (xb_vecN_2xf32*)H;
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_LVN_2XF32_IP(_g, G_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_h, H_r, 2*BBE_SIMD_WIDTH);
      /* g = sqrtf((1+g/h)/2); */
      _a = BBE_RECIP0N_2XF32(_h);
      _b = BBE_CONSTN_2XF32(1);
      BBE_MULSN_2XF32(_b,_a,_h);
      BBE_MULAN_2XF32(_a,_b,_a);
      _t = BBE_MULN_2XF32(_g,_a);
      BBE_MULSN_2XF32(_g,_h,_t);
      BBE_MULAN_2XF32(_t,_a,_g);
      _t = BBE_ADDN_2XF32(BBE_CONSTN_2XF32(1),_t);
      _t = BBE_MULN_2XF32(BBE_CONSTN_2XF32(3),_t);
      _g = BBE_FSQRTN_2XF32(_t);
      BBE_SVN_2XF32_IP(_g, G_w, 2*BBE_SIMD_WIDTH);
    } /* k */
    __Pragma("no_reorder");
#else
    for ( k=0; k<L; k++ ) {
      g = G[2*k+0]; h = H[2*k+0];
      g = sqrtf((1+g/h)/2);
      G[2*k+0] = G[2*k+1] = g;
    } /* k */
#endif
    /* Compute a sine/consine pair, part 4. */
#if 1
    P_r  = (xb_vecN_2xf32*)P;  Q_r  = (xb_vecN_2xf32*)Q;
    R_r  = (xb_vecN_2xf32*)R;  F_r  = (xb_vecN_2xf32*)F;
    G_r  = (xb_vecN_2xf32*)G;  H_r  = (xb_vecN_2xf32*)H;
    SN_w = (xb_vecN_2xf32*)SN; CS_w = (xb_vecN_2xf32*)CS;
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_LVN_2XF32_IP(_p, P_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_q, Q_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_r, R_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_f, F_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_g, G_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_h, H_r, 2*BBE_SIMD_WIDTH);
      /* cs = _makecomplexf(g,0.f); */
      _cs = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0), _g, BBE_SELI_INTERLEAVE_2_EVEN);
      /* sn = crdivf(f,-h*g); */
      _a = BBE_MULN_2XF32(_h, BBE_NEGN_2XF32(_g));
      _b = BBE_RECIP0N_2XF32(_a);
      _c = BBE_CONSTN_2XF32(1);
      BBE_MULSN_2XF32(_c,_a,_b);
      BBE_MULAN_2XF32(_b,_c,_b);
      _sn = BBE_MULN_2XF32(_f, _b);
      BBE_MULSN_2XF32(_f,_a,_sn);
      BBE_MULAN_2XF32(_sn,_f,_b);
      /* if (q<r) {
       *   f = cs; cs = sn; sn = f;
       * } else if (cmagf(p)<=tol*r) {
       *   cs = c1f; sn = c0f;
       * } */
      blt = BBE_OLTN_2XF32(_q,_r);
      _t = BBE_ABSN_2XF32(_p);
      _t = BBE_ADDN_2XF32(_t, BBE_SHFLN_2XF32I(_t, BBE_SHFLI_SWAP_2));
      ble = BBE_OLEN_2XF32T(_t, BBE_MULN_2XF32(_tol,_r), BBE_NOTBN_2(blt));
      _f = _cs; _cs = BBE_MOVN_2XF32T(_sn,_cs,blt); _sn = BBE_MOVN_2XF32T(_f,_sn,blt);
      _cs = BBE_MOVN_2XF32T(_c1f,_cs,ble); _sn = BBE_MOVN_2XF32T(_c0f,_sn,ble);
      BBE_SVN_2XF32_IP(_sn, SN_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_cs, CS_w, 2*BBE_SIMD_WIDTH);
    } /* k */
    __Pragma("no_reorder");
#else
    for ( k=0; k<L; k++ ) {
      p = P[k]; q = Q[2*k+0]; r = R[2*k+0];
      f = F[k]; g = G[2*k+0]; h = H[2*k+0];
      cs = _makecomplexf(g,0.f);
      sn = crdivf(f,-h*g);
      if (q<r) {
        f = cs; cs = sn; sn = f;
      } else if (cmagf(p)<=tol*r) {
        cs = c1f; sn = c0f;
      }
      SN[k] = sn; CS[k] = cs;
    } /* k */
#endif
    /* Apply the rotator to working matrix: A <- A*[cs,sn;-sn',cs'] */
#if 1
    SN_r = (xb_vecN_2xf32*)SN;
    CS_r = (xb_vecN_2xf32*)CS;
    A_r  = (xb_vecN_2xf32*)&A[IDX(0,0,0)];
    A_w0 = (xb_vecN_2xf32*)&A[IDX(0,0,0)];
    A_w1 = (xb_vecN_2xf32*)&A[IDX(0,1,0)];
    A_w2 = (xb_vecN_2xf32*)&A[IDX(1,0,0)];
    A_w3 = (xb_vecN_2xf32*)&A[IDX(1,1,0)];
    for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
      BBE_LVN_2XF32_IP(_sn, SN_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(_cs, CS_r, 2*BBE_SIMD_WIDTH);
      /* a = A[IDX(0,0,k)]; b = A[IDX(0,1,k)];
       * c = A[IDX(1,0,k)]; d = A[IDX(1,1,k)]; */
      BBE_LVN_2XF32_XP(_a, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_c, A_r, L*sz_f32c);
      BBE_LVN_2XF32_XP(_d, A_r, 2*BBE_SIMD_WIDTH-3*L*sz_f32c);
      /* p = csubf(cmulf(a,cs), cmuljf(b,sn)); */
      _p = BBE_MULMN_2XF32(_a,_cs,0,8);
      BBE_MULMASN_2XF32(_p,_a,_cs,1,7);
      BBE_MULMASN_2XF32(_p,_b,_sn,3,8);
      BBE_MULMASN_2XF32(_p,_b,_sn,1,7);
      /* q = caddf(cmulf(a,sn), cmuljf(b,cs)); */
      _q = BBE_MULMN_2XF32(_a,_sn,0,8);
      BBE_MULMASN_2XF32(_q,_b,_cs,0,8);
      BBE_MULMASN_2XF32(_q,_a,_sn,1,7);
      BBE_MULMASN_2XF32(_q,_b,_cs,2,7);
      /* r = csubf(cmulf(c,cs), cmuljf(d,sn)); */
      _r = BBE_MULMN_2XF32(_c,_cs,0,8);
      BBE_MULMASN_2XF32(_r,_c,_cs,1,7);
      BBE_MULMASN_2XF32(_r,_d,_sn,3,8);
      BBE_MULMASN_2XF32(_r,_d,_sn,1,7);
      /* t = caddf(cmulf(c,sn), cmuljf(d,cs)); */
      _t = BBE_MULMN_2XF32(_c,_sn,0,8);
      BBE_MULMASN_2XF32(_t,_c,_sn,1,7);
      BBE_MULMASN_2XF32(_t,_d,_cs,0,8);
      BBE_MULMASN_2XF32(_t,_d,_cs,2,7);
      /* A[IDX(0,0,k)] = p; A[IDX(0,1,k)] = q;
       * A[IDX(1,0,k)] = r; A[IDX(1,1,k)] = t; */
      BBE_SVN_2XF32_IP(_p, A_w0, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_q, A_w1, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_r, A_w2, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(_t, A_w3, 2*BBE_SIMD_WIDTH);
    } /* k */
#else
    for ( k=0; k<L; k++ ) {
      sn = SN[k]; cs = CS[k];
      a = A[IDX(0,0,k)]; b = A[IDX(0,1,k)];
      c = A[IDX(1,0,k)]; d = A[IDX(1,1,k)];
      A[IDX(0,0,k)] = csubf(cmulf(a,cs), cmuljf(b,sn));
      A[IDX(0,1,k)] = caddf(cmulf(a,sn), cmuljf(b,cs));
      A[IDX(1,0,k)] = csubf(cmulf(c,cs), cmuljf(d,sn));
      A[IDX(1,1,k)] = caddf(cmulf(c,sn), cmuljf(d,cs));
    } /* k */
#endif
    /* Accumulate transformations: V <- V*[cs,sn;-sn',cs'] */
    if (V) {
#if 1
      SN_r = (xb_vecN_2xf32*)SN;
      CS_r = (xb_vecN_2xf32*)CS;
      V_r  = (xb_vecN_2xf32*)&V[IDX(0,0,0)];
      V_w0 = (xb_vecN_2xf32*)&V[IDX(0,0,0)];
      V_w1 = (xb_vecN_2xf32*)&V[IDX(0,1,0)];
      V_w2 = (xb_vecN_2xf32*)&V[IDX(1,0,0)];
      V_w3 = (xb_vecN_2xf32*)&V[IDX(1,1,0)];
      for ( k=0; k<L/(BBE_SIMD_WIDTH/4); k++ ) {
        BBE_LVN_2XF32_IP(_sn, SN_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(_cs, CS_r, 2*BBE_SIMD_WIDTH);
        /* a = V[IDX(0,0,k)]; b = V[IDX(0,1,k)];
         * c = V[IDX(1,0,k)]; d = V[IDX(1,1,k)]; */
        BBE_LVN_2XF32_XP(_a, V_r, L*sz_f32c);
        BBE_LVN_2XF32_XP(_b, V_r, L*sz_f32c);
        BBE_LVN_2XF32_XP(_c, V_r, L*sz_f32c);
        BBE_LVN_2XF32_XP(_d, V_r, 2*BBE_SIMD_WIDTH-3*L*sz_f32c);
        /* p = csubf(cmulf(a,cs), cmuljf(b,sn)); */
        _p = BBE_MULMN_2XF32(_a,_cs,0,8);
        BBE_MULMASN_2XF32(_p,_a,_cs,1,7);
        BBE_MULMASN_2XF32(_p,_b,_sn,3,8);
        BBE_MULMASN_2XF32(_p,_b,_sn,1,7);
        /* q = caddf(cmulf(a,sn), cmuljf(b,cs)); */
        _q = BBE_MULMN_2XF32(_a,_sn,0,8);
        BBE_MULMASN_2XF32(_q,_b,_cs,0,8);
        BBE_MULMASN_2XF32(_q,_a,_sn,1,7);
        BBE_MULMASN_2XF32(_q,_b,_cs,2,7);
        /* r = csubf(cmulf(c,cs), cmuljf(d,sn)); */
        _r = BBE_MULMN_2XF32(_c,_cs,0,8);
        BBE_MULMASN_2XF32(_r,_c,_cs,1,7);
        BBE_MULMASN_2XF32(_r,_d,_sn,3,8);
        BBE_MULMASN_2XF32(_r,_d,_sn,1,7);
        /* t = caddf(cmulf(c,sn), cmuljf(d,cs)); */
        _t = BBE_MULMN_2XF32(_c,_sn,0,8);
        BBE_MULMASN_2XF32(_t,_c,_sn,1,7);
        BBE_MULMASN_2XF32(_t,_d,_cs,0,8);
        BBE_MULMASN_2XF32(_t,_d,_cs,2,7);
        /* V[IDX(0,0,k)] = p; V[IDX(0,1,k)] = q;
         * V[IDX(1,0,k)] = r; V[IDX(1,1,k)] = t; */
        BBE_SVN_2XF32_IP(_p, V_w0, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(_q, V_w1, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(_r, V_w2, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(_t, V_w3, 2*BBE_SIMD_WIDTH);
      } /* k */
#else
      for ( k=0; k<L; k++ ) {
        sn = SN[k]; cs = CS[k];
        a = V[IDX(0,0,k)]; b = V[IDX(0,1,k)];
        c = V[IDX(1,0,k)]; d = V[IDX(1,1,k)];
        V[IDX(0,0,k)] = csubf(cmulf(a,cs), cmuljf(b,sn));
        V[IDX(0,1,k)] = caddf(cmulf(a,sn), cmuljf(b,cs));
        V[IDX(1,0,k)] = csubf(cmulf(c,cs), cmuljf(d,sn));
        V[IDX(1,1,k)] = caddf(cmulf(c,sn), cmuljf(d,cs));
      } /* k */
#endif
    } /* V */
  } /* n */
  /* Compute singular values. */
#if 1
  A_r0 = (xb_vecN_2xf32*)&A[IDX(0,0,0)];
  A_r1 = (xb_vecN_2xf32*)&A[IDX(1,0,0)];
  S_w = (xb_vecN_2xf32*)s;
  for ( k=0; k<2*L/(BBE_SIMD_WIDTH/2); k++ ) {
    BBE_LVN_2XF32_IP(_a, A_r0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(_b, A_r0, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(_c, A_r1, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(_d, A_r1, 2*BBE_SIMD_WIDTH);
    _p = BBE_MULN_2XF32(_a,_a); BBE_MULAN_2XF32(_p,_c,_c);
    _q = BBE_MULN_2XF32(_b,_b); BBE_MULAN_2XF32(_q,_d,_d);
    _p = BBE_ADDN_2XF32(_p, BBE_SHFLN_2XF32I(_p, BBE_SHFLI_SWAP_2));
    _q = BBE_ADDN_2XF32(_q, BBE_SHFLN_2XF32I(_q, BBE_SHFLI_SWAP_2));
    _g = BBE_SELN_2XF32I(_q, _p, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
    _h = BBE_FSQRTN_2XF32(_g);
    blt = BBE_OLTN_2XF32(_g,_rmin);
    BBE_CONSTN_2XF32T(_h,0,blt);
    BBE_SVN_2XF32_IP(_h, S_w, 2*BBE_SIMD_WIDTH);
  } /* k */
  __Pragma("no_reorder");
#else
  for ( k=0; k<2*L; k++ ) {
    g = cabs2f(A[IDX(0,0,k)]) + cabs2f(A[IDX(1,0,k)]);
    /* A singular value may be zero or too small to be evaluated to working precision */
    s[IDX(0,0,k)] = ( g<rmin ? 0.f : sqrtf(g) );
  } /* k */
#endif
  if (U) {
    /* Compute first left-singular vectors. */
#if 1
    S_r = (xb_vecN_2xf32*)&s[IDX(0,0,0)];
    A_r = (xb_vecN_2xf32*)&A[IDX(0,0,0)];
    U_w = (xb_vecN_2xf32*)&U[IDX(0,0,0)];
    for ( k=0; k<L/(BBE_SIMD_WIDTH/2); k++ ) {
      BBE_LVN_2XF32_IP(_h, S_r, 2*BBE_SIMD_WIDTH);
      _f = BBE_RECIP0N_2XF32(_h); _g = BBE_CONSTN_2XF32(1);
      BBE_MULSN_2XF32(_g,_h,_f);
      BBE_MULAN_2XF32(_f,_g,_f);
      _h0 = BBE_SHFLN_2XF32I(_h, BBE_SHFLI_DOUBLE_2_LO);
      _h1 = BBE_SHFLN_2XF32I(_h, BBE_SHFLI_DOUBLE_2_HI);
      _f0 = BBE_SHFLN_2XF32I(_f, BBE_SHFLI_DOUBLE_2_LO);
      _f1 = BBE_SHFLN_2XF32I(_f, BBE_SHFLI_DOUBLE_2_HI);
      /* c = ( h==0.f ? c1f : crdivf(A[IDX(0,0,k)],h) );
       * d = ( h==0.f ? c0f : crdivf(A[IDX(1,0,k)],h) ); */
      BBE_LVN_2XF32_XP(_a, A_r, 2*L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, 2*BBE_SIMD_WIDTH-2*L*sz_f32c);
      _c = BBE_MULN_2XF32(_a,_f0);
      BBE_MULSN_2XF32(_a,_c,_h0);
      BBE_MULAN_2XF32(_c,_a,_f0);
      _d = BBE_MULN_2XF32(_b,_f0); 
      BBE_MULSN_2XF32(_b,_d,_h0);
      BBE_MULAN_2XF32(_d,_b,_f0);
      beq = BBE_OEQN_2XF32(_h0, BBE_CONSTN_2XF32(0));
      _c = BBE_MOVN_2XF32T(_c1f,_c,beq);
      _d = BBE_MOVN_2XF32T(_c0f,_d,beq);
      _r = BBE_NEGN_2XF32(BBE_CONJN_2XF32(_d));
      _t = BBE_CONJN_2XF32(_c);
      BBE_SVN_2XF32_XP(_c, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_r, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_d, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_t, U_w, 2*BBE_SIMD_WIDTH-3*L*sz_f32c);
      /* c = ( h==0.f ? c1f : crdivf(A[IDX(0,0,k)],h) );
       * d = ( h==0.f ? c0f : crdivf(A[IDX(1,0,k)],h) ); */
      BBE_LVN_2XF32_XP(_a, A_r, 2*L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, 2*BBE_SIMD_WIDTH-2*L*sz_f32c);
      _c = BBE_MULN_2XF32(_a,_f1);
      BBE_MULSN_2XF32(_a,_c,_h1);
      BBE_MULAN_2XF32(_c,_a,_f1);
      _d = BBE_MULN_2XF32(_b,_f1); 
      BBE_MULSN_2XF32(_b,_d,_h1);
      BBE_MULAN_2XF32(_d,_b,_f1);
      beq = BBE_OEQN_2XF32(_h1, BBE_CONSTN_2XF32(0));
      _c = BBE_MOVN_2XF32T(_c1f,_c,beq);
      _d = BBE_MOVN_2XF32T(_c0f,_d,beq);
      _r = BBE_NEGN_2XF32(BBE_CONJN_2XF32(_d));
      _t = BBE_CONJN_2XF32(_c);
      BBE_SVN_2XF32_XP(_c, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_r, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_d, U_w, L*sz_f32c);
      BBE_SVN_2XF32_XP(_t, U_w, 2*BBE_SIMD_WIDTH-3*L*sz_f32c);
    } /* k */
    __Pragma("no_reorder");
#else
    for ( k=0; k<L; k++ ) {
      h = s[IDX(0,0,k)];
      U[IDX(0,0,k)] = ( h==0.f ? c1f : crdivf(A[IDX(0,0,k)],h) );
      U[IDX(1,0,k)] = ( h==0.f ? c0f : crdivf(A[IDX(1,0,k)],h) );
    } /* k */
#endif
    /* Compute second left-singular vectors. */
#if 1
    S_r = (xb_vecN_2xf32*)&s[IDX(0,1,0)];
    A_r = (xb_vecN_2xf32*)&A[IDX(0,1,0)];
    U_w = (xb_vecN_2xf32*)&U[IDX(0,1,0)];
    for ( k=0; k<L/(BBE_SIMD_WIDTH/2); k++ ) {
      BBE_LVN_2XF32_IP(_h, S_r, 2*BBE_SIMD_WIDTH);
      _f = BBE_RECIP0N_2XF32(_h); _g = BBE_CONSTN_2XF32(1);
      BBE_MULSN_2XF32(_g,_h,_f);
      BBE_MULAN_2XF32(_f,_g,_f);
      _h0 = BBE_SHFLN_2XF32I(_h, BBE_SHFLI_DOUBLE_2_LO);
      _h1 = BBE_SHFLN_2XF32I(_h, BBE_SHFLI_DOUBLE_2_HI);
      _f0 = BBE_SHFLN_2XF32I(_f, BBE_SHFLI_DOUBLE_2_LO);
      _f1 = BBE_SHFLN_2XF32I(_f, BBE_SHFLI_DOUBLE_2_HI);
      /* c = ( h==0.f ? cnegf(_conjf(U[IDX(1,0,k)])) : crdivf(A[IDX(0,1,k)],h) );
       * d = ( h==0.f ?       _conjf(U[IDX(0,0,k)])  : crdivf(A[IDX(1,1,k)],h) ); */
      BBE_LVN_2XF32_XP(_a, A_r, 2*L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, 2*BBE_SIMD_WIDTH-2*L*sz_f32c);
      _c = BBE_MULN_2XF32(_a,_f0);
      BBE_MULSN_2XF32(_a,_c,_h0);
      BBE_MULAN_2XF32(_c,_a,_f0);
      _d = BBE_MULN_2XF32(_b,_f0); 
      BBE_MULSN_2XF32(_b,_d,_h0);
      BBE_MULAN_2XF32(_d,_b,_f0);
      bneq = BBE_UNEQN_2XF32(_h0, BBE_CONSTN_2XF32(0));
      BBE_SVN_2XF32T_XP(_c, U_w, 2*L*sz_f32c, bneq);
      BBE_SVN_2XF32T_XP(_d, U_w, 2*BBE_SIMD_WIDTH-2*L*sz_f32c, bneq);
      /* c = ( h==0.f ? cnegf(_conjf(U[IDX(1,0,k)])) : crdivf(A[IDX(0,1,k)],h) );
       * d = ( h==0.f ?       _conjf(U[IDX(0,0,k)])  : crdivf(A[IDX(1,1,k)],h) ); */
      BBE_LVN_2XF32_XP(_a, A_r, 2*L*sz_f32c);
      BBE_LVN_2XF32_XP(_b, A_r, 2*BBE_SIMD_WIDTH-2*L*sz_f32c);
      _c = BBE_MULN_2XF32(_a,_f1);
      BBE_MULSN_2XF32(_a,_c,_h1);
      BBE_MULAN_2XF32(_c,_a,_f1);
      _d = BBE_MULN_2XF32(_b,_f1); 
      BBE_MULSN_2XF32(_b,_d,_h1);
      BBE_MULAN_2XF32(_d,_b,_f1);
      bneq = BBE_UNEQN_2XF32(_h1, BBE_CONSTN_2XF32(0));
      BBE_SVN_2XF32T_XP(_c, U_w, 2*L*sz_f32c, bneq);
      BBE_SVN_2XF32T_XP(_d, U_w, 2*BBE_SIMD_WIDTH-2*L*sz_f32c, bneq);
    } /* k */
#else
    for ( k=0; k<L; k++ ) {
      h = s[IDX(0,1,k)];
      U[IDX(0,1,k)] = ( h==0.f ? cnegf(_conjf(U[IDX(1,0,k)])) : crdivf(A[IDX(0,1,k)],h) );
      U[IDX(1,1,k)] = ( h==0.f ?       _conjf(U[IDX(0,0,k)])  : crdivf(A[IDX(1,1,k)],h) );
    } /* k */
#endif
  }

  #undef IDX

} /* svd2x2sf() */

size_t svd2x2sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 2==M && 2==N ); 
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (   L*sz_f32c +  /* P: dot products               */
           2*L*sz_f32  +  /* Q: squared norm of 1st column */
           2*L*sz_f32  +  /* R: squared norm of 2nd column */
             L*sz_f32c +  /* F: normalized dot products    */
           2*L*sz_f32 );  /* G: intermediate data          */
}

#else /* HAVE_VFPU */

DISCARD_FUN( void, svd2x2sf, ( void * pScr,
                      complex_float * restrict U,
                      float32_t     * restrict s,
                      complex_float * restrict V,
                      complex_float * restrict A,
                      int L ) )

size_t svd2x2sf_getScratchSize ( int M, int N, int L )
{
  NASSERT( 2==M && 2==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (0);
}

#endif /* HAVE_VFPU */

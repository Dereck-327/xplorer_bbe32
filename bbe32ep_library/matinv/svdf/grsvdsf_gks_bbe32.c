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
    Internal functions of complex Golub-Kahan SVD step implementation for 
    generic matrix sizes.
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#if 0
#include <math.h>
#endif

/* Portable data types. */
#include "NatureDSP_types.h"
/* Complex math functions. */
#include "NatureDSP_Baseband_complex.h"
/* Common utility declarations. */
#include "common.h"
/* SVD common declarations */
#include "svd_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define MIN(a,b)     ((a)<(b)?(a):(b))
#define MAX(a,b)     ((a)>(b)?(a):(b))
#define SIDX(i,k)    ((i)*L+(k))

#define sz_vbn4      sizeof(vboolN_4)
#if 0
/* Set p-th boolean to b
 * Input:
 *   p                         Boolean index, 0<=p<L
 *   b                         Boolean value, zero or non-zero
 * Input/Output:
 *   vb[L/(BBE_SIMD_WIDTH/4)]  Array of vector booleans */
static void vbn4Set( vboolN_4 * vb, int p, uint8_t b )
{
  vboolN vb0,vb1,vbm;
  int m,n;
  m = p/(BBE_SIMD_WIDTH/4); n = p%(BBE_SIMD_WIDTH/4);
  vbm = BBE_EQNX16(BBE_SRAINX16(BBE_SEQNX16(),2), BBE_MOVVA16(n));
  vb0 = BBE_ANDNOTBN(BBE_MOVN_FROMN_4(vb[m]),vbm);
  vb1 = BBE_ANDBN(BBE_NEQNX16(BBE_MOVVA16(b),BBE_MOVVI16(0)),vbm);
  vb[m] = BBE_MOVN_4_FROMN(BBE_ORBN(vb0,vb1));
}

/* Get value of p-th boolean.
 * Input:
 *   p                         Boolean index, 0<=p<L
 *   vb[L/(BBE_SIMD_WIDTH/4)]  Array of vector booleans
 * Return value:
 *   Boolean value, 0 or 1. */
static uint8_t vbn4Get( const vboolN_4 * vb, int p )
{
  xb_vecNx16 v;
  int m,n;
  m = p/(BBE_SIMD_WIDTH/4); n = p%(BBE_SIMD_WIDTH/4);
  v = BBE_MOVNX16T(BBE_MOVVI16(1), BBE_MOVVI16(0), BBE_MOVN_FROMN_4(vb[m]));
  return (uint8_t)BBE_EXTRANX16(v,4*n);
}

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
 * Internal function of Golub-Kahan SVD step implementation.
 * Summarize index data and fetch data for Wilkinson's shift computation.
 * Input:
 *   itsLim       Iterations count limit
 *   N            Number of columns in an input matrix
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_l[L]       Left index of working subblocks
 *   a_k[L]       Right index of working subblocks
 * Input/Output:
 *   a_its[L]     Iteration counters
 * Output:
 *   p_l_lo       Lowest left index over input matrices
 *   p_k_up       Uppermost right index over input matrices
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 *   a_a[L]       Fetched data: D[a_k[0..L-1]-1] 
 *   a_b[L]       Fetched data: D[a_k[0..L-1]] 
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1] 
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2] 
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_its,rng_l,rng_k,a_m,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must not overlap
 *   a_its,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_fetch( 
                    int16_t       * restrict p_l_lo,
                    int16_t       * restrict p_k_up,
                    int16_t       * restrict a_its,
                    vboolN_4      * restrict a_m,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
                    complex_float * restrict a_x,
                    complex_float * restrict a_y,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const int16_t       *          a_l,
              const int16_t       *          a_k,
              const complex_float *          D,
              const complex_float *          F,
              int itsLim, int N, int L )
{
#if 1
  const xb_vecNx16    * restrict ITS_r;
  const long long     * restrict ITS_rs;
        xb_vecNx16    * restrict ITS_w;
        vboolN_4      * restrict M_w;
        xb_vecN_2xf32 * restrict A_w;
        xb_vecN_2xf32 * restrict B_w;
        xb_vecN_2xf32 * restrict X_w;
        xb_vecN_2xf32 * restrict Y_w;
        xb_vecN_2xf32 * restrict C_w;
        xb_vecN_2xf32 * restrict S_w;
  const xb_vecNx16    * restrict L_r;
  const long long     * restrict L_rs;
  const xb_vecNx16    * restrict K_r;
  const long long     * restrict K_rs;
  const xb_vecN_2xf32 * restrict D0_r;
  const xb_vecN_2xf32 * restrict D1_r;
  const xb_vecN_2xf32 * restrict F0_r;
  const xb_vecN_2xf32 * restrict F1_r;

  int n,p;
  int l_lo,l_up,k_lo,k_up;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  static const int16_t ALIGN(32) isel[BBE_SIMD_WIDTH] = {
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
  };

  NASSERT_ALIGN(a_its, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

#if 1
  {
    xb_vecNx16 vl_lo, vl_up, vk_lo, vk_up, vl, vk, its, vitsLim;
    vboolN blt;
    vl_lo = vk_lo = BBE_MOVVA16(N-1);
    vl_up = vk_up = BBE_ZERONX16();
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(its, ITS_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      blt = BBE_LTNX16(its, vitsLim);
      BBE_MINNX16T(vl_lo, vl_lo, vl, blt);
      BBE_MAXNX16T(vl_up, vl_up, vl, blt);
      BBE_MINNX16T(vk_lo, vk_lo, vk, blt);
      BBE_MAXNX16T(vk_up, vk_up, vk, blt);
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(its, ITS_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      its = BBE_SELNX16I(vitsLim, its, BBE_SELI_EXTRACT_LO_HALVES);
      blt = BBE_LTNX16(its, vitsLim);
      BBE_MINNX16T(vl_lo, vl_lo, vl, blt);
      BBE_MAXNX16T(vl_up, vl_up, vl, blt);
      BBE_MINNX16T(vk_lo, vk_lo, vk, blt);
      BBE_MAXNX16T(vk_up, vk_up, vk, blt);
    } /* L */
    vl_lo = BBE_MAXNX16(vl_lo, BBE_ZERONX16());
    vk_lo = BBE_MAXNX16(vk_lo, BBE_ZERONX16());
    l_lo = BBE_RMINA16VNX16(vl_lo); l_up = BBE_RMAXA16VNX16(vl_up);
    k_lo = BBE_RMINA16VNX16(vk_lo); k_up = BBE_RMAXA16VNX16(vk_up);
  }
#else
  l_lo = k_lo = N-1; l_up = k_up = 0;
  for ( p=0; p<L; p++ ) {
    if (a_its[p]<itsLim) {
      if (l_lo>a_l[p]) l_lo = a_l[p];
      if (l_up<a_l[p]) l_up = a_l[p];
      if (k_lo>a_k[p]) k_lo = a_k[p];
      if (k_up<a_k[p]) k_up = a_k[p];
    }
  } /* p */
  if (l_lo<0) l_lo = 0;
  if (k_lo<0) k_lo = 0;
#endif
  if (l_lo<k_up) {
#if 1
    {
      xb_vecN_4x64 sl,sk,sits;
      xb_vecNx16 vl,vk,vn,vits;
      xb_vecN_2xf32 a,b,c,s;
      vboolN bgt,blt,beq;
      vboolN_2 bm,bn;
      vselN sel;
      sel = BBE_MOVVSELNX16(BBE_LVNX16_I((xb_vecNx16*)isel, 0), 0);
      D0_r = (xb_vecN_2xf32*)&D[SIDX(MAX(k_lo,1)-1,0)];
      D1_r = (xb_vecN_2xf32*)&D[SIDX(MAX(k_lo,1),0)];
      F0_r = (xb_vecN_2xf32*)&F[SIDX(MAX(k_lo,1)-1,0)];
      F1_r = (xb_vecN_2xf32*)&F[SIDX(MAX(k_lo,1)-2,0)];
      for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
        vn = BBE_MOVVA16(n);
        /* bn <- n>1 */
        bn = BBE_MOVN_2_FROMN(BBE_GTNX16(vn, BBE_MOVVINT16(1)));
        ITS_rs = (long long*)a_its;
        L_rs = (long long*)a_l;
        K_rs = (long long*)a_k;
        A_w = (xb_vecN_2xf32*)a_a;
        B_w = (xb_vecN_2xf32*)a_b;
        C_w = (xb_vecN_2xf32*)a_c;
        S_w = (xb_vecN_2xf32*)a_s;
        for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
          /* bgt <- a_k[p]>a_l[p] */
          BBE_LSN_4X64_IP(sl, L_rs, sizeof(*L_rs));
          BBE_LSN_4X64_IP(sk, K_rs, sizeof(*K_rs));
          vl = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sl), sel);
          vk = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sk), sel);
          bgt = BBE_GTNX16(vk,vl);
          /* blt <- a_its[p]<itsLim */
          BBE_LSN_4X64_IP(sits, ITS_rs, sizeof(*ITS_rs));
          vits = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sits), sel);
          blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim)); 
          /* beq <- a_k[p]==n */
          beq = BBE_EQNX16(vk,vn); 
          bm = BBE_MOVN_2_FROMN(BBE_ANDBN(bgt, BBE_ANDBN(blt,beq)));
          /* D[SIDX(n-1,p)] */
          BBE_LVN_2XF32_IP(a, D0_r, 2*BBE_SIMD_WIDTH);
          BBE_SVN_2XF32T_IP(a, A_w, 2*BBE_SIMD_WIDTH, bm);
          /* D[SIDX(n,p)] */
          BBE_LVN_2XF32_IP(b, D1_r, 2*BBE_SIMD_WIDTH);
          BBE_SVN_2XF32T_IP(b, B_w, 2*BBE_SIMD_WIDTH, bm);
          /* F[SIDX(n-1,p)] */
          BBE_LVN_2XF32_IP(c, F0_r, 2*BBE_SIMD_WIDTH);
          BBE_SVN_2XF32T_IP(c, C_w, 2*BBE_SIMD_WIDTH, bm);
          /* n>1 ? F[SIDX(n-2,p)] : c0f */
          BBE_LVN_2XF32T_IP(s, F1_r, 2*BBE_SIMD_WIDTH, bn);
          BBE_SVN_2XF32T_IP(s, S_w, 2*BBE_SIMD_WIDTH, bm);
        } /* p */
      } /* n */
    }
#else
    for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_k[p]==n) {
            a_a[p] = D[SIDX(n-1,p)]; 
            a_b[p] = D[SIDX(n,p)]; 
            a_c[p] = F[SIDX(n-1,p)]; 
            a_s[p] = ( n>1 ? F[SIDX(n-2,p)] : c0f );
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
#endif
#if 1
    {
      xb_vecN_4x64 sl,sk,sits;
      xb_vecNx16 vl,vk,vits;
      xb_vecN_2xf32 x,y;
      vboolN bgt,blt,beq;
      vboolN_2 bm;
      vselN sel;
      sel = BBE_MOVVSELNX16(BBE_LVNX16_I((xb_vecNx16*)isel, 0), 0);
      D0_r = (xb_vecN_2xf32*)&D[SIDX(l_lo,0)];
      F0_r = (xb_vecN_2xf32*)&F[SIDX(l_lo,0)];
      for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
        ITS_rs = (long long*)a_its;
        L_rs = (long long*)a_l;
        K_rs = (long long*)a_k;
        X_w = (xb_vecN_2xf32*)a_x;
        Y_w = (xb_vecN_2xf32*)a_y;
        for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
          /* bgt <- a_k[p]>a_l[p] */
          BBE_LSN_4X64_IP(sl, L_rs, sizeof(*L_rs));
          BBE_LSN_4X64_IP(sk, K_rs, sizeof(*K_rs));
          vl = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sl), sel);
          vk = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sk), sel);
          bgt = BBE_GTNX16(vk,vl);
          /* blt <- a_its[p]<itsLim */
          BBE_LSN_4X64_IP(sits, ITS_rs, sizeof(*ITS_rs));
          vits = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sits), sel);
          blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim)); 
          /* beq <- a_l[p]==n */
          beq = BBE_EQNX16(vl, BBE_MOVVA16(n)); 
          bm = BBE_MOVN_2_FROMN(BBE_ANDBN(bgt, BBE_ANDBN(blt,beq)));
          /* D[SIDX(n,p)] */
          BBE_LVN_2XF32_IP(x, D0_r, 2*BBE_SIMD_WIDTH);
          BBE_SVN_2XF32T_IP(x, X_w, 2*BBE_SIMD_WIDTH, bm);
          /* F[SIDX(n,p)] */
          BBE_LVN_2XF32_IP(y, F0_r, 2*BBE_SIMD_WIDTH);
          BBE_SVN_2XF32T_IP(y, Y_w, 2*BBE_SIMD_WIDTH, bm);
        } /* p */
      } /* n */
    }
#else
    for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_l[p]==n) {
            a_x[p] = D[SIDX(n,p)]; 
            a_y[p] = F[SIDX(n,p)]; 
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
#endif
  } /* l_lo, k_up */
#if 1
  M_w = a_m + l_lo*LW4;
  for ( n=l_lo; n<k_up; n++ ) {
    xb_vecNx16 vl,vk,vn,vits,vitsLim;
    vboolN bl,bk,bi;
    vboolN bm0,bm1,bm2,bm3;
    vn = BBE_MOVVA16(n);
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bl = BBE_LENX16(vl,vn);
      bk = BBE_LTNX16(vn,vk);
      bi = BBE_LTNX16(vits, vitsLim);
      BBE_EXTRACTB(bm2, bm0, BBE_ANDBN(bl, BBE_ANDBN(bk,bi)));
      BBE_EXTRACTB(bm1,bm0,bm0);
      BBE_EXTRACTB(bm3,bm2,bm2);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm0), M_w, sz_vbn4);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm1), M_w, sz_vbn4);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm2), M_w, sz_vbn4);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm3), M_w, sz_vbn4);
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bl = BBE_LENX16(vl,vn);
      bk = BBE_LTNX16(vn,vk);
      bi = BBE_LTNX16(vits, vitsLim);
      BBE_EXTRACTB(bm2, bm0, BBE_ANDBN(bl, BBE_ANDBN(bk,bi)));
      BBE_EXTRACTB(bm1,bm0,bm0);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm0), M_w, sz_vbn4);
      BBE_SBN_4_IP(BBE_MOVN_4_FROMN(bm1), M_w, sz_vbn4);
    } /* L */
  } /* n */
#else
  for ( n=l_lo; n<k_up; n++ ) {
    for ( p=0; p<L; p++ ) {
      vbn4Set(a_m+n*LW4, p, a_l[p]<=n && n<a_k[p] && a_its[p]<itsLim);
    } /* p */
  } /* n */
#endif
#if 1
  {
    xb_vecNx16 vl,vk,vits,vitsLim;
    vboolN bgt, blt;
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = ITS_w = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bgt = BBE_GTNX16(vk,vl);
      blt = BBE_LTNX16(vits,vitsLim);
      vits = BBE_ADDNX16(vits, BBE_MOVVINT16(1));
      BBE_SVNX16T_IP(vits, ITS_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN(bgt,blt));
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      vits = BBE_SELNX16I(vitsLim, vits, BBE_SELI_EXTRACT_LO_HALVES);
      bgt = BBE_GTNX16(vk,vl);
      blt = BBE_LTNX16(vits,vitsLim);
      vits = BBE_ADDNX16(vits, BBE_MOVVINT16(1));
      BBE_SVNX16T_IP(vits, ITS_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN(bgt,blt));
    } /* L */
  }
#else
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
      a_its[p]++;
    }
  } /* p */
#endif
  *p_l_lo = l_lo; *p_k_up = k_up;
#else
  int n,p;
  int l_lo,l_up,k_lo,k_up;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  const complex_float c0f = _makecomplexf(0.f,0.f);

  NASSERT_ALIGN(a_its, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  l_lo = k_lo = N-1; l_up = k_up = 0;
  for ( p=0; p<L; p++ ) {
    if (a_its[p]<itsLim) {
      if (l_lo>a_l[p]) l_lo = a_l[p];
      if (l_up<a_l[p]) l_up = a_l[p];
      if (k_lo>a_k[p]) k_lo = a_k[p];
      if (k_up<a_k[p]) k_up = a_k[p];
    }
  } /* p */
  if (l_lo<0) l_lo = 0;
  if (k_lo<0) k_lo = 0;
  if (l_lo<k_up) {
    for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_k[p]==n) {
            a_a[p] = D[SIDX(n-1,p)]; 
            a_b[p] = D[SIDX(n,p)]; 
            a_c[p] = F[SIDX(n-1,p)]; 
            a_s[p] = ( n>1 ? F[SIDX(n-2,p)] : c0f );
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
    for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_l[p]==n) {
            a_x[p] = D[SIDX(n,p)]; 
            a_y[p] = F[SIDX(n,p)]; 
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
  } /* l_lo, k_up */
  for ( n=l_lo; n<k_up; n++ ) {
    for ( p=0; p<L; p++ ) {
      vbn4Set(a_m+n*LW4, p, a_l[p]<=n && n<a_k[p] && a_its[p]<itsLim);
    } /* p */
  } /* n */
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
      a_its[p]++;
    }
  } /* p */
  *p_l_lo = l_lo; *p_k_up = k_up;
#endif
} /* grsvdsf_gks_fetch() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Wilkinson's shift.
 * Input:
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1].  May be reused
 *                as a temporal storage of intermediate results.
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2]. May be reused
 *                as a temporal storage of intermediate results.
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 * Input/Output:
 *   a_a[L]       In:  fetched data: D[a_k[0..L-1]-1] 
 *                Out: element (0,0) of B'*B-mu*I
 *   a_b[L]       In:  fetched data: D[a_k[0..L-1]] 
 *                Out: element (0,1) of B'*B-mu*I
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_x,a_y,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_wilkShift(
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const complex_float *          a_x,
              const complex_float *          a_y,
              int L )
{
  /*
   * MATLAB outline:
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
   */
#if 1
  const xb_vecN_2xf32 * restrict A_r;
        xb_vecN_2xf32 * restrict A_w;
  const xb_vecN_2xf32 * restrict B_r;
        xb_vecN_2xf32 * restrict B_w;
  const xb_vecN_2xf32 * restrict C_r;
  const xb_vecN_2xf32 * restrict S_r;
  const xb_vecN_2xf32 * restrict X_r;
  const xb_vecN_2xf32 * restrict Y_r;
  const xb_vecN_2xf32 * restrict T_r;
        xb_vecN_2xf32 * restrict T_w;
  const xb_vecN_2xf32 * restrict F_r;
  const xb_vecN_2xf32 * restrict G_r;

  float32_t * restrict a_t;
  float32_t * restrict a_f;
  float32_t * restrict a_g;
  int p;
  
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  a_t = (float32_t*)a_s;
  a_g = (float32_t*)a_s + L;
  a_f = (float32_t*)a_c;

#if 1
  {
    xb_vecN_2xf32 f0,f1,g0,g1,r0,r1,s0,s1,t;
    xb_vecN_2xf32 x0,x1,y0,y1,v0,v1,w0,w1;
    A_r = (xb_vecN_2xf32*)a_a;
    B_r = (xb_vecN_2xf32*)a_b;
    C_r = (xb_vecN_2xf32*)a_c;
    S_r = (xb_vecN_2xf32*)a_s;
    T_w = (xb_vecN_2xf32*)a_t;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      /* D[SIDX(a_k[p]-1,p)] */
      BBE_LVN_2XF32_IP(x0, A_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(x1, A_r, 2*BBE_SIMD_WIDTH);
      /* D[SIDX(a_k[p],p)] */
      BBE_LVN_2XF32_IP(y0, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(y1, B_r, 2*BBE_SIMD_WIDTH);
      /* F[SIDX(a_k[p]-1,p)] */
      BBE_LVN_2XF32_IP(v0, C_r, 2*BBE_SIMD_WIDTH); 
      BBE_LVN_2XF32_IP(v1, C_r, 2*BBE_SIMD_WIDTH); 
      /* F[SIDX(a_k[p]-2,p)] */
      BBE_LVN_2XF32_IP(w0, S_r, 2*BBE_SIMD_WIDTH); 
      BBE_LVN_2XF32_IP(w1, S_r, 2*BBE_SIMD_WIDTH); 
      /* t <- real((x-y)*conj(x+y)+(w-v)*conj(w+v)) */
      BBE_ADDSUBN_2XF32(f0,g0,x0,y0); BBE_ADDSUBN_2XF32(f1,g1,x1,y1); 
      BBE_DSELN_2XF32I(f1, f0, f1, f0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELN_2XF32I(g1, g0, g1, g0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_ADDSUBN_2XF32(r0,s0,w0,v0); BBE_ADDSUBN_2XF32(r1,s1,w1,v1);
      BBE_DSELN_2XF32I(r1, r0, r1, r0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELN_2XF32I(s1, s0, s1, s0, BBE_DSELI_DEINTERLEAVE_2);
      t = BBE_MULN_2XF32(f0,g0); BBE_MULAN_2XF32(t,f1,g1);
      BBE_MULAN_2XF32(t,r0,s0); BBE_MULAN_2XF32(t,r1,s1);
      BBE_SVN_2XF32_IP(t, T_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    complex_float v,w,x,y;
    for ( p=0; p<L; p++ ) {
      x = a_a[p]; /* D[SIDX(a_k[p]-1,p)] */
      y = a_b[p]; /* D[SIDX(a_k[p],p)] */
      v = a_c[p]; /* F[SIDX(a_k[p]-1,p)] */
      w = a_s[p]; /* F[SIDX(a_k[p]-2,p)] */
      a_t[p] = cmulj_rf(csubf(x,y), caddf(x,y)) + 
               cmulj_rf(csubf(w,v), caddf(w,v));
    } /* p */
  }
#endif
  /* F[SIDX(a_k[p]-1,p)] */
  vcabsf(a_g,a_c,L);
  /* D[SIDX(a_k[p]-1,p)] */
  vcabsf(a_f,a_a,L);
#if 1
  {
    xb_vecN_2xf32 f,g,r,t;
    F_r = (xb_vecN_2xf32*)a_f;
    G_r = (xb_vecN_2xf32*)a_g;
    T_r = T_w = (xb_vecN_2xf32*)a_t;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(f, F_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(g, G_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(t, T_r, 2*BBE_SIMD_WIDTH);
      /* t <- t/(2*f*g) */
      r = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), BBE_MULN_2XF32(f,g));
      BBE_SVN_2XF32_IP(BBE_DIVN_2XF32(t,r), T_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    float32_t f,g,t;
    for ( p=0; p<L; p++ ) {
      f = a_f[p]; g = a_g[p]; t = a_t[p];
      a_t[p] = t/(2*f*g);
    } /* p */
  }
#endif

#if 1
  {
    xb_vecN_2xf32 r,s,t;
    T_r = T_w = (xb_vecN_2xf32*)a_t;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(t, T_r, 2*BBE_SIMD_WIDTH);
      /* r <- sqrtf(1.f+t*t) */
      s = BBE_CONSTN_2XF32(1); BBE_MULAN_2XF32(s,t,t);
      r = BBE_FSQRTN_2XF32(s);
      /* t <- t<0 ? t-r : t+r */
      BBE_NEGN_2XF32T(r, r, BBE_OLTN_2XF32(t, BBE_CONSTN_2XF32(0)));
      BBE_SVN_2XF32_IP(BBE_ADDN_2XF32(t,r), T_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    float32_t r,t;
    for ( p=0; p<L; p++ ) {
      t = a_t[p];
      r = sqrtf(1.f+t*t); if (t<0) r = -r;
      a_t[p] = t+r;
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_2xf32 f,g,r,s,t;
    F_r = (xb_vecN_2xf32*)a_f;
    G_r = (xb_vecN_2xf32*)a_g;
    T_r = T_w = (xb_vecN_2xf32*)a_t;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(f, F_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(g, G_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(t, T_r, 2*BBE_SIMD_WIDTH);
      /* s <- f/t; magnitude of the divisor is large enough to use a 
       * simplified division sequence with no exponent manipulation. */
      r = BBE_RECIP0N_2XF32(t);
      s = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(s,t,r); BBE_MULAN_2XF32(r,s,r);
      s = BBE_MULN_2XF32(f,r); BBE_MULSN_2XF32(f,s,t); BBE_MULAN_2XF32(s,f,r);
      /* t <- g*(s-g) */
      t = BBE_MULN_2XF32(g, BBE_SUBN_2XF32(s,g));
      BBE_SVN_2XF32_IP(t, T_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    float32_t f,g,t;
    for ( p=0; p<L; p++ ) {
      f = a_f[p]; g = a_g[p]; t = a_t[p];
      /* Magnitude of the divisor is large enough to use a simplified
       * division sequence with no exponent manipulation. */
      NASSERT(!(fabsf(t)<1.f));
      a_t[p] = g*(f/t-g);
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_2xf32 a0,a1,b0,b1,r0,r1,s0,s1;
    xb_vecN_2xf32 u0,u1,y0,y1,z0,z1,t;
    A_w = (xb_vecN_2xf32*)a_a;
    B_r = B_w = (xb_vecN_2xf32*)a_b;
    X_r = (xb_vecN_2xf32*)a_x;
    Y_r = (xb_vecN_2xf32*)a_y;
    T_r = (xb_vecN_2xf32*)a_t;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      /* D[SIDX(a_k[p],p)] */
      BBE_LVN_2XF32_IP(y0, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(y1, B_r, 2*BBE_SIMD_WIDTH);
      /* D[SIDX(a_l[p],p)] */
      BBE_LVN_2XF32_IP(z0, X_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(z1, X_r, 2*BBE_SIMD_WIDTH);
      /* F[SIDX(a_l[p],p)] */
      BBE_LVN_2XF32_IP(u0, Y_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(u1, Y_r, 2*BBE_SIMD_WIDTH);
      /* a <- real((z-y)*conj(z+y))+t */
      BBE_ADDSUBN_2XF32(r0,s0,z0,y0); BBE_ADDSUBN_2XF32(r1,s1,z1,y1);
      BBE_DSELN_2XF32I(r1, r0, r1, r0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELN_2XF32I(s1, s0, s1, s0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_LVN_2XF32_IP(t, T_r, 2*BBE_SIMD_WIDTH);
      BBE_MULAN_2XF32(t,r0,s0); BBE_MULAN_2XF32(t,r1,s1);
      BBE_DSELN_2XF32I(a1, a0, BBE_CONSTN_2XF32(0), t, BBE_DSELI_INTERLEAVE_2);
      BBE_SVN_2XF32_IP(a0, A_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(a1, A_w, 2*BBE_SIMD_WIDTH);
      /* b <- u*conj(z) */
      b0 = BBE_MULMN_2XF32(u0,z0,0,8); BBE_MULMASN_2XF32(b0,u0,z0,2,7);
      b1 = BBE_MULMN_2XF32(u1,z1,0,8); BBE_MULMASN_2XF32(b1,u1,z1,2,7);
      BBE_SVN_2XF32_IP(b0, B_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(b1, B_w, 2*BBE_SIMD_WIDTH);
    } /* p */
  }
#else
  {
    complex_float u,y,z;
    float32_t t;
    for ( p=0; p<L; p++ ) {
      y = a_b[p]; /* D[SIDX(a_k[p],p)] */
      z = a_x[p]; /* D[SIDX(a_l[p],p)] */
      u = a_y[p]; /* F[SIDX(a_l[p],p)] */
      t = a_t[p];
      /* a <- t11-mu; b <- t12 */
      a_a[p] = _makecomplexf( cmulj_rf(csubf(z,y), caddf(z,y)) + t, 0.f );
      a_b[p] = cmuljf(u,z);
    } /* p */
  }
#endif
#else
  float32_t * restrict a_f;
  float32_t * restrict a_g;
  float32_t * restrict a_t;
  complex_float u,v,w,x,y,z;
  float32_t f,g,r,t;
  int p;
  
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  a_t = (float32_t*)a_s;
  for ( p=0; p<L; p++ ) {
    x = a_a[p]; /* D[SIDX(a_k[p]-1,p)] */
    y = a_b[p]; /* D[SIDX(a_k[p],p)] */
    v = a_c[p]; /* F[SIDX(a_k[p]-1,p)] */
    w = a_s[p]; /* F[SIDX(a_k[p]-2,p)] */
    a_t[p] = cmulj_rf(csubf(x,y), caddf(x,y)) + 
             cmulj_rf(csubf(w,v), caddf(w,v));
  } /* p */

  a_g = (float32_t*)a_s + L;
  /* F[SIDX(a_k[p]-1,p)] */
  vcabsf(a_g,a_c,L);

  a_f = (float32_t*)a_c;
  /* D[SIDX(a_k[p]-1,p)] */
  vcabsf(a_f,a_a,L);

  for ( p=0; p<L; p++ ) {
    f = a_f[p]; g = a_g[p]; t = a_t[p];
    a_t[p] = t/(2*f*g);
  } /* p */

  for ( p=0; p<L; p++ ) {
    t = a_t[p];
    r = sqrtf(1.f+t*t); if (t<0) r = -r;
    a_t[p] = t+r;
  } /* p */

  for ( p=0; p<L; p++ ) {
    f = a_f[p]; g = a_g[p]; t = a_t[p];
    /* Magnitude of the divisor is large enough to use a simplified
     * division sequence with no exponent manipulation. */
    NASSERT(!(fabsf(t)<1.f));
    a_t[p] = g*(f/t-g);
  } /* p */

  for ( p=0; p<L; p++ ) {
    y = a_b[p]; /* D[SIDX(a_k[p],p)] */
    z = a_x[p]; /* D[SIDX(a_l[p],p)] */
    u = a_y[p]; /* F[SIDX(a_l[p],p)] */
    t = a_t[p];
    /* a <- t11-mu; b <- t12 */
    a_a[p] = _makecomplexf( cmulj_rf(csubf(z,y), caddf(z,y)) + t, 0.f );
    a_b[p] = cmuljf(u,z);
  } /* p */
#endif
} /* grsvdsf_gks_wilkShift() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Givens's rotating matrix.
 * For real data:
 *   G(a,b) <- [c,s;-s,c]: [a,b]*[c,s;-s,c] == [*,0], c^2+s^2 == 1
 * For complex data:
 *   G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) == [*,0], c*conj(c)+s*conj(s) == 1
 * Temporary:
 *   a_x[L]       Scratch array of L entries
 * Input:
 *   L            Number of matrices
 *   a_a[L]       a-values
 *   a_b[L]       b-values
 * Output:
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *  Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_givens( 
                    complex_float * restrict a_x,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const complex_float *          a_a, 
              const complex_float *          a_b, 
              int L )
{
  /*
   * MATLAB reference code:
   *
   * function [c,s] = givens(a,b)
   * if b~=0
   *   if b*conj(b)>a*conj(a)
   *     t = cdiv(a,b); r = sqrt(1+t*conj(t));
   *     s = 1/r; c = -conj(t)/r;
   *   else
   *     t = cdiv(b,a); r = sqrt(1+t*conj(t));
   *     c = 1/r; s = -conj(t)/r;
   *   end
   * else
   *   c = 1; s = 0;
   * end;
   */
#if 1
  const xb_vecN_2xf32 * restrict X0_r;
  const xb_vecN_2xf32 * restrict X1_r;
  const xb_vecN_2xf32 * restrict C_r;
        xb_vecN_2xf32 * restrict C_w;
  const xb_vecN_2xf32 * restrict S_r;
        xb_vecN_2xf32 * restrict S_w;
  const xb_vecN_2xf32 * restrict A_r;
  const xb_vecN_2xf32 * restrict B_r;

  int p;

  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

#if 1
  {
    xb_vecN_2xf32 a,b,f,g,u,v,w;
    vboolN_2 blt,beq;

    const xb_vecN_2xf32 c1f = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0), 
                                              BBE_CONSTN_2XF32(1),
                                              BBE_SELI_INTERLEAVE_2_LO);
    A_r = (xb_vecN_2xf32*)a_a;
    B_r = (xb_vecN_2xf32*)a_b;
    C_w = (xb_vecN_2xf32*)a_c;
    S_w = (xb_vecN_2xf32*)a_s;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* u <- cabs2f(a); */
      BBE_LVN_2XF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
      u = BBE_MULN_2XF32(a,a);
      u = BBE_ADDN_2XF32(u, BBE_SHFLN_2XF32I(u, BBE_SHFLI_SWAP_2));
      /* v <- cabs2f(b); */
      BBE_LVN_2XF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
      v = BBE_MULN_2XF32(b,b);
      v = BBE_ADDN_2XF32(v, BBE_SHFLN_2XF32I(v, BBE_SHFLI_SWAP_2));
      /* f <- (u<v || isnan(u) || isnan(v)) ? a : b; */
      blt = BBE_ULTN_2XF32(u,v);
      f = BBE_MOVN_2XF32T(a,b,blt);
      BBE_SVN_2XF32_IP(f, C_w, 2*BBE_SIMD_WIDTH);
      /* g <- (u<v || isnan(u) || isnan(v)) ? b : a; */
      g = BBE_MOVN_2XF32T(b,a,blt);
      /* if (cmagf(b)==0) g <- c1f; */
      w = BBE_ABSN_2XF32(b);
      w = BBE_ADDN_2XF32(w, BBE_SHFLN_2XF32I(w, BBE_SHFLI_SWAP_2));
      beq = BBE_OEQN_2XF32(w, BBE_CONSTN_2XF32(0));
      g = BBE_MOVN_2XF32T(c1f,g,beq);
      BBE_SVN_2XF32_IP(g, S_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    complex_float a,b,f,g;
    float32_t u,v;
    const complex_float c1f = _makecomplexf(1.f,0.f);
    for ( p=0; p<L; p++ ) {
      a = a_a[p]; b = a_b[p];
      u = cabs2f(a); v = cabs2f(b);
      if (u<v || isnan(u) || isnan(v)) {
        f = a; g = b;
      } else {
        f = b; g = a;
      }
      if (cmagf(b)==0) {
        g = c1f;
      }
      a_c[p] = f; a_s[p] = g;
    } /* p */
  }
#endif

  vdividecxf(a_x,a_c,a_s,L);

#if 1
  {
    xb_vecN_2xf32 f0,f1,g0,g1,r,s,s0,s1,t0,t1;
    X0_r = X1_r = (xb_vecN_2xf32*)a_x;
    C_w = (xb_vecN_2xf32*)a_c;
    S_w = (xb_vecN_2xf32*)a_s;
    __Pragma("concurrent");
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      /* r <- 1/sqrtf(1.f+cabs2f(t)); */
      BBE_LVN_2XF32_IP(t0, X0_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(t1, X0_r, 2*BBE_SIMD_WIDTH);
      BBE_DSELN_2XF32I(s1, s0, t1, t0, BBE_DSELI_DEINTERLEAVE_2);
      s = BBE_CONSTN_2XF32(1);
      BBE_MULAN_2XF32(s,s0,s0); BBE_MULAN_2XF32(s,s1,s1);
      r = BBE_RSQRTN_2XF32(s);
      /* g <- _makecomplexf(r,0.f); */
      BBE_DSELN_2XF32I(g1, g0, BBE_CONSTN_2XF32(0), r, BBE_DSELI_INTERLEAVE_2);
      BBE_SVN_2XF32_IP(g0, S_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(g1, S_w, 2*BBE_SIMD_WIDTH);
      /* f <- _conjf(rcmulf(-r,t)); */
      BBE_LVN_2XF32_IP(t0, X1_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(t1, X1_r, 2*BBE_SIMD_WIDTH);
      f0 = BBE_MULMN_2XF32(g0,t0,1,4); f1 = BBE_MULMN_2XF32(g1,t1,1,4);
      BBE_SVN_2XF32_IP(f0, C_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(f1, C_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    complex_float f,g,t;
    float32_t r;
    for ( p=0; p<L; p++ ) {
      t = a_x[p];
      r = 1/sqrtf(1.f+cabs2f(t));
      f = _conjf(rcmulf(-r,t));
      g = _makecomplexf(r,0.f);
      a_c[p] = f; a_s[p] = g;
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_2xf32 a,b,f,g,u,v;
    vboolN_2 blt;
    A_r = (xb_vecN_2xf32*)a_a;
    B_r = (xb_vecN_2xf32*)a_b;
    C_r = C_w = (xb_vecN_2xf32*)a_c;
    S_r = S_w = (xb_vecN_2xf32*)a_s;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* u <- cabs2f(a); */
      BBE_LVN_2XF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
      u = BBE_MULN_2XF32(a,a);
      u = BBE_ADDN_2XF32(u, BBE_SHFLN_2XF32I(u, BBE_SHFLI_SWAP_2));
      /* v <- cabs2f(b); */
      BBE_LVN_2XF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
      v = BBE_MULN_2XF32(b,b);
      v = BBE_ADDN_2XF32(v, BBE_SHFLN_2XF32I(v, BBE_SHFLI_SWAP_2));
      /* c <- (u<v || isnan(u) || isnan(v)) ? f : g;
       * s <- (u<v || isnan(u) || isnan(v)) ? g : f */
      blt = BBE_ULTN_2XF32(u,v);
      BBE_LVN_2XF32_IP(f, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(g, S_r, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32F_IP(g, C_w, 2*BBE_SIMD_WIDTH, blt);
      BBE_SVN_2XF32F_IP(f, S_w, 2*BBE_SIMD_WIDTH, blt);
    } /* p */
  }
#else
  {
    complex_float a,b,f,g;
    float32_t u,v;
    for ( p=0; p<L; p++ ) {
      a = a_a[p]; b = a_b[p];
      f = a_c[p]; g = a_s[p];
      u = cabs2f(a); v = cabs2f(b);
      if (u>=v) {
        a_c[p] = g; a_s[p] = f; 
      }
    } /* p */
  }
#endif
#else
  complex_float a,b,f,g,t;
  float32_t r,u,v;
  int p;

  const complex_float c1f = _makecomplexf(1.f,0.f);

  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  for ( p=0; p<L; p++ ) {
    a = a_a[p]; b = a_b[p];
    u = cabs2f(a); v = cabs2f(b);
    if (u<v || isnan(u) || isnan(v)) {
      f = a; g = b;
    } else {
      f = b; g = a;
    }
    if (cmagf(b)==0) {
      g = c1f;
    }
    a_c[p] = f; a_s[p] = g;
  } /* p */

  vdividecxf(a_x,a_c,a_s,L);

  for ( p=0; p<L; p++ ) {
    t = a_x[p];
    r = 1/sqrtf(1.f+cabs2f(t));
    f = _conjf(rcmulf(-r,t));
    g = _makecomplexf(r,0.f);
    a_c[p] = f; a_s[p] = g;
  } /* p */

  for ( p=0; p<L; p++ ) {
    a = a_a[p]; b = a_b[p];
    f = a_c[p]; g = a_s[p];
    u = cabs2f(a); v = cabs2f(b);
    if (u>=v) {
      a_c[p] = g; a_s[p] = f; 
    }
  } /* p */
#endif
} /* grsvdsf_gks_givens() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 1st updating step of SVD step iteration: B <- B*G(a,b).
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_l[L]       Left index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n-1,n) Out: B(n,n)  
 *   a_b[L]       In: B(n-1,n+1) Out: B(n+1,n)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_l,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_l
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_step1(
                    complex_float * restrict D,
                    complex_float * restrict F,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const int16_t       *          a_l,
              const vboolN_4      *          a_m,
              int n, int L )
{
#if 1
        xb_vecN_4xcf32 * restrict Fm_w;
  const xb_vecN_4xcf32 * restrict Fn_r;
        xb_vecN_4xcf32 * restrict Fn_w;
  const xb_vecN_4xcf32 * restrict Dn_r;
  const xb_vecN_4xcf32 * restrict Dq_r;
        xb_vecN_4xcf32 * restrict Dq_w;
  const long long      * restrict L_rs;
  const xb_vecN_4xcf32 * restrict A_r;
        xb_vecN_4xcf32 * restrict A_w;
  const xb_vecN_4xcf32 * restrict B_r;
        xb_vecN_4xcf32 * restrict B_w;
  const xb_vecN_4xcf32 * restrict C_r;
  const xb_vecN_4xcf32 * restrict S_r;
  const vboolN_4       * restrict M_r;

  int p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  static const int16_t ALIGN(32) isel[BBE_SIMD_WIDTH] = {
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
  };

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
#if 1
  {
    xb_vecN_4xcf32 a,b,c,f,s;
    xb_vecNx16 vn,vl;
    xb_vecN_4x64 sl;
    vboolN_4 bm,bn;
    vselN sel = BBE_MOVVSELNX16(BBE_LVNX16_I((xb_vecNx16*)isel, 0), 0);
    vn = BBE_MOVVA16(n);
    Fm_w = (xb_vecN_4xcf32*)&F[SIDX(n-1,0)];
    A_r = (xb_vecN_4xcf32*)a_a;
    B_r = (xb_vecN_4xcf32*)a_b;
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    L_rs = (long long*)a_l;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      BBE_LSN_4X64_IP(sl, L_rs, sizeof(*L_rs));
      vl = BBE_SHFLNX16(BBE_MOVNX16_FROMN_4X64(sl), sel);
      bn = BBE_MOVN_4_FROMN(BBE_GTNX16(vn,vl));
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      /* if (n>a_l[p]) {
       *   F[SIDX(n-1,p)] <- csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p]));
       * }
       */
      BBE_LVN_4XCF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      f = BBE_MULN_4XCF32(c,a); BBE_MULSN_4XCF32(f,s,b);
      BBE_SVN_4XCF32T_IP(f, Fm_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN_4(bn,bm));
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        /* B <- B*G(a,b) */
        if (n>a_l[p]) {
          F[SIDX(n-1,p)] = csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p]));
        }
      } /* a_m */
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_4xcf32 a,b,c,s,x,y,z;
    vboolN_4 bm;
    Dn_r = (xb_vecN_4xcf32*)&D[SIDX(n,0)];
    Dq_r = (xb_vecN_4xcf32*)&D[SIDX(q,0)];
    Fn_r = (xb_vecN_4xcf32*)&F[SIDX(n,0)];
    A_w = (xb_vecN_4xcf32*)a_a;
    B_w = (xb_vecN_4xcf32*)a_b;
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* x <- D[SIDX(n,p)]; y <- F[SIDX(n,p)]; z <- D[SIDX(q,p)]; */
      BBE_LVN_4XCF32_IP(x, Dn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(y, Fn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(z, Dq_r, 2*BBE_SIMD_WIDTH);
      /* a <- csubf(cmulf(a_c[p],x), cmulf(a_s[p],y)); */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      a = BBE_MULN_4XCF32(c,x); BBE_MULSN_4XCF32(a,s,y);
      /* b <- cnegf(cmulf(a_s[p],z)); */
      b = BBE_MULMN_4XCF32(s,z,3,8);
      BBE_MULMASN_4XCF32(b,s,z,2,7);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(a, A_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(b, B_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    complex_float x,y,z;
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
        a_a[p] = csubf(cmulf(a_c[p],x), cmulf(a_s[p],y));
        a_b[p] = cnegf(cmulf(a_s[p],z));
      } /* a_m */
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_4xcf32 c,d,f,s,x,y,z;
    vboolN_4 bm;
    Dn_r = (xb_vecN_4xcf32*)&D[SIDX(n,0)];
    Dq_r = Dq_w = (xb_vecN_4xcf32*)&D[SIDX(q,0)];
    Fn_r = Fn_w = (xb_vecN_4xcf32*)&F[SIDX(n,0)];
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* x <- D[SIDX(n,p)]; y <- F[SIDX(n,p)]; z <- D[SIDX(q,p)]; */
      BBE_LVN_4XCF32_IP(x, Dn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(y, Fn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(z, Dq_r, 2*BBE_SIMD_WIDTH);
      /* F[SIDX(n,p)] <- caddf(cmuljf(x,a_s[p]), cmuljf(y,a_c[p])); */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      f = BBE_MULJN_4XCF32(x,s); BBE_MULJAN_4XCF32(f,y,c);
      /* D[SIDX(q,p)] <- cmuljf(z,a_c[p]); */
      d = BBE_MULJN_4XCF32(z,c);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(f, Fn_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(d, Dq_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
  }
#else
  {
    complex_float x,y,z;
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
        F[SIDX(n,p)] = caddf(cmuljf(x,a_s[p]), cmuljf(y,a_c[p])); 
        D[SIDX(q,p)] = cmuljf(z,a_c[p]);
      } /* a_m */
    } /* p */
  }
#endif
#else
  complex_float x,y,z;
  int p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      /* B <- B*G(a,b) */
      if (n>a_l[p]) {
        F[SIDX(n-1,p)] = csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p]));
      }
    } /* a_m */
  } /* p */
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
      a_a[p] = csubf(cmulf(a_c[p],x), cmulf(a_s[p],y));
      a_b[p] = cnegf(cmulf(a_s[p],z));
    } /* a_m */
  } /* p */
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
      F[SIDX(n,p)] = caddf(cmuljf(x,a_s[p]), cmuljf(y,a_c[p])); 
      D[SIDX(q,p)] = cmuljf(z,a_c[p]);
    } /* a_m */
  } /* p */
#endif
} /* grsvdsf_gks_step1() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 2nd updating step of SVD step iteration: B <- G(a,b)'*B.
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_k[L]       Right index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n,n) Out: B(n,n+1)  
 *   a_b[L]       In: B(n+1,n) Out: B(n,n+2)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_k,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_k
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_step2(
                    complex_float * restrict D,  
                    complex_float * restrict F,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const int16_t       *          a_k,
              const vboolN_4      *          a_m,
              int n, int L )
{
#if 1
        xb_vecN_4xcf32 * restrict Dn_w;
  const xb_vecN_4xcf32 * restrict Dq_r;
        xb_vecN_4xcf32 * restrict Dq_w;
  const xb_vecN_4xcf32 * restrict Fn_r;
        xb_vecN_4xcf32 * restrict Fn_w;
  const xb_vecN_4xcf32 * restrict Fq_r;
        xb_vecN_4xcf32 * restrict Fq_w;
  const long long      * restrict K_rs;
  const xb_vecN_4xcf32 * restrict A_r;
        xb_vecN_4xcf32 * restrict A_w;
  const xb_vecN_4xcf32 * restrict B_r;
        xb_vecN_4xcf32 * restrict B_w;
  const xb_vecN_4xcf32 * restrict C_r;
  const xb_vecN_4xcf32 * restrict S_r;
  const vboolN_4       * restrict M_r;

  int p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  static const int16_t ALIGN(32) isel[BBE_SIMD_WIDTH] = {
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
  };

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;

#if 1
  {
    xb_vecN_4xcf32 a,b,c,d,s;
    vboolN_4 bm;
    Dn_w = (xb_vecN_4xcf32*)&D[SIDX(n,0)];
    A_r = (xb_vecN_4xcf32*)a_a;
    B_r = (xb_vecN_4xcf32*)a_b;
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* D[SIDX(n,p)] = csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p])); */
      BBE_LVN_4XCF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      d = BBE_MULN_4XCF32(c,a); BBE_MULSN_4XCF32(d,s,b);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(d, Dn_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        /* B <- G(a,b).'*B  */
        D[SIDX(n,p)] = csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p]));
      } /* a_m */
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_4xcf32 c,d,f,s,y,z;
    vboolN_4 bm;
    Dq_r = Dq_w = (xb_vecN_4xcf32*)&D[SIDX(q,0)];
    Fn_r = Fn_w = (xb_vecN_4xcf32*)&F[SIDX(n,0)];
    A_w = (xb_vecN_4xcf32*)a_a;
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* y <- F[SIDX(n,p)]; z <- D[SIDX(q,p)]; */
      BBE_LVN_4XCF32_IP(y, Fn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(z, Dq_r, 2*BBE_SIMD_WIDTH);
      /* a_a[p] <- F[SIDX(n,p)] <- csubf(cmulf(a_c[p],y), cmulf(a_s[p],z)); */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      f = BBE_MULN_4XCF32(c,y); BBE_MULSN_4XCF32(f,s,z);
      /* D[SIDX(q,p)] <- caddf(cmuljf(y,a_s[p]), cmuljf(z,a_c[p])); */
      d = BBE_MULJN_4XCF32(y,s); BBE_MULJAN_4XCF32(d,z,c);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(f, Fn_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(f, A_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(d, Dq_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    complex_float y,z;
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
        a_a[p] = F[SIDX(n,p)] = csubf(cmulf(a_c[p],y), cmulf(a_s[p],z));
        D[SIDX(q,p)] = caddf(cmuljf(y,a_s[p]), cmuljf(z,a_c[p]));
      } /* a_m */
    } /* p */
  }
#endif
#if 1
  {
    xb_vecN_4xcf32 c,b,f,s,y,z,w;
    xb_vecN_4x64 sk;
    xb_vecNx16 vn,vk;
    vboolN_4 bm,bn,bmn;
    vselN sel = BBE_MOVVSELNX16(BBE_LVNX16_I((xb_vecNx16*)isel, 0), 0);
    vn = BBE_MOVVA16(n+1);
    Fn_r = (xb_vecN_4xcf32*)&F[SIDX(n,0)];
    Fq_r = Fq_w = (xb_vecN_4xcf32*)&F[SIDX(q,0)];
    Dq_r = (xb_vecN_4xcf32*)&D[SIDX(q,0)];
    K_rs = (long long*)a_k;
    B_w = (xb_vecN_4xcf32*)a_b;
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      BBE_LSN_4X64_IP(sk, K_rs, sizeof(*K_rs));
      vk = BBE_SHFLNX16(BBE_MOVNX16_FROM64(sk), sel);
      bn = BBE_MOVN_4_FROMN(BBE_LTNX16(vn,vk));
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      bmn = BBE_ANDBN_4(bm,bn);
      /* y <- F[SIDX(n,p)]; z <- D[SIDX(q,p)]; */
      BBE_LVN_4XCF32_IP(y, Fn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(z, Dq_r, 2*BBE_SIMD_WIDTH);
      /* if (n<a_k[p]-1) {
       *   w <- F[SIDX(q,p)];
       *   a_b[p] <- cnegf(cmulf(a_s[p],w)); 
       *   F[SIDX(q,p)] <- cmuljf(w,a_c[p]);
       * } */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32T_IP(w, Fq_r, 2*BBE_SIMD_WIDTH, bmn);
      b = BBE_MULMN_4XCF32(s,w,3,8); BBE_MULMASN_4XCF32(b,s,w,2,7);
      f = BBE_MULJN_4XCF32(w,c);
      BBE_SVN_4XCF32T_IP(b, B_w, 2*BBE_SIMD_WIDTH, bmn);
      BBE_SVN_4XCF32T_IP(f, Fq_w, 2*BBE_SIMD_WIDTH, bmn);
    } /* p */
  }
#else
  {
    complex_float y,z,w;
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        if (n<a_k[p]-1) {
          y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
          w = F[SIDX(q,p)]; a_b[p] = cnegf(cmulf(a_s[p],w)); 
          F[SIDX(q,p)] = cmuljf(w,a_c[p]);
        }
      } /* a_m */
    } /* p */
  }
#endif
#else
  complex_float y,z,w;
  int p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      /* B <- G(a,b).'*B  */
      D[SIDX(n,p)] = csubf(cmulf(a_c[p],a_a[p]), cmulf(a_s[p],a_b[p]));
    } /* a_m */
  } /* p */
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
      a_a[p] = F[SIDX(n,p)] = csubf(cmulf(a_c[p],y), cmulf(a_s[p],z));
      D[SIDX(q,p)] = caddf(cmuljf(y,a_s[p]), cmuljf(z,a_c[p]));
    } /* a_m */
  } /* p */
  for ( p=0; p<L; p++ ) {
    if (vbn4Get(a_m+n*LW4, p)) {
      if (n<a_k[p]-1) {
        y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
        w = F[SIDX(q,p)]; a_b[p] = cnegf(cmulf(a_s[p],w)); 
        F[SIDX(q,p)] = cmuljf(w,a_c[p]);
      }
    } /* a_m */
  } /* p */
#endif
} /* grsvdsf_gks_step2() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Accumulate transformations: U <- U*conj(G(a,b)) or V <- V*G(a,b),
 * where U (MxN) and V (NxN) are matrices comprised of left- and right
 * singular vectors. For real data, a single function is used for both
 * U and V updates. For complex data there are two separate functions
 * for matrices U and V.
 * Input:
 *   n            Current position
 *   M,N          Matrix dimensions
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   W[M*N][L]    Matrix of left-singular (U) or right-singular (V) 
 *                orthonormal vectors
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   W,a_c,a_s,a_m
 *      Must not overlap
 *   W,a_c,a_s
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_accumLeft(
                    complex_float * restrict W,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const vboolN_4      *          a_m,
              int n, int M, int N, int L )
{
#if 1
  const xb_vecN_4xcf32 * restrict Wn_r;
        xb_vecN_4xcf32 * restrict Wn_w;
  const xb_vecN_4xcf32 * restrict Wp_r;
        xb_vecN_4xcf32 * restrict Wp_w;
  const xb_vecN_4xcf32 * restrict C_r;
  const xb_vecN_4xcf32 * restrict S_r;
  const vboolN_4       * restrict M_r;

  xb_vecN_4xcf32 c,s,x,y,w,z;
  vboolN_4 bm;
  int i,p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( i=0; i<M; i++ ) {
    Wn_r = Wn_w = (xb_vecN_4xcf32*)&W[SIDX(i*N+n,0)];
    Wp_r = Wp_w = (xb_vecN_4xcf32*)&W[SIDX(i*N+q,0)];
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* x <- W[SIDX(i*N+n,p)]; y <- W[SIDX(i*N+q,p)]; */
      BBE_LVN_4XCF32_IP(x, Wn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(y, Wp_r, 2*BBE_SIMD_WIDTH);
      /* W[SIDX(i*N+n,p)] <- csubf(cmuljf(x,a_c[p]), cmuljf(y,a_s[p])); */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      z = BBE_MULJN_4XCF32(x,c); BBE_MULJSN_4XCF32(z,y,s);
      /* W[SIDX(i*N+q,p)] <- caddf(cmulf(x,a_s[p]), cmulf(y,a_c[p])); */
      w = BBE_MULN_4XCF32(x,s); BBE_MULAN_4XCF32(w,y,c);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(z, Wn_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(w, Wp_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
  } /* i */
#else
  complex_float x,y;
  int i,p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( i=0; i<M; i++ ) {
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        x = W[SIDX(i*N+n,p)]; y = W[SIDX(i*N+q,p)];
        /* W <- W*conj(G(a,b)) */
        W[SIDX(i*N+n,p)] = csubf(cmuljf(x,a_c[p]), cmuljf(y,a_s[p]));
        W[SIDX(i*N+q,p)] = caddf(cmulf(x,a_s[p]), cmulf(y,a_c[p]));
      } /* a_m */
    } /* p */
  } /* i */
#endif
} /* grsvdsf_gks_accumLeft() */

void grsvdsf_gks_accumRight(
                    complex_float * restrict W,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const vboolN_4      *          a_m,
              int n, int N, int L )
{
#if 1
  const xb_vecN_4xcf32 * restrict Wn_r;
        xb_vecN_4xcf32 * restrict Wn_w;
  const xb_vecN_4xcf32 * restrict Wp_r;
        xb_vecN_4xcf32 * restrict Wp_w;
  const xb_vecN_4xcf32 * restrict C_r;
  const xb_vecN_4xcf32 * restrict S_r;
  const vboolN_4       * restrict M_r;

  xb_vecN_4xcf32 c,s,x,y,w,z;
  vboolN_4 bm;
  int i,p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( i=0; i<N; i++ ) {
    Wn_r = Wn_w = (xb_vecN_4xcf32*)&W[SIDX(i*N+n,0)];
    Wp_r = Wp_w = (xb_vecN_4xcf32*)&W[SIDX(i*N+q,0)];
    C_r = (xb_vecN_4xcf32*)a_c;
    S_r = (xb_vecN_4xcf32*)a_s;
    M_r = a_m + n*LW4;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/4); p++ ) {
      /* x <- W[SIDX(i*N+n,p)]; y <- W[SIDX(i*N+q,p)]; */
      BBE_LVN_4XCF32_IP(x, Wn_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(y, Wp_r, 2*BBE_SIMD_WIDTH);
      /* W[SIDX(i*N+n,p)] <- csubf(cmulf(x,a_c[p]), cmulf(y,a_s[p])); */
      BBE_LVN_4XCF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_4XCF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      z = BBE_MULN_4XCF32(x,c); BBE_MULSN_4XCF32(z,y,s);
      /* W[SIDX(i*N+q,p)] <- caddf(cmuljf(x,a_s[p]), cmuljf(y,a_c[p])); */
      w = BBE_MULJN_4XCF32(x,s); BBE_MULJAN_4XCF32(w,y,c);
      BBE_LBN_4_IP(bm, M_r, sz_vbn4);
      BBE_SVN_4XCF32T_IP(z, Wn_w, 2*BBE_SIMD_WIDTH, bm);
      BBE_SVN_4XCF32T_IP(w, Wp_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
  } /* i */
#else
  complex_float x,y;
  int i,p,q;
  int LW4 = L/(BBE_SIMD_WIDTH/4);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( i=0; i<N; i++ ) {
    for ( p=0; p<L; p++ ) {
      if (vbn4Get(a_m+n*LW4, p)) {
        x = W[SIDX(i*N+n,p)]; y = W[SIDX(i*N+q,p)];
        /* W <- W*G(a,b) */
        W[SIDX(i*N+n,p)] = csubf(cmulf(x,a_c[p]), cmulf(y,a_s[p]));
        W[SIDX(i*N+q,p)] = caddf(cmuljf(x,a_s[p]), cmuljf(y,a_c[p]));
      } /* a_m */
    } /* p */
  } /* i */
#endif
} /* grsvdsf_gks_accumRight() */

#endif /* HAVE_VFPU */

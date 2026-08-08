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
    NatureDSP_Baseband library. Fitting and Interpolation Routines
    Complex Polynomial Roots Finding
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <string.h>
#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Real and complex arithmetic primitives optimized for BBEN VFPU */
#include "vfpu_math.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_fit.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#define BBEN_2   (BBE_SIMD_WIDTH/2)
#define BBEN_4   (BBE_SIMD_WIDTH/4)

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )

#if HAVE_VFPU

#if 0

#include <math.h>

static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
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

/* Cheap estimation of complex number's magnitude, single precision
 * floating-point. */
static float32_t cmagf( complex_float x )
{
  return ( fabsf(crealf(x)) + fabsf(cimagf(x)) );
}

/* Negation of a complex number, single precision floating-point. */
static complex_float cnegf( complex_float x )
{
  return ( _makecomplexf(-crealf(x),-cimagf(x)) );

} /* cnegf() */

#endif

/* Polynomial coefficients normalization: f[0..N-1] = c[1..N]/c[0]
 * Input:
 *   N       Polynomial degree
 *   c[N+1]  Polynomial coefficients
 * Output:
 *   f[N]    Normalized polynomial coeffcicients
 * Restrictions:
 *   c[0]    Either real or imaginary component of the leading coefficient
 *           must be non-zero. */
static void normCoef( complex_float * f, const complex_float * c, int N );

/* A variant of balancing algorithm tailored for the companion matrix
 * Input:
 *   N          Polynomial degree
 *   iterNum    Max number of balancing iterations
 * Input/Output:
 *   f[N]       Polynomial coeffs ordered in descending poweres. The leading 
 *              coeff (at power N) equals 1 and is not stored
 * Output:
 *   b[N]       Optional array of scale factors stored as integral powers of 2
 *   d[N-1]     1st subdiagonal of balamced companion matrix
 * Restrictions:
 *   f[N-1]!=0  The free term must be non-zero
 *   N>1        Polynomial degree must be at least 2 */
static void balCompan( int * b, float32_t * d, complex_float * f, int N, int iterNum );

/*---------------------------------------------------------------------------
Find Roots of a Polynomial with Complex/Real Coefficients

Given input array c[N+1] containing N+1 complex/real coefficients of a 
degree N polynomial, find its N roots (possibly with repetitions) and
store results to the output array r[N].

Representation: IEEE-754 Std single precision floating-point data

Parameters:
Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by functions 
            [r]prootf_getScratchSize(N)
Input:
  N         Polynomial degree
  c[N+1]    Polynomial coefficients in descending powers order, i.e. 
            p(x) = c[0]*x^N + ... + c[N-1]*x + c[N]
Output:
  r[N]      Roots of the polynomial

Restrictions:
  pScr,r,c  Must not overlap and must be aligned on 32-byte boundary
  N>0       Polynomial degree must be positive
  c[0]!=0   The leading coefficient must be non-zero
---------------------------------------------------------------------------*/
/* Complex polynomial coefficients */
void prootf (
          void * pScr, 
          complex_float * restrict r,
    const complex_float * restrict c,
    int N )
{
#if 1
  float32_t *D;
  complex_float *F,*C;
  xtcomplexfloat * restrict Cw;
  const xtcomplexfloat * restrict _c = (xtcomplexfloat*)c;
  xtcomplexfloat * restrict _r = (xtcomplexfloat*)r;
  const xtcomplexfloat * restrict Fr;

  const xtcomplexfloat c0f = BBE_CONSTCF32(0);

  int SD,SF,SC;
  int i,stp;

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(r   , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(c   , 2*BBE_SIMD_WIDTH);
  NASSERT(N>0);

  NASSERT(!((float32_t)BBE_CREALCF32(_c[0])==0.f) ||
          !((float32_t)BBE_CIMAGCF32(_c[0])==0.f));

  SD = (N-1+BBEN_2-1)/BBEN_2*BBEN_2;
  SF = (N+BBEN_4-1)/BBEN_4*BBEN_4;
  SC = (N*(N+3)/2-1+BBEN_4-1)/BBEN_4*BBEN_4;

  {
    void * p = pScr;
    /* Partition the scratch area. */
    D = (float32_t    *)p; p = D+SD;
    F = (complex_float*)p; p = F+SF;
    C = (complex_float*)p; p = C+SC;
    /* Make sure that scratch arrays fit into the requested space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)prootf_getScratchSize(N) );
  }
  /* Normalize trailing coeffs by the leading coeff. */
  normCoef(F,c,N);
  /* Truncate trailing zeros and deflate the problem. */
  Fr = (xtcomplexfloat*)&F[0];
  while (N>0 && BBE_MOVAB1(XT_OEQ_S(IT_CMAGCF32(Fr[N-1]), XT_CONST_S(0)))) {
    _r[N-1] = c0f; N--;
  }
  if (N>1) {
    /* Apply a variant of balancing algorithm tailored for the companion matrix. */
    balCompan(0,D,F,N,5);
    /* Construct the companion matrix. */
    memset(C, 0, SC*sz_f32c);
    Fr = (xtcomplexfloat*)&F[0];
    Cw = (xtcomplexfloat*)&C[HIDX(N-1,N-1)];
    stp = 2;
    for ( i=0; i<N; i++ ) {
      *Cw = BBE_NEGCF32(*Fr++);
      Cw -= stp++;
    }
    Cw = (xtcomplexfloat*)&C[HIDX(1,0)];
    stp = N;
    for ( i=0; i<N-1; i++ ) {
      *Cw = BBE_CMPLXF32(XT_CONST_S(0), D[i]);
      Cw += stp--;
    }
    /* Compute eigenvalues of the companion matrix by QR algorithm. */
    eigen_hqr_f(r,C,0,0,N-1,N);
  } else {
    _r[0] = BBE_NEGCF32(*(xtcomplexfloat*)&F[0]);
  }
#else
  float32_t *D;
  complex_float *F,*C;
  int SD,SF,SC;
  int i;
  const complex_float c0f = _makecomplexf(0.f,0.f);

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(r   , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(c   , 2*BBE_SIMD_WIDTH);
  NASSERT(N>0);
  NASSERT(!(crealf(c[0])==0.f && cimagf(c[0])==0.f));

  SD = (N-1+BBEN_2-1)/BBEN_2*BBEN_2;
  SF = (N+BBEN_4-1)/BBEN_4*BBEN_4;
  SC = (N*(N+3)/2-1+BBEN_4-1)/BBEN_4*BBEN_4;

  {
    void * p = pScr;
    /* Partition the scratch area. */
    D = (float32_t    *)p; p = D+SD;
    F = (complex_float*)p; p = F+SF;
    C = (complex_float*)p; p = C+SC;
    /* Make sure that scratch arrays fit into the requested space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)prootf_getScratchSize(N) );
  }
  /* Normalize trailing coeffs by the leading coeff. */
  normCoef(F,c,N);
  /* Truncate trailing zeros and deflate the problem. */
  while (N>0 && cmagf(F[N-1])==0.f) {
    r[N-1] = c0f; N--;
  }
  if (N>1) {
    /* Apply a variant of balancing algorithm tailored for the companion matrix. */
    balCompan(0,D,F,N,5);
    /* Construct the companion matrix. */
    memset(C, 0, SC*sz_f32c);
    for ( i=0; i<N-1; i++ ) { 
      C[HIDX(i+1,i)] = _makecomplexf(D[i],0.f);
      C[HIDX(i,N-1)] = cnegf(F[N-(i+1)]);
    }
    C[HIDX(N-1,N-1)] = cnegf(F[0]);
    /* Compute eigenvalues of the companion matrix by QR algorithm. */
    eigen_hqr_f(r,C,0,0,N-1,N);
  } else {
    r[0] = cnegf(F[0]);
  }
#endif
} /* prootf() */

/* Return the scratch area size, in bytes. */
size_t prootf_getScratchSize ( int N ) 
{ 
  NASSERT(N>0);
  return ( (N-1        +BBEN_2-1)/BBEN_2*BBEN_2*sz_f32  +  /* D: 1st subdiagonal of the companion matrix  */
           (N          +BBEN_4-1)/BBEN_4*BBEN_4*sz_f32c +  /* F: balanced polynomial coeffs               */
           (N*(N+3)/2-1+BBEN_4-1)/BBEN_4*BBEN_4*sz_f32c ); /* C: companion matrix, elements below the 1st *
                                                            *    subdiagonal are not stored               */
}

/* Polynomial coefficients normalization: f[0..N-1] = c[1..N]/c[0]
 * Input:
 *   N       Polynomial degree
 *   c[N+1]  Polynomial coefficients
 * Output:
 *   f[N]    Normalized polynomial coeffcicients
 * Restrictions:
 *   c[0]    Either real or imaginary component of the leading coefficient
 *           must be non-zero. */
void normCoef( complex_float * f, const complex_float * c, int N )
{
#if 1
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

  xb_vecN_2xf32 * restrict fw;
  const xb_vecN_2xf32 * restrict cr;
  valign fva, cva;
  xtfloat _g,_h,yr,yi;
  xb_vecN_2xf32 g,h,s,x,z;
  int n, nb;
  yr = BBE_CREALCF32(*(xtcomplexfloat*)&c[0]);
  yi = BBE_CIMAGCF32(*(xtcomplexfloat*)&c[0]);
  NASSERT(!((float32_t)yr==0.f) || !((float32_t)yi==0.f));
  if (BBE_MOVAB1(XT_OLT_S(yi,yr))) {
    _h = XT_DIV_S(yi,yr);
    _g = yr; XT_MADD_S(_g,_h,yi);
    g = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_g),0);
    h = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_h),0);
    fw = (xb_vecN_2xf32*)&f[0];
    cr = (xb_vecN_2xf32*)&c[1];
    fva = BBE_ZALIGN();
    cva = BBE_LAN_2XF32_PP(cr);
    for ( n=0; n<N/(BBE_SIMD_WIDTH/4); n++ ) {
      BBE_LAN_2XF32_IP(x, cva, cr); 
      s = x; BBE_MULMASN_2XF32(s,h,x,2,1);
      z = BBE_DIVN_2XF32(s,g);
      BBE_SAN_2XF32_IP(z, fva, fw);
    }
    if (N&(BBE_SIMD_WIDTH/4-1)) {
      __Pragma("frequency_hint frequent");
      nb = (N&(BBE_SIMD_WIDTH/4-1))*sz_f32c;
      BBE_LAVN_2XF32_XP(x, cva, cr, nb);
      s = x; BBE_MULMASN_2XF32(s,h,x,2,1);
      z = BBE_DIVN_2XF32(s,g);
      BBE_SAVN_2XF32_XP(z, fva, fw, nb);
    }
    BBE_SAVN_2XF32POS_FP(fva, fw);
  } else {
    _h = XT_DIV_S(yr,yi);
    _g = yi; XT_MADD_S(_g,_h,yr);
    g = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_g),0);
    h = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_h),0);
    h = BBE_SELN_2XF32I(BBE_CONSTN_2XF32(1), h, BBE_SELI_INTERLEAVE_2_EVEN);
    fw = (xb_vecN_2xf32*)&f[0];
    cr = (xb_vecN_2xf32*)&c[1];
    fva = BBE_ZALIGN();
    cva = BBE_LAN_2XF32_PP(cr);
    for ( n=0; n<N/(BBE_SIMD_WIDTH/4); n++ ) {
      BBE_LAN_2XF32_IP(x, cva, cr); 
      s = BBE_MULMN_2XF32(h,x,2,11);
      BBE_MULMASN_2XF32(s,h,x,0,4);
      z = BBE_DIVN_2XF32(s,g);
      BBE_SAN_2XF32_IP(z, fva, fw);
    }
    if (N&(BBE_SIMD_WIDTH/4-1)) {
      __Pragma("frequency_hint frequent");
      nb = (N&(BBE_SIMD_WIDTH/4-1))*sz_f32c;
      BBE_LAVN_2XF32_XP(x, cva, cr, nb);
      s = BBE_MULMN_2XF32(h,x,2,11);
      BBE_MULMASN_2XF32(s,h,x,0,4);
      z = BBE_DIVN_2XF32(s,g);
      BBE_SAVN_2XF32_XP(z, fva, fw, nb);
    }
    BBE_SAVN_2XF32POS_FP(fva, fw);
  }
#elif 0
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

  xb_vecN_2xf32 * restrict fw;
  const xb_vecN_2xf32 * restrict cr;
  valign fva, cva;
  xtfloat yr,yi;
  xtfloat p,q,_g,_h,_r;
  xb_vecN_2xf32 g,h,r,s,x,z;
  int n, nb;
  yr = BBE_CREALCF32(*(xtcomplexfloat*)&c[0]);
  yi = BBE_CIMAGCF32(*(xtcomplexfloat*)&c[0]);
  NASSERT(!((float32_t)XT_ABS_S(yr)<FLT_MIN) || !((float32_t)XT_ABS_S(yi)<FLT_MIN));
  if (BBE_MOVAB1(XT_OLT_S(yi,yr))) {
    p = XT_RECIP0_S(yr);
    q = XT_CONST_S(1); XT_MSUB_S(q,p,yr); XT_MADD_S(p,q,p);
    _h = XT_MUL_S(p,yi); q = yi; XT_MSUB_S(q,_h,yr); XT_MADD_S(_h,q,p);
    _g = yr; XT_MADD_S(_g,_h,yi);
    _r = XT_RECIP0_S(_g);
    q = XT_CONST_S(1); XT_MSUB_S(q,_r,_g); XT_MADD_S(_r,q,_r);
    g = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_g),0);
    h = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_h),0);
    r = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_r),0);
    fw = (xb_vecN_2xf32*)&f[0];
    cr = (xb_vecN_2xf32*)&c[1];
    fva = BBE_ZALIGN();
    cva = BBE_LAN_2XF32_PP(cr);
    __Pragma("no_unroll");
    for ( n=0; n<N/(BBE_SIMD_WIDTH/4); n++ ) {
      BBE_LAN_2XF32_IP(x, cva, cr); 
      s = x; BBE_MULMASN_2XF32(s,h,x,2,1);
      z = BBE_MULN_2XF32(s,r);
      BBE_MULSN_2XF32(s,z,g);
      BBE_MULAN_2XF32(z,s,r);
      BBE_SAN_2XF32_IP(z, fva, fw);
    }
    if (N&(BBE_SIMD_WIDTH/4-1)) {
      __Pragma("frequency_hint frequent");
      nb = (N&(BBE_SIMD_WIDTH/4-1))*sz_f32c;
      BBE_LAVN_2XF32_XP(x, cva, cr, nb);
      s = x; BBE_MULMASN_2XF32(s,h,x,2,1);
      z = BBE_MULN_2XF32(s,r);
      BBE_MULSN_2XF32(s,z,g);
      BBE_MULAN_2XF32(z,s,r);
      BBE_SAVN_2XF32_XP(z, fva, fw, nb);
    }
    BBE_SAVN_2XF32POS_FP(fva, fw);
  } else {
    p = XT_RECIP0_S(yi);
    q = XT_CONST_S(1); XT_MSUB_S(q,p,yi); XT_MADD_S(p,q,p);
    _h = XT_MUL_S(p,yr); q = yr; XT_MSUB_S(q,_h,yi); XT_MADD_S(_h,q,p);
    _g = yi; XT_MADD_S(_g,_h,yr);
    _r = XT_RECIP0_S(_g);
    q = XT_CONST_S(1); XT_MSUB_S(q,_r,_g); XT_MADD_S(_r,q,_r);
    g = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_g),0);
    h = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_h),0);
    r = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(_r),0);
    fw = (xb_vecN_2xf32*)&f[0];
    cr = (xb_vecN_2xf32*)&c[1];
    fva = BBE_ZALIGN();
    cva = BBE_LAN_2XF32_PP(cr);
    __Pragma("no_unroll");
    for ( n=0; n<N/(BBE_SIMD_WIDTH/4); n++ ) {
      BBE_LAN_2XF32_IP(x, cva, cr); 
      s = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x, BBE_SHFLI_SWAP_2));
      BBE_MULMASN_2XF32(s,h,x,0,4);
      z = BBE_MULN_2XF32(s,r);
      BBE_MULSN_2XF32(s,z,g);
      BBE_MULAN_2XF32(z,s,r);
      BBE_SAN_2XF32_IP(z, fva, fw);
    }
    if (N&(BBE_SIMD_WIDTH/4-1)) {
      __Pragma("frequency_hint frequent");
      nb = (N&(BBE_SIMD_WIDTH/4-1))*sz_f32c;
      BBE_LAVN_2XF32_XP(x, cva, cr, nb);
      s = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x, BBE_SHFLI_SWAP_2));
      BBE_MULMASN_2XF32(s,h,x,0,4);
      z = BBE_MULN_2XF32(s,r);
      BBE_MULSN_2XF32(s,z,g);
      BBE_MULAN_2XF32(z,s,r);
      BBE_SAVN_2XF32_XP(z, fva, fw, nb);
    }
    BBE_SAVN_2XF32POS_FP(fva, fw);
  }
#else
  int i;
  NASSERT(!(fabsf(crealf(c[0]))<FLT_MIN) || !(fabsf(cimagf(c[0]))<FLT_MIN));
  for ( i=0; i<N; i++ ) {
    f[i] = _cdivf(c[i+1],c[0]);
  }
#endif
} /* normCoef() */

/* A variant of balancing algorithm tailored for the companion matrix
 * Input:
 *   N          Polynomial degree
 *   iterNum    Max number of balancing iterations
 * Input/Output:
 *   f[N]       Polynomial coeffs ordered in descending poweres. The leading 
 *              coeff (at power N) equals 1 and is not stored
 * Output:
 *   b[N]       Optional array of scale factors stored as integral powers of 2
 *   d[N-1]     1st subdiagonal of balamced companion matrix
 * Restrictions:
 *   f[N-1]!=0  The free term must be non-zero
 *   N>1        Polynomial degree must be at least 2 */
void balCompan( int * b, float32_t * d, complex_float * f, int N, int iterNum )
{
  /*
   * See "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors"
   * by B.N. Partlett and C. Reinsch, Handbook for Automatie Computation,
   * Vol.II Linear Algebra, Contribution II/11.
   * MATLAB reference code:
   *
   *   function [b,d,f] = balCompan(f,iterNum)
   *   rmin = realmin(class(f));
   *   emin = log2(rmin);
   *   N = length(f);
   *   d = cast(ones(1,N-1),'like',f);
   *   b = zeros(1,N);
   *   iterCnt = 0;
   *   while true
   *     conv = true;
   *     for i=1:N-1
   *       r = cmag(f(N+1-i)); if i>1, r = r+d(i-1); end;
   *       c = d(i);
   *       assert(c>0);
   *       if r>0 % Condition is false iff f[N]==0 (the free term is zero)
   *         if r>=rmin, [u,p] = log2(r); else u = 0.5; p = emin; end;
   *         if c>=rmin, [v,q] = log2(c); else v = 0.5; q = emin; end;
   *         e = floor((p-q+(u>v))/2);
   *         assert((p==emin)||(q==emin)||((2^(2*e-1)<r/c)&&(r/c<=2^(2*e+1))));
   *         if pow2(r,-e)+pow2(c,e)<0.95*(r+c)
   *           conv = false;
   *           b(i) = b(i)+e;
   *           if (i>1), d(i-1) = d(i-1)*2^-e; end;
   *           f(N+1-i) = f(N+1-i)*2^-e;
   *           d(i) = d(i)*2^e;
   *         end
   *       end
   *     end
   *     r = d(N-1);
   *     c = 0; for j=2:N, c = c+cmag(f(j)); end;
   *     assert((r>0)&&(c>0));
   *     if r>=rmin, [u,p] = log2(r); else u = 0.5; p = emin; end;
   *     if c>=rmin, [v,q] = log2(c); else v = 0.4; q = emin; end;
   *     e = floor((p-q+(u>v))/2);
   *     assert((p==emin)||(q==emin)||((2^(2*e-1)<r/c)&&(r/c<=2^(2*e+1))));
   *     if pow2(r,-e)+pow2(c,e)<0.95*(r+c)
   *       conv = false;
   *       b(N) = b(N)+e;
   *       d(N-1) = d(N-1)*2^-e;
   *       f(2:N) = f(2:N)*2^e;
   *     end
   *     iterCnt = iterCnt + 1;
   *     if conv || iterCnt>=iterNum, break; end;
   *   end
   */
#if 1
  const xtfloat        * restrict dr;
        xtfloat        * restrict dw;
  const xtcomplexfloat * restrict fr = (xtcomplexfloat*)f;
        xtcomplexfloat * restrict fw;
  const xtfloat        * restrict _fr;
        xtfloat        * restrict _fw;
  xtcomplexfloat f0;
  float32_t d0,d1;
  float32_t c,r,g,h,s,t;
  uint32_t ru,cu;
  uint32_t p,q,u,v;
  int32_t e;
  vbool1 blt, bconv;
  int iterCnt=0;
  int i,j;

  NASSERT(N>1);
  NASSERT(!((float32_t)IT_CMAGCF32(fr[N-1])==0.f));

  for ( i=0; i<N-1; i++ ) d[i] = XT_CONST_S(1);

  do {
    bconv = BBE_MOVBA1(1);
    dr = dw = (xtfloat*)&d[0];
    fr = fw = (xtcomplexfloat*)&f[N-1];
    /*--------------------------------------------------------*
     * Process rows/cols 0                                    */
    xtfloat_loadip(d0, dr, sz_f32);
    xtcomplexfloat_loadip(f0, fr, -(int)sz_f32c);
    r = IT_CMAGCF32(f0);
    ru = XT_RFR(r); cu = XT_RFR(d0);
    p = ru>>23; u = ru & ((1U<<23)-1); XT_MOVEQZ(u,0,p);
    q = cu>>23; v = cu & ((1U<<23)-1); XT_MOVEQZ(v,0,q);
    e = (int32_t)(p-q+(u>v))>>1;
    /* When balancing a companion matrix, peak exponents are never 
     * encountered. Thus there is no need to take care of subnormal
     * factor of 2^-127. */
    g = XT_WFR((uint32_t)(e+127)<<23);
    h = XT_WFR((uint32_t)(127-e)<<23);
    s = XT_MUL_S(r,h); XT_MADD_S(s,d0,g);
    t = XT_ADD_S(r,d0); t = XT_MUL_S(0.95f,t);
    blt = XT_OLT_S(s,t);
    IT_RCMULCF32T(f0,h,f0,blt);
    BBE_MULF32T(d0,d0,g,blt);
    xtcomplexfloat_storeip(f0, fw, -(int)sz_f32c);
    /*--------------------------------------------------------*
     * Process rows/cols 1..N-2                               */
    for ( i=1; i<N-1; i++ ) {
      xtcomplexfloat_loadip(f0, fr, -(int)sz_f32c);
      xtfloat_loadip(d1, dr, sz_f32);
      r = XT_ADD_S(IT_CMAGCF32(f0), d0);
      ru = XT_RFR(r); cu = XT_RFR(d1);
      p = ru>>23; u = ru & ((1U<<23)-1); XT_MOVEQZ(u,0,p);
      q = cu>>23; v = cu & ((1U<<23)-1); XT_MOVEQZ(v,0,q);
      e = (int32_t)(p-q+(u>v))>>1;
      g = XT_WFR((uint32_t)(e+127)<<23);
      h = XT_WFR((uint32_t)(127-e)<<23);
      s = XT_MUL_S(r,h); XT_MADD_S(s,d1,g);
      t = XT_ADD_S(r,d1); t = XT_MUL_S(0.95f,t);
      blt = XT_OLT_S(s,t);
      BBE_MULF32T(d0,d0,h,blt);
      IT_RCMULCF32T(f0,h,f0,blt);
      BBE_MULF32T(d1,d1,g,blt);
      xtfloat_storeip(d0, dw, sz_f32);
      xtcomplexfloat_storeip(f0, fw, -(int)sz_f32c);
      bconv = BBE_OPERATOR_ANDB1(bconv,BBE_NOTB1(blt));
      d0 = d1;
    } /* i */
    __Pragma("no_reorder");
    /*--------------------------------------------------------*
     * Process the last row/col                               */
    _fr = (xtfloat*)fr+2;
    __Pragma("loop_count min=3");
    for ( c=0.f, j=0; j<2*(N-1); j++ ) { t = *_fr++; c += (float32_t)XT_ABS_S(t); };
    ru = XT_RFR(d0); cu = XT_RFR(c);
    p = ru>>23; u = ru & ((1U<<23)-1); XT_MOVEQZ(u,0,p);
    q = cu>>23; v = cu & ((1U<<23)-1); XT_MOVEQZ(v,0,q);
    e = (int32_t)(p-q+(u>v))>>1;
    g = XT_WFR((uint32_t)(e+127)<<23);
    h = XT_WFR((uint32_t)(127-e)<<23);
    s = XT_MUL_S(d0,h); XT_MADD_S(s,c,g);
    t = XT_ADD_S(d0,c); t = XT_MUL_S(0.95f,t);
    blt = XT_OLT_S(s,t);
    BBE_MULF32T(d0,d0,h,blt);
    xtfloat_storeip(d0, dw, sz_f32);
    if (BBE_MOVAB1(blt)) {
      _fr = _fw = (xtfloat*)_fr-1;
      for ( c=0.f, j=0; j<2*(N-1); j++ ) *_fw-- = (float32_t)*_fr--*g;
    }
    bconv = BBE_OPERATOR_ANDB1(bconv,BBE_NOTB1(blt));
    __Pragma("no_reorder");
  } while (++iterCnt<iterNum && !BBE_MOVAB1(bconv));
#else
  union ufloat32uint32 r,c,g,h;
  uint32_t u,v;
  int e,p,q;
  int conv,iterCnt=0;
  int i,j;

  NASSERT(N>0);
  NASSERT(!(cmagf(f[N-1])==0.f));

  for ( i=0; i<N-1; i++ ) d[i] = 1.f;
  if (b) for ( i=0; i<N; i++ ) b[i] = 0;

  do {
    conv = 1;
    /*--------------------------------------------------------*
     * Process rows/cols 0..N-2                               */
    for ( i=0; i<N-1; i++ ) {
      r.f = cmagf(f[N-1-i]); if (i>0) r.f += d[i-1];
      c.f = d[i];
      NASSERT(!(r.f==0.f) && !(c.f==0.f));
      p = (int)(r.u>>23); u = ( p>0 ? (r.u&((1U<<23)-1)) : 0 );
      q = (int)(c.u>>23); v = ( q>0 ? (c.u&((1U<<23)-1)) : 0 );
      e = (p-q+(u>v))>>1;
      /* When balancing a companion matrix, peak exponents are never 
       * encountered. Thus there is no need to take care of subnormal
       * factor of 2^-127. */
      NASSERT((e>-127 && e<127) || !isfinite(r.f) || !isfinite(c.f));
      g.u = (uint32_t)(127+e)<<23;
      h.u = (uint32_t)(127-e)<<23;
      if (r.f*h.f+c.f*g.f < 0.95f*(r.f+c.f)) {
        conv = 0;
        if (b) b[i] += e;
        if (i>0) d[i-1] *= h.f;
        f[N-1-i] = rcmulf(h.f, f[N-1-i]);
        d[i] *= g.f;
      }
    } /* i */
    /*--------------------------------------------------------*
     * Process the last row/col                               */
    r.f = d[N-2];
    for ( c.f=0.f, j=1; j<N; j++ ) c.f += cmagf(f[j]);
    NASSERT(!(r.f==0.f) && !(c.f==0.f));
    p = (int)(r.u>>23); u = ( p>0 ? (r.u&((1U<<23)-1)) : 0 );
    q = (int)(c.u>>23); v = ( q>0 ? (c.u&((1U<<23)-1)) : 0 );
    e = (p-q+(u>v))>>1;
    NASSERT((e>-127 && e<127) || !isfinite(r.f) || !isfinite(c.f));
    g.u = (uint32_t)(127+e)<<23;
    h.u = (uint32_t)(127-e)<<23;
    if (r.f*h.f+c.f*g.f < 0.95f*(r.f+c.f)) {
      conv = 0;
      if (b) b[N-1] += e;
      d[N-2] *= h.f;
      for ( j=1; j<N; j++ ) f[j] = rcmulf(g.f, f[j]);
    }
  } while (++iterCnt<iterNum && !conv);
#endif
} /* balCompan() */

#else /* HAVE_VFPU */

DISCARD_FUN( void, prootf, ( void * pScr, 
                             complex_float * restrict r,
                       const complex_float * restrict c,
                       int N ) )

size_t prootf_getScratchSize ( int N ) { return (0); }

#endif /* HAVE_VFPU */

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
    Real Polynomial Roots Finding
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_fit.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#define sz_f32   sizeof(float32_t)

#define BBEN_2   (BBE_SIMD_WIDTH/2)

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )

#if HAVE_VFPU

#include <math.h>

#if 0

static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
}

#endif

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
static void balCompan( int * b, float32_t * d, float32_t * f, int N, int iterNum );

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
/* Real polynomial coefficients */
void rprootf ( 
          void * pScr, 
          complex_float * restrict r,
    const float32_t     * restrict c,
    int N )
{
#if 1
  float32_t *D,*F,*C;
  float32_t * restrict Cw;
  xtcomplexfloat * restrict _r = (xtcomplexfloat*)r;

  const xtcomplexfloat c0f = BBE_CONSTCF32(0);

  int SD,SF,SC;
  int i,stp;

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(r   , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(c   , 2*BBE_SIMD_WIDTH);
  NASSERT(N>0);
  NASSERT(!(c[0]==0.f));

  SD = (N-1+BBEN_2-1)/BBEN_2*BBEN_2;
  SF = (N+BBEN_2-1)/BBEN_2*BBEN_2;
  SC = (N*(N+3)/2-1+BBEN_2-1)/BBEN_2*BBEN_2;

  {
    void * p = pScr;
    /* Partition the scratch area. */
    D = (float32_t*)p; p = D+SD;
    F = (float32_t*)p; p = F+SF;
    C = (float32_t*)p; p = C+SC;
    /* Make sure that scratch arrays fit into the requested space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)rprootf_getScratchSize(N) );
  }

  /* Normalize trailing coeffs by the leading coeff. */
  for ( i=0; i<N; i++ ) {
    F[i] = c[i+1]/c[0];
  }
  /* Truncate trailing zeros and deflate the problem. */
  while (N>0 && 0.f==F[N-1]) {
    _r[N-1] = c0f; N--;
  }
  if (N>1) {
    /* Apply a variant of balancing algorithm tailored for the companion matrix. */
    balCompan(0,D,F,N,5);
    /* Construct the companion matrix in upper-Hessenberg form. */
    memset(C, 0, SC*sz_f32);
    Cw = &C[HIDX(N-1,N-1)];
    stp = 2;
    for ( i=0; i<N; i++ ) {
      *Cw = -F[i]; Cw -= stp++;
    }
    Cw = &C[HIDX(1,0)];
    stp = N;
    for ( i=0; i<N-1; i++ ) {
      *Cw = D[i]; Cw += stp--;
    }
    /* Compute eigenvalues of the companion matrix by QR algorithm. */
    reigen_hqr_f(r,C,0,0,N-1,N);
  } else {
    _r[0] = BBE_CMPLXF32(XT_CONST_S(0),-F[0]);
  }
#else
  float32_t *D,*F,*C;
  int SD,SF,SC;
  int i;
  const complex_float c0f = _makecomplexf(0.f,0.f);

  NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(r   , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(c   , 2*BBE_SIMD_WIDTH);
  NASSERT(N>0);
  NASSERT(!(c[0]==0.f));

  SD = (N-1+BBEN_2-1)/BBEN_2*BBEN_2;
  SF = (N+BBEN_2-1)/BBEN_2*BBEN_2;
  SC = (N*(N+3)/2-1+BBEN_2-1)/BBEN_2*BBEN_2;

  {
    void * p = pScr;
    /* Partition the scratch area. */
    D = (float32_t*)p; p = D+SD;
    F = (float32_t*)p; p = F+SF;
    C = (float32_t*)p; p = C+SC;
    /* Make sure that scratch arrays fit into the requested space. */
    NASSERT( (uint8_t*)p - (uint8_t*)pScr <= (int)rprootf_getScratchSize(N) );
  }

  /* Normalize trailing coeffs by the leading coeff. */
  for ( i=0; i<N; i++ ) {
    F[i] = c[i+1]/c[0];
  }
  /* Truncate trailing zeros and deflate the problem. */
  while (N>0 && 0.f==F[N-1]) {
    r[N-1] = c0f; N--;
  }
  if (N>1) {
    /* Apply a variant of balancing algorithm tailored for the companion matrix. */
    balCompan(0,D,F,N,5);
    /* Construct the companion matrix. */
    memset(C, 0, SC*sz_f32);
    for ( i=0; i<N-1; i++ ) { 
      C[HIDX(i+1,i)] = D[i];
      C[HIDX(i,N-1)] = -F[N-(i+1)];
    }
    C[HIDX(N-1,N-1)] = -F[0];
    /* Compute eigenvalues of the companion matrix by QR algorithm. */
    reigen_hqr_f(r,C,0,0,N-1,N);
  } else {
    r[0] = _makecomplexf(-F[0],0.f);
  }
#endif
} /* rprootf() */

/* Return the scratch area size, in bytes. */
size_t rprootf_getScratchSize ( int N ) 
{ 
  NASSERT(N>0);
  return ( (N-1        +BBEN_2-1)/BBEN_2*BBEN_2*sz_f32 +  /* D: 1st subdiagonal of the companion matrix  */
           (N          +BBEN_2-1)/BBEN_2*BBEN_2*sz_f32 +  /* F: balanced polynomial coeffs               */
           (N*(N+3)/2-1+BBEN_2-1)/BBEN_2*BBEN_2*sz_f32 ); /* C: companion matrix, elements below the 1st *
                                                           *    subdiagonal are not stored               */
}

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
void balCompan( int * b, float32_t * d, float32_t * f, int N, int iterNum )
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
  const xtfloat * restrict dr;
        xtfloat * restrict dw;
  const xtfloat * restrict fr;
        xtfloat * restrict fw;
  float32_t f0,d0,d1;
  float32_t c,r,g,h,s,t;
  uint32_t ru,cu;
  uint32_t p,q,u,v;
  int32_t e;
  vbool1 blt, bconv;
  int iterCnt=0;
  int i,j;

  NASSERT(N>1);
  NASSERT(!(f[N-1]==0.f));

  for ( i=0; i<N-1; i++ ) d[i] = XT_CONST_S(1);

  do {
    bconv = BBE_MOVBA1(1);
    dr = dw = (xtfloat*)&d[0];
    fr = fw = (xtfloat*)&f[N-1];
    /*--------------------------------------------------------*
     * Process rows/cols 0                                    */
    xtfloat_loadip(d0, dr, sz_f32);
    xtfloat_loadip(f0, fr, -(int)sz_f32);
    r = XT_ABS_S(f0);
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
    BBE_MULF32T(f0,f0,h,blt);
    BBE_MULF32T(d0,d0,g,blt);
    xtfloat_storeip(f0, fw, -(int)sz_f32);
    /*--------------------------------------------------------*
     * Process rows/cols 1..N-2                               */
    for ( i=1; i<N-1; i++ ) {
      xtfloat_loadip(f0, fr, -(int)sz_f32);
      xtfloat_loadip(d1, dr, sz_f32);
      r = XT_ADD_S(XT_ABS_S(f0), d0);
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
      BBE_MULF32T(f0,f0,h,blt);
      BBE_MULF32T(d1,d1,g,blt);
      xtfloat_storeip(d0, dw, sz_f32);
      xtfloat_storeip(f0, fw, -(int)sz_f32);
      bconv = BBE_OPERATOR_ANDB1(bconv,BBE_NOTB1(blt));
      d0 = d1;
    } /* i */
    __Pragma("no_reorder");
    /*--------------------------------------------------------*
     * Process the last row/col                               */
    __Pragma("loop_count min=1");
    for ( c=0.f, j=1; j<N; j++ ) c += (float32_t)XT_ABS_S(*++fr);
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
      fw = (xtfloat*)fr;
      for ( j=1; j<N; j++ ) *fw-- = (float32_t)*fr--*g;
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
  NASSERT(!(f[N-1]==0.f));

  for ( i=0; i<N-1; i++ ) d[i] = 1.f;
  if (b) for ( i=0; i<N; i++ ) b[i] = 0;

  do {
    conv = 1;
    /*--------------------------------------------------------*
     * Process rows/cols 0..N-2                               */
    for ( i=0; i<N-1; i++ ) {
      r.f = fabsf(f[N-1-i]); if (i>0) r.f += d[i-1];
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
        f[N-1-i] *= h.f;
        d[i] *= g.f;
      }
    } /* i */
    /*--------------------------------------------------------*
     * Process the last row/col                               */
    r.f = d[N-2];
    for ( c.f=0.f, j=1; j<N; j++ ) c.f += fabsf(f[j]);
    NASSERT(!(r.f==0.f) && !(c.f==0.f));
    p = (int)(r.u>>23); u = ( p>0 ? (r.u&((1U<<23)-1)) : 0 );
    q = (int)(c.u>>23); v = ( q>0 ? (c.u&((1U<<23)-1)) : 0 );
    e = (p-q+(u>v))>>1;
    NASSERT((e>-127 && e<127) || !isfinite(r.f) || !isfinite(c.f));
    g.u = (uint32_t)(127+e)<<23;
    h.u = (uint32_t)(127-e)<<23;
    if (e>0) {
      g.u = (uint32_t)(e+127)<<23;
      h.f = 1.f/g.f;
    } else {
      h.u = (uint32_t)(127-e)<<23;
      g.f = 1.f/h.f;
    }
    if (r.f*h.f+c.f*g.f < 0.95f*(r.f+c.f)) {
      conv = 0;
      if (b) b[N-1] += e;
      d[N-2] *= h.f;
      for ( j=1; j<N; j++ ) f[j] *= g.f;
    }
  } while (++iterCnt<iterNum && !conv);
#endif
} /* balCompan() */

#else /* HAVE_VFPU */
DISCARD_FUN( void, rprootf, ( void * pScr, 
                              complex_float * restrict r,
                        const float32_t     * restrict c,
                        int N ) )

size_t rprootf_getScratchSize ( int N ) { return (0); }

#endif /* HAVE_VFPU */

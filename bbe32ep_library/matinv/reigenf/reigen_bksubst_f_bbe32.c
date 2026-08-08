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
    Backsubstitution to determine eigenvectors of a quasi-triangular real form.
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <math.h>
#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#include "vfpu_math.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define sz_f32    sizeof(float32_t)

#define EPS       FLT_EPSILON

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )

#define MAX(a,b)    ( (a)>(b) ? (a) : (b) )

#if 0
static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
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
#endif

/* 
 * Determine eigenvectors of a quasi-triangular real form T by in-place 
 * backsubstitution process. Real part of a complex eigenvector is stored
 * in the first vector of the conjugate pair, and the imaginary part is
 * stored in the second vector. 
 * Input:
 *   N               Matrix size
 *   e[N]            Eigenvalues, either real or in conjugate pairs
 * Input/Output:
 *   T[N*(N+3)/2-1]  NxN quasi-triangular form T with 2x2 blocks on the main
 *                   diagonal (in); N column eigenvectors (out). Zeros below
 *                   the first subdiagonal are not stored
 */

void reigen_bksubst_f ( float32_t     * restrict T, 
                  const complex_float * restrict e,
                  int N )
#if 0
{
  float32_t small,big;
  float32_t a,b,p,q,r,ra,s=0.f,sa,t=0.f;
  float32_t vr,vi,w,wr,wi,x,y,z=0.f;
  complex_float c;
  int en,na,i,j,m,exp;

  small = FLT_MIN;
  big = 1e3f;
  /* Process eigenpairs from the last to the first. */
  for ( en=N-1; en>=0; en-- ) {
    na = en-1;
    p = crealf(e[en]); q = cimagf(e[en]);
    if (q==0.f) {
      /* Real eigenvector */
      m = en; T[HIDX(en,en)] = 1.f;
      if (en<N-1) T[HIDX(en+1,en)] = 0.f;
      for ( i=na; i>=0; i-- ) {
        w = T[HIDX(i,i)]-p;
        for ( r=0.f, j=m; j<=en; j++ ) {
          r += T[HIDX(i,j)]*T[HIDX(j,en)];
        }
        wr = crealf(e[i]); wi = cimagf(e[i]);
        if (wi>0.f) {
          m = i;
          /* Solve a 2x2 system for a,b: [w x;y z]*[a,b] == [-r;-s]. */
          x = T[HIDX(i,i+1)]; y = T[HIDX(i+1,i)];
          q = (wr-p)*(wr-p) + wi*wi; a = (x*s-z*r)/q;
          b = ( fabsf(x)>fabsf(z) ? (-r-w*a)/x : (-s-y*a)/z );
          T[HIDX(i,en)] = a; T[HIDX(i+1,en)] = b;
          /* Conditionally rescale the vector to prevent overflow. */
          t = MAX(fabsf(a), fabsf(b));
        } else if (wi==0.f) {
          m = i;
          T[HIDX(i,en)] = ( w!=0.f ? -r/w : -r/MAX(EPS*fabsf(wr),small) );
          t = fabsf( T[HIDX(i,en)] );
        } else {
          /* Carry data to the next iteration of the i-loop. */
          z = w; s = r; continue;
        }
        /* Conditionally rescale the vector to prevent overflow. */
        if (t>=big) {
          exp = ilogbf(t);
          for ( j=i; j<=en; j++ ) {
            T[HIDX(j,en)] = ldexpf(T[HIDX(j,en)], -exp);
          }
        }
      } /* i */
    }
    else if (q<0.f) {
      /* Complex vector associated with p-1j*q */
      m = na;
      if ( fabsf(T[HIDX(en,na)]) > fabsf(T[HIDX(na,en)]) ) {
        T[HIDX(na,na)] = -(T[HIDX(en,en)]-p)/T[HIDX(en,na)];
        T[HIDX(na,en)] = -q/T[HIDX(en,na)];
      } else {
        c = _cdivf( _makecomplexf( -T[HIDX(na,en)], 0.f ), 
                    _makecomplexf( T[HIDX(na,na)]-p, q ) );
        T[HIDX(na,na)] = crealf(c); T[HIDX(na,en)] = cimagf(c);
      }
      T[HIDX(en,na)] = 1.f; T[HIDX(en,en)] = 0.f;
      if (en<N-1) T[HIDX(en+1,en)] = 0.f;
      for ( i=na-1; i>=0; i-- ) {
        wr = crealf(e[i]); wi = cimagf(e[i]); w = T[HIDX(i,i)]-p; 
        ra = T[HIDX(i,en)] * T[HIDX(en,na)]; sa = 0.f;
        for ( j=m; j<=na; j++ ) {
          ra += T[HIDX(i,j)] * T[HIDX(j,na)];
          sa += T[HIDX(i,j)] * T[HIDX(j,en)];
        }
        if (wi>0.f) {
          m = i;
          /* Solve a complex-valued 2x2 system for a+1j*b and c+1j*d:
           * [w+1j*q,x;y,z+1j*q]*[a+1j*b;c+1j*d] == [-ra-1j*sa;-r-1j*s] */
          x = T[HIDX(i,i+1)]; y = T[HIDX(i+1,i)];
          vr = (wr-p)*(wr-p) + wi*wi - q*q;
          vi = (wr-p)*2*q;
          if ( vr==0.f && vi==0.f ) {
            vr = MAX( small, EPS*(fabsf(wr)+fabsf(wi)) );
          }
          c =_cdivf( _makecomplexf( x*r-z*ra+q*sa, x*s-z*sa-q*ra ), _makecomplexf(vr,vi) );
          T[HIDX(i,na)] = crealf(c); T[HIDX(i,en)] = cimagf(c);
          if ( fabsf(x) > (fabsf(z)+fabsf(q)) ) {
            T[HIDX(i+1,na)] = (-ra-w*T[HIDX(i,na)]+q*T[HIDX(i,en)])/x;
            T[HIDX(i+1,en)] = (-sa-w*T[HIDX(i,en)]-q*T[HIDX(i,na)])/x;
          } else {
            c = _cdivf( _makecomplexf( -r-y*T[HIDX(i,na)], -s-y*T[HIDX(i,en)] ), 
                        _makecomplexf(z,q) );
            T[HIDX(i+1,na)] = crealf(c); T[HIDX(i+1,en)] = cimagf(c);
          }
          t = MAX( MAX(fabsf(T[HIDX(i,na)]), fabsf(T[HIDX(i+1,na)])),
                   MAX(fabsf(T[HIDX(i,en)]), fabsf(T[HIDX(i+1,en)])) );
        } else if (wi==0.f) {
          m = i;
          c = _cdivf( _makecomplexf(-ra,-sa), _makecomplexf(w,q) );
          T[HIDX(i,na)] = crealf(c); T[HIDX(i,en)] = cimagf(c);
          t = MAX( fabsf(T[HIDX(i,na)]), fabsf(T[HIDX(i,en)]) );
        } else {
          /* Carry data to the next iteration of the i-loop. */
          z = w; r = ra; s = sa; continue;
        }
        /* Conditionally rescale the vector to prevent overflow. */
        if (t>=big) {
          exp = ilogbf(t);
          for ( j=i; j<=en; j++ ) {
            T[HIDX(j,na)] = ldexpf(T[HIDX(j,na)], -exp);
            T[HIDX(j,en)] = ldexpf(T[HIDX(j,en)], -exp);
          }
        }
      } /* i */
    } /* q */
  } /* en */

} /* reigen_bksubst_f() */
#else
{
    const xb_vecNx16 exp_mask = BBE_MOVVINX16(BBE_FLOAT_MASK_EXP);

    xtcomplexfloat * restrict E = (xtcomplexfloat*)e;
    const xtfloat  * restrict T_r0;
    const xtfloat  * restrict T_r1;
          xtfloat  * restrict T_w;

    float32_t small, big;
    float32_t a, b, p, q, r = 0.f, ra, s = 0.f, sa, t = 0.f;
    float32_t vr, vi, w, wr, wi, x, y, z = 0.f;
    float32_t d1, d2, d3, r0, r1;
    xtcomplexfloat _c, c0, c1;
    vbool1 bless;
    int en, na, i, j, m, step;

    small = FLT_MIN;
    big = 1e3f;
    /* Process eigenpairs from the last to the first. */
    for (en = N - 1; en >= 0; en--) {
        na = en - 1;
        p = BBE_CREALCF32(E[en]); q = BBE_CIMAGCF32(E[en]);
        if (q == 0.f) {
            /* Real eigenvector */
            m = en; T[HIDX(en, en)] = XT_CONST_S(1);
            if (en<N - 1) T[HIDX(en + 1, en)] = XT_CONST_S(0);
            for (i = na; i >= 0; i--) {
                w = XT_SUB_S(T[HIDX(i, i)], p);

                T_r0 = (xtfloat *)&T[HIDX(m, en)];
                T_r1 = (xtfloat *)&T[HIDX(i, m)];
                step = N - m;
                r = XT_CONST_S(0);
                for (j = m; j <= en; j++) {
                    xtfloat_loadxp(a, T_r0, step--*sz_f32);
                    xtfloat_loadip(b, T_r1, sz_f32);
                    XT_MADD_S(r, b, a);
                }

                wr = BBE_CREALCF32(E[i]); wi = BBE_CIMAGCF32(E[i]);
                if (wi>0.f) {
                    m = i;
                    /* Solve a 2x2 system for a,b: [w x;y z]*[a,b] == [-r;-s]. */
                    x = T[HIDX(i, i + 1)]; y = T[HIDX(i + 1, i)];
                    q = XT_MUL_S( XT_SUB_S(wr, p), XT_SUB_S(wr,p));  XT_MADD_S(q, wi, wi);
                    a = XT_MUL_S(x, s);  XT_MSUB_S(a, z, r);
                    a = IT_FDIVF32(a,q,1);

                    bless = XT_OLT_S(XT_ABS_S(z), XT_ABS_S(x));
                    d1 = r; XT_MOVF_S(d1, s, bless);
                    d2 = w; XT_MOVF_S(d2, y, bless);
                    d3 = x; XT_MOVF_S(d3, z, bless);
                    b = XT_NEG_S(d1); XT_MSUB_S(b, d2, a);
                    b = IT_FDIVF32(b, d3, 1);

                    T[HIDX(i, en)] = a; T[HIDX(i + 1, en)] = b;
                    /* Conditionally rescale the vector to prevent overflow. */
                    t = XT_MAX_S(XT_ABS_S(a), XT_ABS_S(b));
                } else if (wi == 0.f) {
                    m = i;
                    a = XT_NEG_S(r);
                    if ( BBE_MOVAB1(XT_OEQ_S(w, 0.f)) ){
                        b = IT_FDIVF32(a, XT_MAX_S(XT_MUL_S(EPS, XT_ABS_S(wr)), small), 0);
                    } else{
                        b = IT_FDIVF32(a, w, 1);
                    }
                    T[HIDX(i, en)] = b;
                    t = XT_ABS_S(T[HIDX(i, en)]);
                } else {
                    /* Carry data to the next iteration of the i-loop. */
                    z = w; s = r; continue;
                }
                /* Conditionally rescale the vector to prevent overflow. */
                if (t >= big) {
                    float32_t exp;
                    step = N - i;
                    xb_vecNx16 _t = BBE_MOVNX16_FROMF32(t);
                    exp = BBE_MOVF32_FROMN_2XF32(BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMNX16(BBE_ANDNX16(_t, exp_mask))));
                    T_r0 = T_w = (xtfloat *)&T[HIDX(i, en)];
                    for (j = i; j <= en; j++) {
                        xtfloat_loadxp(a, T_r0, step*sz_f32);
                        xtfloat_storexp(XT_MUL_S(a, exp), T_w, step--*sz_f32);
                    }
                }
            } /* i */
        }
        else if (q<0.f) {
            /* Complex vector associated with p-1j*q */
            m = na;
            a = XT_ABS_S(T[HIDX(na, en)]); b = XT_ABS_S(T[HIDX(en, na)]);
            if ( BBE_MOVAB1(XT_OLT_S(a, b)) ) {
                T[HIDX(na, na)] = IT_FDIVF32(XT_NEG_S(XT_SUB_S(T[HIDX(en, en)], p)), T[HIDX(en, na)], 1);
                T[HIDX(na, en)] = IT_FDIVF32(XT_NEG_S(q), T[HIDX(en, na)], 1);
            } else {
                c0 = BBE_CMPLXF32(0.f, XT_NEG_S(T[HIDX(na, en)]));
                c1 = BBE_CMPLXF32(q, XT_SUB_S(T[HIDX(na, na)], p));
                _c = IT_DIVCF32(c0,c1, 1);                
                T[HIDX(na, na)] = BBE_CREALCF32(_c); T[HIDX(na, en)] = BBE_CIMAGCF32(_c);
            }
            T[HIDX(en, na)] = XT_CONST_S(1); T[HIDX(en, en)] = XT_CONST_S(0);
            if (en<N - 1) T[HIDX(en + 1, en)] = XT_CONST_S(0);
            for (i = na - 1; i >= 0; i--) {
                wr = BBE_CREALCF32(E[i]); wi = BBE_CIMAGCF32(E[i]);
                w = XT_SUB_S(T[HIDX(i, i)], p);
                ra = XT_MUL_S(T[HIDX(i, en)], T[HIDX(en, na)]); sa = XT_CONST_S(0);

                T_r0 = (xtfloat *)&T[HIDX(m, na)];
                T_r1 = (xtfloat *)&T[HIDX(i, m)];
                step = N - m;
                for (j = m; j <= na; j++) {
                    xtfloat_loadip(a, T_r1, sz_f32);
                    b = xtfloat_loadi(T_r0, sz_f32); // +1
                    XT_MADD_S(sa, a, b);
                    xtfloat_loadxp(b, T_r0, step--*sz_f32); // +0
                    XT_MADD_S(ra, a, b);
                }

                if (wi>0.f) {
                    m = i;
                    /* Solve a complex-valued 2x2 system for a+1j*b and c+1j*d:
                    * [w+1j*q,x;y,z+1j*q]*[a+1j*b;c+1j*d] == [-ra-1j*sa;-r-1j*s] */
                    x = T[HIDX(i, i + 1)]; y = T[HIDX(i + 1, i)];
                    vr = XT_MUL_S(XT_SUB_S(wr, p), XT_SUB_S(wr, p));
                    XT_MADD_S(vr, wi, wi); XT_MSUB_S(vr, q, q);
                    vi = XT_MUL_S(XT_MUL_S(XT_SUB_S(wr, p), 2), q);
                    if (vr == 0.f && vi == 0.f) {
                        vr = XT_MAX_S(small, XT_MUL_S(EPS, XT_ADD_S(XT_ABS_S(wr), XT_ABS_S(wi)) ));
                    }
                    a = XT_MUL_S(x, s); XT_MSUB_S(a, z, sa); XT_MSUB_S(a, q, ra);
                    b = XT_MUL_S(x, r); XT_MSUB_S(b, z, ra); XT_MADD_S(b, q, sa);
                    c0 = BBE_CMPLXF32(a, b);
                    c1 = BBE_CMPLXF32(vi, vr);
                    _c = IT_DIVCF32(c0, c1, 1);
                    T[HIDX(i, na)] = BBE_CREALCF32(_c); T[HIDX(i, en)] = BBE_CIMAGCF32(_c);
                    a = XT_ABS_S(x); b = XT_ADD_S(XT_ABS_S(z), XT_ABS_S(q));
                    if ( BBE_MOVAB1(XT_OLT_S(b, a)) ) {
                        a = XT_NEG_S(ra); XT_MSUB_S(a, w, T[HIDX(i, na)]); XT_MADD_S(a, q, T[HIDX(i, en)]);
                        b = XT_NEG_S(sa); XT_MSUB_S(b, w, T[HIDX(i, en)]); XT_MSUB_S(b, q, T[HIDX(i, na)]);
                        // div a,b by x
                        vbool1 blt = XT_OLT_S(XT_ABS_S(x), FLT_MIN);
                        BBE_MULF32T(x, x, 8388608.f, blt);
                        d2 = XT_RECIP0_S(x);
                        d1 = XT_CONST_S(1); XT_MSUBN_S(d1, d2, x); XT_MADDN_S(d2, d1, d2);
                        r0 = XT_MUL_S(a, d2); r1 = XT_MUL_S(b, d2);
                        XT_MSUB_S(a, x, r0); XT_MADD_S(r0, d2, a);
                        XT_MSUB_S(b, x, r1); XT_MADD_S(r1, d2, b);
                        BBE_MULF32T(r0, r0, 8388608.f, blt);
                        BBE_MULF32T(r1, r1, 8388608.f, blt);
                        T[HIDX(i + 1, na)] = r0;
                        T[HIDX(i + 1, en)] = r1;
                    } else {
                        a = XT_NEG_S(s); XT_MSUB_S(a, y, T[HIDX(i, en)]);
                        b = XT_NEG_S(r); XT_MSUB_S(b, y, T[HIDX(i, na)]);
                        c0 = BBE_CMPLXF32(a, b);
                        c1 = BBE_CMPLXF32(q,z);
                        _c = IT_DIVCF32(c0, c1,1);
                        T[HIDX(i + 1, na)] = BBE_CREALCF32(_c); T[HIDX(i + 1, en)] = BBE_CIMAGCF32(_c);
                    }
                    t = XT_MAX_S(XT_MAX_S(XT_ABS_S(T[HIDX(i, na)]), XT_ABS_S(T[HIDX(i + 1, na)])),
                        XT_MAX_S(XT_ABS_S(T[HIDX(i, en)]), XT_ABS_S(T[HIDX(i + 1, en)])));
                } else if (wi == 0.f) {
                    m = i;
                    c0 = BBE_CMPLXF32(sa, ra); c0 = BBE_NEGCF32(c0);
                    c1 = BBE_CMPLXF32(q,w);
                    _c = IT_DIVCF32(c0, c1, 1);
                    T[HIDX(i, na)] = BBE_CREALCF32(_c); T[HIDX(i, en)] = BBE_CIMAGCF32(_c);
                    t = XT_MAX_S(XT_ABS_S(T[HIDX(i, na)]), XT_ABS_S(T[HIDX(i, en)]));
                } else {
                    /* Carry data to the next iteration of the i-loop. */
                    z = w; r = ra; s = sa; continue;
                }
                /* Conditionally rescale the vector to prevent overflow. */
                if (t >= big) {
                    float32_t exp;
                    step = N - i;
                    xb_vecNx16 _t = BBE_MOVNX16_FROMF32(t);
                    exp = BBE_MOVF32_FROMN_2XF32(BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMNX16(BBE_ANDNX16(_t, exp_mask))));
                    T_r0 = T_w = (xtfloat *)&T[HIDX(i, na)];
                    for (j = i; j <= en; j++) {
                        b = xtfloat_loadi(T_r0, sz_f32); // +1
                        xtfloat_loadxp(a, T_r0, step*sz_f32); // +0
                        xtfloat_storei(XT_MUL_S(b, exp), T_w, sz_f32);
                        xtfloat_storexp(XT_MUL_S(a, exp), T_w, step--*sz_f32);
                    }
                }
            } /* i */
        } /* q */
    } /* en */

} /* reigen_bksubst_f() */
#endif

#endif /* HAVE_VFPU */

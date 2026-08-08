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
    Backsubstitution to determine eigenvectors of a triangular complex form.
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

#if HAVE_VFPU

#define EPS       FLT_EPSILON

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
* in compact packed format. Compactness implies that zeros below the
* first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )

#define sz_f32c  sizeof(complex_float)

#if 0
#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#define MAX(a,b)    ( (a)>(b) ? (a) : (b) )

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

/* Multiply a complex floating-point number by an integral power
 * of two, single precision. */
static complex_float cldexpf( complex_float x, int y )
{
  return ( _makecomplexf(ldexpf(crealf(x),y), ldexpf(cimagf(x),y)) );
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
 * Determine eigenvectors of a triangular complex form T by in-place 
 * backsubstitution process. 
 * Input:
 *   N               Matrix size
 * Input/Output:
 *   T[N*(N+3)/2-1]  NxN triangular form T (in); N column eigenvectors (out).
 *                   Zeros below the first subdiagonal are not stored. Values
 *                   on the subdiagonal are not specified for both input and
 *                   output.
 */

void eigen_bksubst_f ( complex_float * restrict T, int N )
#if 0
{
  float32_t small,big,w;
  complex_float x,y,z;
  int en,i,j,exp;

  const complex_float c0f = _makecomplexf(0.f,0.f);
  const complex_float c1f = _makecomplexf(1.f,0.f);
  
  small = FLT_MIN;
  big = 1e3f;
  /* Process eigenpairs from the last to the first. */
  for ( en=N-1; en>=0; en-- ) {
    x = T[HIDX(en,en)];
    T[HIDX(en,en)] = c1f;
    for ( i=en-1; i>=0; i-- ) {
      y = csubf(x, T[HIDX(i,i)]);
      if (cabs2f(y)==0) {
        y = _makecomplexf( MAX(EPS*cmagf(T[HIDX(i,i)]), small), 0.f );
      }
      for ( z=c0f, j=i+1; j<=en; j++ ) {
        z = caddf(z, cmulf(T[HIDX(i,j)], T[HIDX(j,en)]));
      }
      T[HIDX(i,en)] = _cdivf(z,y);
      /* Conditionally rescale the vector to prevent overflow. */
      w = cmagf(T[HIDX(i,en)]);
      if (w>=big) {
        exp = ilogbf(w);
        for ( j=i; j<=en; j++ ) {
          T[HIDX(j,en)] = cldexpf(T[HIDX(j,en)], -exp);
        }
      }
    }
  } /* en */

} /* eigen_bksubst_f() */
#else
{
    const xtcomplexfloat * restrict T_r0;
    const xtcomplexfloat * restrict T_r1;
          xtcomplexfloat * restrict _T = (xtcomplexfloat *)T;
          xtcomplexfloat * restrict T_w;
    const float32_t small = FLT_MIN;
    const float32_t big = 1e3f;
    const xb_vecNx16 exp_mask = BBE_MOVVINX16(BBE_FLOAT_MASK_EXP);
    const xtcomplexfloat _c0f = BBE_CONSTCF32(0);
    const xtcomplexfloat _c1f = BBE_CMPLXF32(XT_CONST_S(0), XT_CONST_S(1));

    xtcomplexfloat _x, _y, _z, z0, z1, z2, z3;
    float32_t  w;

    int en, i, j, step;

    /* Process eigenpairs from the last to the first. */
    for (en = N - 1; en >= 0; en--) {
        _x = _T[HIDX(en, en)];
        _T[HIDX(en, en)] = _c1f;
        for (i = en - 1; i >= 0; i--) {
            _y = BBE_OPERATOR_SUBCF32(_x, _T[HIDX(i, i)]);
            if (BBE_MOVAB1(XT_OEQ_S(IT_ABS2CF32(_y), 0.f))){
                _y = BBE_CMPLXF32(XT_CONST_S(0), XT_MAX_S(XT_MUL_S(EPS, IT_CMAGCF32(_T[HIDX(i, i)])), small));
            }
            T_r0 = &_T[HIDX(i, i + 1)];
            T_r1 = &_T[HIDX(i+1, en)];
            step = N - (i+1);
            for (z0 = z1 = z2 = z3 = _c0f, j = 0; j < (en - i)>>1; j++) {
                xtcomplexfloat t0, t1;
                xtcomplexfloat_loadip(t0, T_r0, sz_f32c);
                xtcomplexfloat_loadxp(t1, T_r1, step--*sz_f32c);
                BBE_MULMASCF32(z0, t0, t1, 0, 4);
                BBE_MULMASCF32(z1, t0, t1, 1, 11);
                xtcomplexfloat_loadip(t0, T_r0, sz_f32c);
                xtcomplexfloat_loadxp(t1, T_r1, step--*sz_f32c);
                BBE_MULMASCF32(z2, t0, t1, 0, 4);
                BBE_MULMASCF32(z3, t0, t1, 1, 11);
            }
            if ((en - i) & 1){
                xtcomplexfloat t0, t1;
                xtcomplexfloat_loadip(t0, T_r0, sz_f32c);
                xtcomplexfloat_loadxp(t1, T_r1, step--*sz_f32c);
                BBE_MULMASCF32(z0, t0, t1, 0, 4);
                BBE_MULMASCF32(z1, t0, t1, 1, 11);
            }
            _z = BBE_OPERATOR_ADDCF32(BBE_OPERATOR_ADDCF32(z0, z1), BBE_OPERATOR_ADDCF32(z2,z3) );
            _T[HIDX(i, en)] = IT_DIVCF32(_z, _y, 1);
            /* Conditionally rescale the vector to prevent overflow. */
            w = IT_CMAGCF32(_T[HIDX(i, en)]);
            if (w >= big) {
                float32_t exp;
                xb_vecNx16 t = BBE_MOVNX16_FROMF32(w);
                t = BBE_ANDNX16(t, exp_mask);
                exp = XT_RECIP0_S(BBE_MOVF32_FROMN_2XF32(BBE_MOVN_2XF32_FROMNX16(t)));
                step = N - i;
                T_r0 = T_w = &_T[HIDX(i, en)];
                for (j = i; j <= en; j++) {
                    xtcomplexfloat t;
                    xtcomplexfloat_loadxp(t, T_r0, step*sz_f32c);
                    t = IT_RCMULCF32(exp, t);
                    xtcomplexfloat_storexp(t, T_w, step--*sz_f32c);
                }
                __Pragma("no_reorder");
            }
        }
    } /* en */

} /* eigen_bksubst_f() */
#endif

#endif /* HAVE_VFPU */

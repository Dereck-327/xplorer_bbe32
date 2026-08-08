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
    NatureDSP_Baseband library.
    Real and complex arithmetic primitives optimized for BBEN VFPU
    IntegrIT, 2006-2017
*/

#ifndef __VFPU_MATH_H
#define __VFPU_MATH_H

#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"

#if HAVE_VFPU

/* Vector real by complex multiply. */
#define IT_RCMULN_4XCF32(x,y)   BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MULN_2XF32((x), BBE_MOVN_2XF32_FROMN_4XCF32(y)))
/* Sclalar real by complex multiply */
#define IT_RCMULCF32(a,b)       BBE_MOVCF32_FROMN_2XF32(BBE_MULMN_2XF32(BBE_MOVN_2XF32_FROMF32(a),BBE_MOVN_2XF32_FROMCF32(b),0,4))

/* Sclalar real by complex multiply-add */
#define IT_RCMULACF32(a,b,c)    __IT_RCMULACF32(&(a),&(b),&(c))
ATTRIBUTE_ALWAYS_INLINE
inline_ void __IT_RCMULACF32( void * pa, const void * pb, const void * pc )
{
  xb_vecN_2xf32 a,b,c;
  a = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pa);
  b = BBE_MOVN_2XF32_FROMF32(*(xtfloat*)pb);
  c = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pc);
  BBE_MULMASN_2XF32(a,b,c,0,4);
  *(xtcomplexfloat*)pa = BBE_MOVCF32_FROMN_2XF32(a);

} /* __IT_RCMULACF32() */

/* Sclalar real by complex multiply-subtract */
#define IT_RCMULSCF32(a,b,c)    __IT_RCMULSCF32(&(a),&(b),&(c))
ATTRIBUTE_ALWAYS_INLINE
inline_ void __IT_RCMULSCF32( void * pa, const void * pb, const void * pc )
{
  xb_vecN_2xf32 a,b,c;
  a = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pa);
  b = BBE_MOVN_2XF32_FROMF32(*(xtfloat*)pb);
  c = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pc);
  BBE_MULMASN_2XF32(a,b,c,3,4);
  *(xtcomplexfloat*)pa = BBE_MOVCF32_FROMN_2XF32(a);

} /* __IT_RCMULSCF32() */

/* Sclalar real by complex multiply, predicated */
#define IT_RCMULCF32T(a,b,c,d)    __IT_RCMULCF32T(&(a),&(b),&(c),&(d))
ATTRIBUTE_ALWAYS_INLINE
inline_ void __IT_RCMULCF32T( void * pa, const void * pb, const void * pc, const void * pd )
{
  xb_vecN_2xf32 a,b,c;
  vboolN_2 d;
  a = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pa);
  b = BBE_MOVN_2XF32_FROMF32(*(xtfloat*)pb);
  c = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)pc);
  d = BBE_MOVN_2_FROMN(BBE_EXT0IB(BBE_MOVN_FROM1(*(vbool1*)pd),4));
  BBE_MULMN_2XF32T(a,b,c,0,4,d);
  *(xtcomplexfloat*)pa = BBE_MOVCF32_FROMN_2XF32(a);
  
} /* __IT_RCMULCF32T() */

/* Complex scalar division. If real and imagonary components of the
 * divisor (second argument y) may appeat subnormal simultaneously,
 * then immediate input argumemnt scaleDivsr should be non-zero to
 * prevent overflow. Note that accuracy still degrades if divisor is
 * subnormal and dividend is close to subnormal). */
#define IT_DIVCF32(x,y,scaleDivsr)   __IT_DIVCF32(&(x),&(y),(scaleDivsr))
ATTRIBUTE_ALWAYS_INLINE
inline_ xtcomplexfloat __IT_DIVCF32( const void * px, const void * py, int scaleDivsr )
{
  /*
   * Based on cdiv code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   *
   * static complex_float _cdivf( complex_float x, complex_float y )
   * {
   *   float32_t xr,xi,yr,yi,zr,zi,h;
   *   xr = crealf(x); xi = cimagf(x);
   *   yr = crealf(y); yi = cimagf(y);
   *   if (fabsf(yr)>fabsf(yi)) {
   *     h = yi/yr; yr = h*yi+yr;
   *     zr = (xr+h*xi)/yr;
   *     zi = (xi-h*xr)/yr;
   *   } else {
   *     h = yr/yi; yi = h*yr+yi;
   *     zr = (h*xr+xi)/yi;
   *     zi = (h*xi-xr)/yi;
   *   }
   *   return ( _makecomplexf(zr,zi) );
   * }
   */
#if 1
  xb_vecN_2xf32 x0,x1,y0,y1,w0,w1,z0,z1;
  xb_vecN_2xf32 f,g,h;
  vboolN_2 bgt,blt,blt0,blt1;
  x0 = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)px);
  x1 = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x0,BBE_SHFLI_SWAP_2));
  y0 = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)py);
  y1 = BBE_SHFLN_2XF32I(y0,BBE_SHFLI_SWAP_2);

  w0 = BBE_SHFLN_2XF32I(BBE_ABSN_2XF32(y0), BBE_SHFLI_DUPLICATE_2_EVEN);
  w1 = BBE_SHFLN_2XF32I(BBE_ABSN_2XF32(y0), BBE_SHFLI_DUPLICATE_2_ODD);
  bgt = BBE_OGTN_2XF32(w0,w1);

  if (scaleDivsr) {
    blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(y0), FLT_MIN);
    blt0 = BBE_OLTN_2XF32(w0, FLT_MIN);
    blt1 = BBE_OLTN_2XF32(w1, FLT_MIN);
    BBE_MULN_2XF32T(y0,8388608.f,y0,blt);
    BBE_MULN_2XF32T(y1,8388608.f,y1,blt);
  }

  g = BBE_RECIP0N_2XF32(y0);
  f = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(f,g,y0);
  BBE_MULAN_2XF32(g,f,g);
  f = y1; h = BBE_MULN_2XF32(f,g);
  BBE_MULSN_2XF32(f,h,y0);
  BBE_MULAN_2XF32(h,f,g);
  w0 = x0; BBE_MULMASN_2XF32(w0,x0,h,2,2);
  w1 = x1; BBE_MULMASN_2XF32(w1,x1,h,1,7);

  if (scaleDivsr) {
    BBE_MULN_2XF32T(w0,8388608.f,w0,blt0);
    BBE_MULN_2XF32T(w1,8388608.f,w1,blt1);
  }

  f = y0; BBE_MULMASN_2XF32(f,h,y1,0,12);
  g = BBE_RECIP0N_2XF32(f);
  h = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(h,f,g);
  BBE_MULAN_2XF32(g,h,g);

  z0 = BBE_MULMN_2XF32(w0,g,0,8);
  BBE_MULMASN_2XF32(w0,z0,f,3,8);
  BBE_MULMASN_2XF32(z0,w0,g,0,8);

  z1 = BBE_MULMN_2XF32(w1,g,0,13);
  BBE_MULMASN_2XF32(w1,z1,f,3,13);
  BBE_MULMASN_2XF32(z1,w1,g,0,13);

  return BBE_MOVCF32_FROMN_2XF32(BBE_MOVN_2XF32T(z0,z1,bgt));
#else
  xb_vecN_2xf32 x0,x1,y0,y1,w0,w1,z0,z1;
  xb_vecN_2xf32 f,g,h;
  vboolN_2 bgt;
  x0 = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)px);
  x1 = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x0,BBE_SHFLI_SWAP_2));
  y0 = BBE_MOVN_2XF32_FROMCF32(*(xtcomplexfloat*)py);
  y1 = BBE_SHFLN_2XF32I(y0,BBE_SHFLI_SWAP_2);
  bgt = BBE_OGTN_2XF32(BBE_ABSN_2XF32(y0),BBE_ABSN_2XF32(y1));
  bgt = BBE_MOVN_2_FROMN(BBE_EXT0IB(BBE_MOVN_FROMN_2(bgt),4));

  if (scaleDivsr) {
    /* This variant correctly handles subnormals in the divisor. */
    vboolN_2 blt;
    xb_vecN_2xf32 s = y0;
    blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(s), FLT_MIN);
    BBE_MULN_2XF32T(s,s,8388608.f,blt);
    g = BBE_RECIP0N_2XF32(s);
    f = BBE_CONSTN_2XF32(1);
    BBE_MULSN_2XF32(f,g,s);
    BBE_MULAN_2XF32(g,f,g);
    f = y1; h = BBE_MULN_2XF32(f,g);
    BBE_MULSN_2XF32(f,h,s);
    BBE_MULAN_2XF32(h,f,g);
    BBE_MULN_2XF32T(h,h,8388608.f,blt);
  } else {
    /* This variant overflows when real and imaginary parts 
     * of the divisor are subnormal. */
    g = BBE_RECIP0N_2XF32(y0);
    f = BBE_CONSTN_2XF32(1);
    BBE_MULSN_2XF32(f,g,y0);
    BBE_MULAN_2XF32(g,f,g);
    f = y1; h = BBE_MULN_2XF32(f,g);
    BBE_MULSN_2XF32(f,h,y0);
    BBE_MULAN_2XF32(h,f,g);
  }

  f = y0; BBE_MULMASN_2XF32(f,h,y1,0,12);
  w0 = x0; BBE_MULMASN_2XF32(w0,x0,h,2,2);
  w1 = x1; BBE_MULMASN_2XF32(w1,x1,h,1,7);

  g = BBE_RECIP0N_2XF32(f);
  h = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(h,f,g);
  BBE_MULAN_2XF32(g,h,g);

  z0 = BBE_MULMN_2XF32(w0,g,0,8);
  BBE_MULMASN_2XF32(w0,z0,f,3,8);
  BBE_MULMASN_2XF32(z0,w0,g,0,8);

  z1 = BBE_MULMN_2XF32(w1,g,0,13);
  BBE_MULMASN_2XF32(w1,z1,f,3,13);
  BBE_MULMASN_2XF32(z1,w1,g,0,13);

  return BBE_MOVCF32_FROMN_2XF32(BBE_MOVN_2XF32T(z0,z1,bgt));
#endif
} /* __IT_DIVCF32() */

/* Complex vector division. If real and imagonary components of the
 * divisor (second argument y) may appeat subnormal simultaneously,
 * then immediate input argumemnt scaleDivsr should be non-zero to
 * prevent overflow. Note that accuracy still degrades if divisor is
 * subnormal and dividend is close to subnormal). */
#define IT_DIVN_4XCF32(x,y,scaleDivsr)   __IT_DIVN_4XCF32(&(x),&(y),(scaleDivsr))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_4xcf32 __IT_DIVN_4XCF32( const void * px, const void * py, int scaleDivsr )
{
  /*
   * Based on cdiv code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   *
   * static complex_float _cdivf( complex_float x, complex_float y )
   * {
   *   float32_t xr,xi,yr,yi,zr,zi,h;
   *   xr = crealf(x); xi = cimagf(x);
   *   yr = crealf(y); yi = cimagf(y);
   *   if (fabsf(yr)>fabsf(yi)) {
   *     h = yi/yr; yr = h*yi+yr;
   *     zr = (xr+h*xi)/yr;
   *     zi = (xi-h*xr)/yr;
   *   } else {
   *     h = yr/yi; yi = h*yr+yi;
   *     zr = (h*xr+xi)/yi;
   *     zi = (h*xi-xr)/yi;
   *   }
   *   return ( _makecomplexf(zr,zi) );
   * }
   */
#if 1
  xb_vecN_2xf32 x0,x1,y0,y1,w0,w1,z0,z1;
  xb_vecN_2xf32 f,g,h;
  vboolN_2 bgt,blt,blt0,blt1;
  x0 = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px);
  x1 = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x0,BBE_SHFLI_SWAP_2));
  y0 = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)py);
  y1 = BBE_SHFLN_2XF32I(y0,BBE_SHFLI_SWAP_2);

  w0 = BBE_SHFLN_2XF32I(BBE_ABSN_2XF32(y0), BBE_SHFLI_DUPLICATE_2_EVEN);
  w1 = BBE_SHFLN_2XF32I(BBE_ABSN_2XF32(y0), BBE_SHFLI_DUPLICATE_2_ODD);
  bgt = BBE_OGTN_2XF32(w0,w1);

  if (scaleDivsr) {
    blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(y0), FLT_MIN);
    blt0 = BBE_OLTN_2XF32(w0, FLT_MIN);
    blt1 = BBE_OLTN_2XF32(w1, FLT_MIN);
    BBE_MULN_2XF32T(y0,8388608.f,y0,blt);
    BBE_MULN_2XF32T(y1,8388608.f,y1,blt);
  }

  g = BBE_RECIP0N_2XF32(y0);
  f = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(f,g,y0);
  BBE_MULAN_2XF32(g,f,g);
  f = y1; h = BBE_MULN_2XF32(f,g);
  BBE_MULSN_2XF32(f,h,y0);
  BBE_MULAN_2XF32(h,f,g);
  w0 = x0; BBE_MULMASN_2XF32(w0,x0,h,2,2);
  w1 = x1; BBE_MULMASN_2XF32(w1,x1,h,1,7);

  if (scaleDivsr) {
    BBE_MULN_2XF32T(w0,8388608.f,w0,blt0);
    BBE_MULN_2XF32T(w1,8388608.f,w1,blt1);
  }

  f = y0; BBE_MULMASN_2XF32(f,h,y1,0,12);
  g = BBE_RECIP0N_2XF32(f);
  h = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(h,f,g);
  BBE_MULAN_2XF32(g,h,g);

  z0 = BBE_MULMN_2XF32(w0,g,0,8);
  BBE_MULMASN_2XF32(w0,z0,f,3,8);
  BBE_MULMASN_2XF32(z0,w0,g,0,8);

  z1 = BBE_MULMN_2XF32(w1,g,0,13);
  BBE_MULMASN_2XF32(w1,z1,f,3,13);
  BBE_MULMASN_2XF32(z1,w1,g,0,13);

  return BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MOVN_2XF32T(z0,z1,bgt));
#else
  xb_vecN_2xf32 x0,x1,y0,y1,w0,w1,z0,z1;
  xb_vecN_2xf32 f,g,h;
  vboolN_2 bgt;
  x0 = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px);
  x1 = BBE_CONJN_2XF32(BBE_SHFLN_2XF32I(x0,BBE_SHFLI_SWAP_2));
  y0 = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)py);
  y1 = BBE_SHFLN_2XF32I(y0,BBE_SHFLI_SWAP_2);

  w0 = BBE_SELN_2XF32I(y0,y0,BBE_SELI_INTERLEAVE_2_EVEN);
  w1 = BBE_SELN_2XF32I(y1,y1,BBE_SELI_INTERLEAVE_2_EVEN);
  bgt = BBE_OGTN_2XF32(BBE_ABSN_2XF32(w0),BBE_ABSN_2XF32(w1));

  if (scaleDivsr) {
    /* This variant correctly handles subnormals in the divisor. */
    vboolN_2 blt;
    xb_vecN_2xf32 s = y0;
    blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(s), FLT_MIN);
    BBE_MULN_2XF32T(s,s,8388608.f,blt);
    g = BBE_RECIP0N_2XF32(s);
    f = BBE_CONSTN_2XF32(1);
    BBE_MULSN_2XF32(f,g,s);
    BBE_MULAN_2XF32(g,f,g);
    f = y1; h = BBE_MULN_2XF32(f,g);
    BBE_MULSN_2XF32(f,h,s);
    BBE_MULAN_2XF32(h,f,g);
    BBE_MULN_2XF32T(h,h,8388608.f,blt);
  } else {
    /* This variant overflows when real and imaginary parts 
     * of the divisor are subnormal. */
    g = BBE_RECIP0N_2XF32(y0);
    f = BBE_CONSTN_2XF32(1);
    BBE_MULSN_2XF32(f,g,y0);
    BBE_MULAN_2XF32(g,f,g);
    f = y1; h = BBE_MULN_2XF32(f,g);
    BBE_MULSN_2XF32(f,h,y0);
    BBE_MULAN_2XF32(h,f,g);
  }

  f = y0; BBE_MULMASN_2XF32(f,h,y1,0,12);
  w0 = x0; BBE_MULMASN_2XF32(w0,x0,h,2,2);
  w1 = x1; BBE_MULMASN_2XF32(w1,x1,h,1,7);

  g = BBE_RECIP0N_2XF32(f);
  h = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(h,f,g);
  BBE_MULAN_2XF32(g,h,g);

  z0 = BBE_MULMN_2XF32(w0,g,0,8);
  BBE_MULMASN_2XF32(w0,z0,f,3,8);
  BBE_MULMASN_2XF32(z0,w0,g,0,8);

  z1 = BBE_MULMN_2XF32(w1,g,0,13);
  BBE_MULMASN_2XF32(w1,z1,f,3,13);
  BBE_MULMASN_2XF32(z1,w1,g,0,13);

  return BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MOVN_2XF32T(z0,z1,bgt));
#endif
} /* __IT_DIVCF32() */

/* Complex vector division for data stored in real vectors. If real
 * and imagonary components of the divisor (second argument b) may
 * appeat subnormal simultaneously, then immediate input argumemnt
 * scaleDivsr should be non-zero to prevent overflow. Note that 
 * accuracy still degrades if divisor is subnormal and dividend 
 * is close to subnormal). */
#define IT_CDIVN_2XF32(a,b,scaleDivsr)     __IT_CDIVN_2XF32(&(a),&(b),(scaleDivsr))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_2xf32 __IT_CDIVN_2XF32( const void * pa, const void * pb, int scaleDivsr )
{
  xb_vecN_4xcf32 a,b,c;
  a = BBE_MOVN_4XCF32_FROMN_2XF32(*(xb_vecN_2xf32*)pa);
  b = BBE_MOVN_4XCF32_FROMN_2XF32(*(xb_vecN_2xf32*)pb);
  c = IT_DIVN_4XCF32(a,b,scaleDivsr);
  return BBE_MOVN_2XF32_FROMN_4XCF32(c);
}

/* Scalar complex by real division */
#define IT_CRDIVCF32(x,y)   __IT_CRDIVCF32(&(x),(y))
ATTRIBUTE_ALWAYS_INLINE
inline_ xtcomplexfloat __IT_CRDIVCF32( const void * px, xtfloat sy )
{
  xb_vecN_2xf32 x,y,r,t;
  x = BBE_MOVN_2XF32_FROMN_4XCF32(*(xtcomplexfloat*)px);
  y = BBE_REPN_2XF32(BBE_MOVN_2XF32_FROMF32(sy),0);
  r = BBE_RECIP0N_2XF32(y);
  t = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(t,r,y);
  BBE_MULAN_2XF32(r,t,r);
  t = BBE_MULN_2XF32(r,x);
  BBE_MULSN_2XF32(x,y,t);
  BBE_MULAN_2XF32(t,r,x);
  return BBE_MOVCF32_FROMN_2XF32(t);

} /* __IT_CRDIVN_4XCF32() */

/* Complex by real vector division */
#define IT_CRDIVN_4XCF32(x,y)   __IT_CRDIVN_4XCF32(&(x),&(y))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_4xcf32 __IT_CRDIVN_4XCF32( const void * px, const void * py )
{
  xb_vecN_2xf32 x,y,r,t;
  x = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px);
  y = *(xb_vecN_2xf32*)py;
  r = BBE_RECIP0N_2XF32(y);
  t = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(t,r,y);
  BBE_MULAN_2XF32(r,t,r);
  t = BBE_MULN_2XF32(r,x);
  BBE_MULSN_2XF32(x,y,t);
  BBE_MULAN_2XF32(t,r,x);
  return BBE_MOVN_4XCF32_FROMN_2XF32(t);

} /* __IT_CRDIVN_4XCF32() */

/* Fast real scalar division based on reciprocal approximation. If
 * argument b may be subnormal, set scaleDivsr to non-zero to
 * prevent overflow (accuracy still degrades if divisor is 
 * subnormal and dividend is close to subnormal).
 * Accuracy is not worse than 1 ULP for normal arguments. */
ATTRIBUTE_ALWAYS_INLINE
inline_ xtfloat IT_FDIVF32( xtfloat a, xtfloat b, int scaleDivsr )
{
  float32_t r,t;
  /* t <- ~1/b */
  if (scaleDivsr) {
    /* This variant correctly handles subnormals in the divisor. */
    vbool1 blt = XT_OLT_S(XT_ABS_S(b), FLT_MIN);
    BBE_MULF32T(b,b,8388608.f,blt);
    t = XT_RECIP0_S(b);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = XT_CONST_S(1); XT_MSUBN_S(r,t,b); XT_MADDN_S(t,r,t);
    /* r <- ~a/b */
    r = XT_MUL_S(a,t);
    /* Modified Newton-Raphson iteration for a quotient */
    XT_MSUB_S(a,b,r); XT_MADD_S(r,t,a);
    BBE_MULF32T(r,r,8388608.f,blt);
  } else {
    /* This variant overflows when the divisor is subnormal. */
    t = XT_RECIP0_S(b);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = XT_CONST_S(1); XT_MSUBN_S(r,t,b); XT_MADDN_S(t,r,t);
    /* r <- ~a/b */
    r = XT_MUL_S(a,t);
    /* Modified Newton-Raphson iteration for a quotient */
    XT_MSUB_S(a,b,r); XT_MADD_S(r,t,a);
  }
  return (r);

} /* IT_FDIVF32() */

/* Fast real vector division based on reciprocal approximation. If
 * argument b may be subnormal, set scaleDivsr to non-zero to
 * prevent overflow (accuracy still degrades if divisor is 
 * subnormal and dividend is close to subnormal).
 * Accuracy is not worse than 1 ULP for normal arguments. */
#define IT_FDIVN_2XF32(a,b,scaleDivsr)     __IT_FDIVN_2XF32(&(a),&(b),(scaleDivsr))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_2xf32 __IT_FDIVN_2XF32( const void * pa, const void * pb, int scaleDivsr )
{
  xb_vecN_2xf32 a,b,r,t;
  a = *(xb_vecN_2xf32*)pa;
  b = *(xb_vecN_2xf32*)pb;
  /* t <- ~1/b */
  if (scaleDivsr) {
    /* This variant correctly handles subnormals in the divisor. */
    vboolN_2 blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(b), FLT_MIN);
    BBE_MULN_2XF32T(b,b,8388608.f,blt);
    t = BBE_RECIP0N_2XF32(b);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(r,t,b); BBE_MULAN_2XF32(t,r,t);
    /* r <- ~a/b */
    r = BBE_MULN_2XF32(a,t);
    /* Modified Newton-Raphson iteration for a quotient */
    BBE_MULSN_2XF32(a,b,r); BBE_MULAN_2XF32(r,t,a);
    BBE_MULN_2XF32T(r,r,8388608.f,blt);
  } else {
    /* This variant overflows when the divisor is subnormal. */
    t = BBE_RECIP0N_2XF32(b);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(r,t,b); BBE_MULAN_2XF32(t,r,t);
    /* r <- ~a/b */
    r = BBE_MULN_2XF32(a,t);
    /* Modified Newton-Raphson iteration for a quotient */
    BBE_MULSN_2XF32(a,b,r); BBE_MULAN_2XF32(r,t,a);
  }
  return (r);

} /* __IT_FDIVN_2XF32() */

/* Predicated fast real vector division based on reciprocal 
 * approximation. If argument b may be subnormal, set scaleDivsr
 * to non-zero to prevent overflow (accuracy still degrades if 
 * divisor is  subnormal and dividend is close to subnormal).
 * Accuracy is not worse than 1 ULP for normal arguments. */
#define IT_FDIVN_2XF32T(a,b,c,p,scaleDivsr)     __IT_FDIVN_2XF32T(&(a),&(b),&(c),&(p),scaleDivsr)
ATTRIBUTE_ALWAYS_INLINE
inline_ void __IT_FDIVN_2XF32T( void * pa, const void * pb, const void * pc, const void * pp, int scaleDivsr )
{
  xb_vecN_2xf32 a,b,c,r,t;
  vboolN_2 bp;
  a = *(xb_vecN_2xf32*)pa;
  b = *(xb_vecN_2xf32*)pb;
  c = *(xb_vecN_2xf32*)pc;
  bp = *(vboolN_2*)pp;
  /* t <- ~1/b */
  if (scaleDivsr) {
    /* This variant correctly handles subnormals in the divisor. */
    vboolN_2 blt = BBE_OLTN_2XF32T(BBE_ABSN_2XF32(c), FLT_MIN, bp);
    BBE_MULN_2XF32T(c,c,8388608.f,blt);
    t = BBE_CONSTN_2XF32(0); BBE_RECIP0N_2XF32T(t,c,bp);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32T(r,t,c,bp); BBE_MULAN_2XF32T(t,r,t,bp);
    /* r <- ~a/b */
    BBE_MULN_2XF32T(a,b,t,bp);
    /* Modified Newton-Raphson iteration for a quotient */
    BBE_MULSN_2XF32T(b,c,a,bp); BBE_MULAN_2XF32T(a,t,b,bp);
    BBE_MULN_2XF32T(a,a,8388608.f,blt);
  } else {
    /* This variant overflows when the divisor is subnormal. */
    t = BBE_CONSTN_2XF32(0); BBE_RECIP0N_2XF32T(t,c,bp);
    /* Newton-Raphson refinement iteration for a reciprocal */
    r = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32T(r,t,c,bp); BBE_MULAN_2XF32T(t,r,t,bp);
    /* r <- ~a/b */
    BBE_MULN_2XF32T(a,b,t,bp);
    /* Modified Newton-Raphson iteration for a quotient */
    BBE_MULSN_2XF32T(b,c,a,bp); BBE_MULAN_2XF32T(a,t,b,bp);
  }
  *(xb_vecN_2xf32*)pa = a;

} /* __IT_FDIVN_2XF32T() */

/* Absolute value of a scalar complex number */
#define IT_ABSCF32(x)  __IT_ABSCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xtfloat __IT_ABSCF32(const void * px)
{
  /*
   * Based on cabs code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   *
   * float32_t _cabsf( complex_float x )
   * {
   *   float32_t t,xr,xi;
   *   xr = fabsf(crealf(x));
   *   xi = fabsf(cimagf(x));
   *   if (xi>0) {
   *     if (xi>xr) { t = xr; xr = xi; xi = t; }
   *     t = ( xi==xr ? 1 : xi/xr );  // Avoid Inf/Inf
   *     return ( xr*sqrtf(1+t*t) );
   *   } else {
   *     return ( xi==xi ? xr : xi );
   *   }
   * }
   */

  xtfloat xr,xi,a,b,r,t;
  vbool1 ble,bgt,beq,bun;
  xr = XT_ABS_S(BBE_CREALCF32(*(xtcomplexfloat*)px));
  xi = XT_ABS_S(BBE_CIMAGCF32(*(xtcomplexfloat*)px));
  ble = XT_ULE_S(xi, XT_CONST_S(0));
  bgt = XT_OLT_S(xr,xi);
  beq = XT_OEQ_S(xr,xi);
  bun = XT_UN_S(xi,xi);
  a = xi; XT_MOVT_S(a,xr,bgt);
  b = xr; XT_MOVT_S(b,xi,bgt);
  r = XT_RECIP0_S(b);
  t = XT_CONST_S(1);
  XT_MSUB_S(t,r,b);
  XT_MADD_S(r,t,r);
  t = XT_MUL_S(a,r);
  XT_MSUB_S(a,b,t);
  XT_MADD_S(t,a,r);
  BBE_CONSTF32T(t,1,beq);
  r = XT_CONST_S(1);
  XT_MADD_S(r,t,t);
  r = XT_MUL_S(b, XT_SQRT_S(r));
  t = xr; XT_MOVT_S(t,xi,bun);
  XT_MOVT_S(r,t,ble);
  return (r);

} /* __IT_ABSCF32() */

/* Absolute value of each complex number in a vector. Results are
 * replicated twice in the output vector, in positions that correspond
 * to real/imaginary parts of the respective complex element in the
 * input vector. */
#define IT_ABSN_4XCF32(x)  __IT_ABSN_4XCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_2xf32 __IT_ABSN_4XCF32(const void * px)
{
  /*
   * Based on cabs code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   *
   * float32_t _cabsf( complex_float x )
   * {
   *   float32_t t,xr,xi;
   *   xr = fabsf(crealf(x));
   *   xi = fabsf(cimagf(x));
   *   if (xi>0) {
   *     if (xi>xr) { t = xr; xr = xi; xi = t; }
   *     t = ( xi==xr ? 1 : xi/xr );  // Avoid Inf/Inf
   *     return ( xr*sqrtf(1+t*t) );
   *   } else {
   *     return ( xi==xi ? xr : xi );
   *   }
   * }
   */

  xb_vecN_2xf32 x,xr,xi,a,b,r,t;
  vboolN_2 ble,bgt,beq,bun;
  x = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px);
  xr = BBE_ABSN_2XF32(BBE_SHFLN_2XF32I(x, BBE_SHFLI_DUPLICATE_2_EVEN));
  xi = BBE_ABSN_2XF32(BBE_SHFLN_2XF32I(x, BBE_SHFLI_DUPLICATE_2_ODD));
  ble = BBE_ULEN_2XF32(xi, BBE_CONSTN_2XF32(0));
  bgt = BBE_OLTN_2XF32(xr,xi);
  beq = BBE_OEQN_2XF32(xr,xi);
  bun = BBE_UNN_2XF32(xi,xi);
  a = BBE_MOVN_2XF32T(xr,xi,bgt);
  b = BBE_MOVN_2XF32T(xi,xr,bgt);
  r = BBE_RECIP0N_2XF32(b);
  t = BBE_CONSTN_2XF32(1);
  BBE_MULSN_2XF32(t,r,b);
  BBE_MULAN_2XF32(r,t,r);
  t = BBE_MULN_2XF32(a,r);
  BBE_MULSN_2XF32(a,b,t);
  BBE_MULAN_2XF32(t,a,r);
  BBE_CONSTN_2XF32T(t,1,beq);
  r = BBE_CONSTN_2XF32(1);
  BBE_MULAN_2XF32(r,t,t);
  r = BBE_MULN_2XF32(b, BBE_SQRTN_2XF32(r));
  t = BBE_MOVN_2XF32T(xi,xr,bun);
  r = BBE_MOVN_2XF32T(t,r,ble);
  return (r);

} /* __IT_ABSN_4XCF32() */

/* Vector squared magnitude of a complex scalar. */
#define IT_ABS2N_4XCF32(x)  __IT_ABS2N_4XCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_2xf32 __IT_ABS2N_4XCF32(const void * px)
{
    xb_vecN_2xf32 x;
    x = BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px);
    x = BBE_MULN_2XF32(x, x);
    x = BBE_ADDN_2XF32(x, BBE_SHFLN_2XF32I(x, BBE_SHFLI_SWAP_2));
    return x;
} /* __IT_ABS2N_4XCF32() */

/* Squared magnitude of a complex scalar. */
#define IT_ABS2CF32(x)  __IT_ABS2CF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xtfloat __IT_ABS2CF32(const void * px)
{
  xtcomplexfloat x,y;
  x = *(xtcomplexfloat*)px;
  y = BBE_MULMCF32(x,x,0,0);
  BBE_MULMASCF32(y,x,x,0,3);
  return (BBE_CREALCF32(y));

} /* __IT_ABS2CF32() */

/* Vector sum of real/imaginary absolute values for a complex scalar */
#define IT_CMAGN_4XCF32(x)  __IT_CMAGN_4XCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_2xf32 __IT_CMAGN_4XCF32(const void * px){
    xb_vecN_2xf32 x;
    x = BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(*(xb_vecN_4xcf32*)px));
    x = BBE_ADDN_2XF32(x, BBE_SHFLN_2XF32I(x, BBE_SHFLI_SWAP_2));
    return x;
} /* __IT_CMAGN_4XCF32() */

/* Sum of real/imaginary absolute values for a complex scalar */
#define IT_CMAGCF32(a)    XT_ADD_S(XT_ABS_S(BBE_CREALCF32(a)),XT_ABS_S(BBE_CIMAGCF32(a)))

/* Vector complex floating-point square root, single precision. */
#define IT_SQRTN_4XCF32(x)    __IT_SQRTN_4XCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xb_vecN_4xcf32 __IT_SQRTN_4XCF32(const void * px)
{
    /*
    * Based on csqrt code from "Similarity Reduction of a General Matrix to
    * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
    * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
    *
    * complex_float _csqrtf( complex_float x )
    * {
    *   float32_t xr,xi,h;
    *   xr = crealf(x); xi = cimagf(x);
    *   h = sqrtf((fabsf(xr) + _cabsf(x))/2);
    *   if (xi!=0) xi /= 2*h;
    *   if (xr<0) {
    *     xr = ( xi<0 ? -xi : xi );
    *     xi = ( xi<0 ? -h : h );
    *   } else {
    *     xr = h;
    *   }
    *   return ( _makecomplexf(xr,xi) );
    * }
    */

    xb_vecN_4xcf32 x;
    xb_vecN_2xf32 xr, xi, f, h, t, y, z;
    vboolN_2 brn, bin, bnz;
    x = (*(xb_vecN_4xcf32*)px);
    xr = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMN_4XCF32(x), BBE_SHFLI_DUPLICATE_2_EVEN);
    xi = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMN_4XCF32(x), BBE_SHFLI_DUPLICATE_2_ODD);
    brn = BBE_OLTN_2XF32(xr, BBE_CONSTN_2XF32(0));
    bin = BBE_OLTN_2XF32(xi, BBE_CONSTN_2XF32(0));
    bnz = BBE_UNEQN_2XF32(xi, BBE_CONSTN_2XF32(0));
    t = BBE_ADDN_2XF32(BBE_ABSN_2XF32(xr), IT_ABSN_4XCF32(x));
    t = BBE_MULN_2XF32(t, BBE_CONSTN_2XF32(3));
    h = BBE_SQRTN_2XF32(t);
    t = BBE_MULN_2XF32(h, BBE_CONSTN_2XF32(2));
    f = BBE_CONSTN_2XF32(0); IT_FDIVN_2XF32T(f, xi, t, bnz, 0);
    y = BBE_SELN_2XF32I(f, h, BBE_SELI_INTERLEAVE_2_EVEN);
    z = BBE_SELN_2XF32I(h, f, BBE_SELI_INTERLEAVE_2_EVEN);
    BBE_NEGN_2XF32T(z, z, bin);
    y = BBE_MOVN_2XF32T(z, y, brn);
    return BBE_MOVN_4XCF32_FROMN_2XF32(y);

} /* __IT_SQRTN_4XCF32() */

/* Scalar complex floating-point square root, single precision. Returns the
 * square root that lies in the right half of the complex plane. */
#define IT_SQRTCF32(x)    __IT_SQRTCF32(&(x))
ATTRIBUTE_ALWAYS_INLINE
inline_ xtcomplexfloat __IT_SQRTCF32( const void * px )
{
  /*
   * Based on csqrt code from "Similarity Reduction of a General Matrix to
   * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
   * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
   *
   * complex_float _csqrtf( complex_float x ) 
   * {
   *   float32_t xr,xi,h;
   *   xr = crealf(x); xi = cimagf(x);
   *   h = sqrtf((fabsf(xr) + _cabsf(x))/2);
   *   if (xi!=0) xi /= 2*h;
   *   if (xr<0) {
   *     xr = ( xi<0 ? -xi : xi );
   *     xi = ( xi<0 ? -h : h );
   *   } else {
   *     xr = h;
   *   }
   *   return ( _makecomplexf(xr,xi) );
   * }
   */

  xb_vecN_4xcf32 x;
  xb_vecN_2xf32 xr,xi,f,h,t,y,z;
  vboolN_2 brn,bin,bnz;
  x = BBE_MOVN_4XCF32_FROMCF32(*(xtcomplexfloat*)px);
  xr = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMN_4XCF32(x), BBE_SHFLI_DUPLICATE_2_EVEN);
  xi = BBE_SHFLN_2XF32I(BBE_MOVN_2XF32_FROMN_4XCF32(x), BBE_SHFLI_DUPLICATE_2_ODD);
  brn = BBE_OLTN_2XF32(xr, BBE_CONSTN_2XF32(0));
  bin = BBE_OLTN_2XF32(xi, BBE_CONSTN_2XF32(0));
  bnz = BBE_UNEQN_2XF32(xi, BBE_CONSTN_2XF32(0));
  t = BBE_ADDN_2XF32(BBE_ABSN_2XF32(xr), IT_ABSN_4XCF32(x));
  t = BBE_MULN_2XF32(t, BBE_CONSTN_2XF32(3));
  h = BBE_SQRTN_2XF32(t);
  t = BBE_MULN_2XF32(h, BBE_CONSTN_2XF32(2));
  f = BBE_CONSTN_2XF32(0); IT_FDIVN_2XF32T(f,xi,t,bnz,0);
  y = BBE_SELN_2XF32I(f, h, BBE_SELI_INTERLEAVE_2_EVEN);
  z = BBE_SELN_2XF32I(h, f, BBE_SELI_INTERLEAVE_2_EVEN);
  BBE_NEGN_2XF32T(z,z,bin);
  y = BBE_MOVN_2XF32T(z,y,brn);
  return BBE_MOVCF32_FROMN_2XF32(y);

} /* __IT_SQRTCF32() */

#endif /* HAVE_VFPU */

#endif /* __VFPU_MATH_H */

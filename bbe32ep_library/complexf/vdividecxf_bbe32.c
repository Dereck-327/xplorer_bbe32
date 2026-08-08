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
  NatureDSP_Baseband library. Complex Math functions
    Division of Complex Numbers
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, vdividecxf, ( complex_float   * restrict z,
                          const complex_float   * restrict x,
                          const complex_float   * restrict y,
                          int N ))
#else

#ifndef BBE_MOVN_2XF32T
#define BBE_MOVN_2XF32T(a, b, c) \
        BBE_MOVN_2XF32_FROMNX16( BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(a),BBE_MOVNX16_FROMN_2XF32(b),BBE_MOVN_FROMN_2(c)) )
#endif

/*-------------------------------------------------------------------------
Division of Complex Numbers

Representation:
vdividecx         16-bit signed fixed-point format (Q15)
vdividecxf,cdivf  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of real
and imaginary components of the division result. Fixed-point format for the
fractional part is 16-bit Q(15-exp), where exp denotes the respective 
exponential value. Full division result can be restored in 48-bit Q31 format
by sign extending the fractional part to 64 bits and shifting it to the left
by 16+exp bit positions

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number          (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number (see *)
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                               (fixed-point functions)
          |         | defined |

* - quotent of infinite number to finite number may have NaN in its real or 
    imaginary part depending on signs of input arguments

Input domain for 'fast' version vfastdividecxf is limited by usage of direct 
formula without normalization of inputs, so input domain is:
|real(x)*real(y)|<Inf
|imag(x)*imag(y)|<Inf
|real(y)*real(y)|<Inf
|imag(y)*imag(y)|<Inf
1.1755e-038 < |real(x)*real(y) + imag(x)*imag(y)| < Inf
1.1755e-038 < |real(y)*real(y) + imag(y)*imag(y)| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
vdividecx         1.5 LSB of the fractional part
vdividecxf,cdivf  2 ULP
vfastdividecxf    2 ULP

Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdividecx:
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[2*N]  Exponential part of quotients, -14..16. Exponential values for
          real and imaginary components are interleaved, with exponential
          parts of real components going first (at even indices)
vdividecxf,cdivf
z[N]      Quotients
Restrictions:
x,y,z,fract,exp   Aligned on 32-byte boundary
x,y,z,fract,exp   Must not overlap
N                 Multiple of 8 (vdividecx) or 4 (vdividecxf, vfastdividecxf)
---------------------------------------------------------------------------*/

void vdividecxf ( complex_float   * restrict z,
            const complex_float   * restrict x,
            const complex_float   * restrict y,
            int N )
{
#if 1
  static const uint32_t cnst_tbl[] =
  {
      /* Constants used in computation */
      (126<<16)|(3*126+2*127),/* exponents used for      */
                              /* normalization of values */
      0x7F800000,  /* +INF */
      0,
      (127+23)<<23,/* 2^23 */
      0x007FFFFF   /* max subnormal number */
  };
  const xb_vecNx16    *restrict TBL;
  const xb_vecN_2xf32 *restrict X;
  const xb_vecN_2xf32 *restrict Y;
  const xb_vecN_2xf32 *restrict Zrd;
        xb_vecN_2xf32 *restrict Zwr;
  xb_vecN_2xf32 a,  b,  c,  d,
                aa, ab, ac, ad;
  xb_vecN_2xf32 cc, dd, den, v, w;
  xb_vecN_2xf32 xmax, ymax, t1, t2;
  xb_vecNx16 ex, ey, sy, sy2, zeroi, tmp, e1, e2;
  xb_vecN_2xf32 plusInf, maxDenorm, pow2_23, zerof;
  vboolN_2 b_yinf, b_dendenorm, b_denzero;
  valign X_va, Y_va, Zwr_va, Zrd_va;
  int n, nZ;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 4) == 0);

  TBL = (const xb_vecNx16 *)cnst_tbl;

  /*
   * First stage: make normalization of denominator
   * and save intermediate results to the output.
   */

  X   = (const xb_vecN_2xf32 *)x;
  Y   = (const xb_vecN_2xf32 *)y;
  Zwr = (      xb_vecN_2xf32 *)z;
  X_va = BBE_LAN_2XF32_PP(X);
  Y_va = BBE_LAN_2XF32_PP(Y);
  Zwr_va = BBE_ZALIGN();
  nZ = N*sizeof(float32_t);

  BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  e1 = BBE_REPNX16(tmp, 1);  e2 = BBE_REPNX16(tmp, 0);
  BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  plusInf = BBE_MOVN_2XF32_FROMNX16(tmp);
  BBE_LPNX16_IP(tmp, TBL, -2*(int)sizeof(int16_t));  zeroi = BBE_REPNX16(tmp, 0);
    
  for ( n=0; n<(N+BBE_SIMD_WIDTH/2-1)/(BBE_SIMD_WIDTH/2); n++ )
  {
    /* a=x.re, b=x.im */
    BBE_LAVN_2XF32_XP( t1, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_LAVN_2XF32_XP( t2, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_DSELN_2XF32I(b, a, t2, t1, BBE_DSELI_DEINTERLEAVE_2);
    /* c=y.re, d=y.im */
    BBE_LAVN_2XF32_XP( t1, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_LAVN_2XF32_XP( t2, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_DSELN_2XF32I(d, c, t2, t1, BBE_DSELI_DEINTERLEAVE_2);
    /* Find maximum absolute values of x and y */
    aa = BBE_ABSN_2XF32(a);    ab = BBE_ABSN_2XF32(b);
    ac = BBE_ABSN_2XF32(c);    ad = BBE_ABSN_2XF32(d);
    xmax = BBE_MAXN_2XF32(aa, ab);
    ymax = BBE_MAXN_2XF32(ac, ad);

    /* Get the exponent of xmax and ymax
     * and calculate sy=3*126-ey-(MAX(126, MAX(ey, ex))) */
    ex = BBE_MOVNX16_FROMN_2XF32(xmax);
    ey = BBE_MOVNX16_FROMN_2XF32(ymax);
    ex = BBE_SRLINX16(ex, 7);
    ey = BBE_SRLINX16(ey, 7);
    sy2 = BBE_MAXNX16(ex, ey);
    sy2 = BBE_MAXNX16(e1, sy2);
    sy = BBE_SUBNX16(e2, ey);
    sy = BBE_SUBNX16(sy, sy2);/* -132...252 */

    /* Implementaion of 
     * cc = ldexpf(c, sy)            
     * dd = ldexpf(d, sy) */
    sy2 = BBE_SRLINX16(sy,  1);  sy = BBE_SUBNX16(sy, sy2);
    sy2 = BBE_SLLINX16(sy2, 7);
    sy = BBE_SLLINX16(sy, 7);
    sy  = BBE_SELNX16I(sy , zeroi, BBE_SELI_INTERLEAVE_1_ODD);
    sy2 = BBE_SELNX16I(sy2, zeroi, BBE_SELI_INTERLEAVE_1_ODD);

    t1 = BBE_MOVN_2XF32_FROMNX16(sy2);
    t2 = BBE_MOVN_2XF32_FROMNX16(sy);
    cc = c*t1;  cc = cc*t2;
    dd = d*t1;  dd = dd*t2;
    /*--------------------*/

    /* Save intermediate results to the output */
    BBE_SAVN_2XF32_XP( cc, Zwr_va, Zwr, nZ );
    BBE_SAVN_2XF32_XP( dd, Zwr_va, Zwr, nZ );
    nZ -= 2*BBE_SIMD_WIDTH;
  }
  BBE_SAN_2XF32POS_FP( Zwr_va, Zwr );


  /*
   * Second stage: perform division with the use of normalized denominator
   */

  X   = (const xb_vecN_2xf32 *)x;
  Y   = (const xb_vecN_2xf32 *)y;
  Zrd = (const xb_vecN_2xf32 *)z;
  Zwr = (      xb_vecN_2xf32 *)z;
  X_va   = BBE_LAN_2XF32_PP(X);
  Y_va   = BBE_LAN_2XF32_PP(Y);
  Zrd_va = BBE_LAN_2XF32_PP(Zrd);
  Zwr_va = BBE_ZALIGN();
  nZ = N*sizeof(float32_t);

  for ( n=0; n<(N+BBE_SIMD_WIDTH/2-1)/(BBE_SIMD_WIDTH/2); n++ )
  {
    /* Reload constants used in computations */
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  plusInf = BBE_MOVN_2XF32_FROMNX16(tmp);
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  zerof = BBE_MOVN_2XF32_FROMNX16(tmp);
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  pow2_23 = BBE_MOVN_2XF32_FROMNX16(tmp);
    BBE_LPNX16_XP(tmp, TBL, -6*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  maxDenorm = BBE_MOVN_2XF32_FROMNX16(tmp);

    /* a=x.re, b=x.im */
    BBE_LAVN_2XF32_XP( t1, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_LAVN_2XF32_XP( t2, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_DSELN_2XF32I(b, a, t2, t1, BBE_DSELI_DEINTERLEAVE_2);
    /* c=y.re, d=y.im */
    BBE_LAVN_2XF32_XP( t1, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_LAVN_2XF32_XP( t2, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_DSELN_2XF32I(d, c, t2, t1, BBE_DSELI_DEINTERLEAVE_2);

    /* Find maximum absolute values of x */
    aa = BBE_ABSN_2XF32(a);    ab = BBE_ABSN_2XF32(b);
    xmax=aa;
    BBE_MAXN_2XF32T(xmax,xmax,ab,BBE_OEQN_2XF32(xmax,xmax));
    
    /* Load previously computed normalized denominator */
    BBE_LAVN_2XF32_XP( cc, Zrd_va, Zrd, nZ );
    BBE_LAVN_2XF32_XP( dd, Zrd_va, Zrd, nZ );
    nZ -= 2*BBE_SIMD_WIDTH;

    /* v=a*cc+b*dd; */
    t1 = BBE_MULN_2XF32(a,  cc);
    t2 = BBE_MULN_2XF32(b,  dd);
    v  = BBE_ADDN_2XF32(t1, t2);
    /* w=b*cc-a*dd; */
    t1 = BBE_MULN_2XF32(b,  cc);
    t2 = BBE_MULN_2XF32(a,  dd);
    w  = BBE_SUBN_2XF32(t1, t2);
    /* den=c*cc+d*dd; */
    t1  = BBE_MULN_2XF32(c, cc);
    t2  = BBE_MULN_2XF32(d, dd);
    den = BBE_ADDN_2XF32(t1, t2);

    /* prepare numerator and denominator for division */
    b_dendenorm = BBE_OLEN_2XF32(den, maxDenorm);
    t1 = BBE_MULN_2XF32(v  , pow2_23); v   = BBE_MOVN_2XF32T(t1, v  , b_dendenorm);
    t1 = BBE_MULN_2XF32(w  , pow2_23); w   = BBE_MOVN_2XF32T(t1, w  , b_dendenorm);
    t1 = BBE_MULN_2XF32(den, pow2_23); den = BBE_MOVN_2XF32T(t1, den, b_dendenorm);
    b_yinf = BBE_ANDNOTBN_2(BBE_OEQN_2XF32(den,plusInf), BBE_UEQN_2XF32(xmax,plusInf));

    /* divide v, w by den */
    b_denzero = BBE_OEQN_2XF32(den, zerof);
    t1 = BBE_RECIPN_2XF32(den);
    t1 = BBE_MOVN_2XF32T(plusInf, t1, b_denzero);
    v  = BBE_MULN_2XF32(t1, v);
    w  = BBE_MULN_2XF32(t1, w);
    BBE_CONSTN_2XF32T(v, 0, b_yinf);
    BBE_CONSTN_2XF32T(w, 0, b_yinf);

    /* z.re=v; z.im=w */
    BBE_DSELN_2XF32I(t2, t1, w, v, BBE_DSELI_INTERLEAVE_2);
    BBE_SAVN_2XF32_XP( t1, Zwr_va, Zwr, (uint8_t*)z + N*sizeof(complex_float) - (uint8_t*)Zwr );
    BBE_SAVN_2XF32_XP( t2, Zwr_va, Zwr, (uint8_t*)z + N*sizeof(complex_float) - (uint8_t*)Zwr );
  }
  BBE_SAN_2XF32POS_FP( Zwr_va, Zwr );

#else
  static const union ufloat32uint32 cnst_tbl[] =
  {
      /* Constants used in computation */
      {(126<<16)|(3*126+2*127)},/*   exponents         */
                                /*   for normalization */
      {0x7F800000},  /* +INF */
      {0},
      {(127+23)<<23}, /* 2^23 */
      {0x007FFFFF}  /* max subnormal number */
  };
  const xb_vecN_2xf32 *restrict TBL;
  const xb_vecN_2xf32 *restrict X;
  const xb_vecN_2xf32 *restrict Y;
        xb_vecN_2xf32 *restrict Z;
  xb_vecN_2xf32 a,  b,  c,  d,
                aa, ab, ac, ad;
  xb_vecN_2xf32 cc, dd, den, v, w;
  xb_vecN_2xf32 xmax, ymax, t1, t2;
  xb_vecNx16 ex, ey, sy, sy2, zeroi, tmp, e1, e2;
  xb_vecN_2xf32 plusInf, maxDenorm, pow2_23, zerof;
  vboolN_2 b_yinf, b_dendenorm, b_denzero;
  valign X_va, Y_va, Z_va;
  int n;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 4) == 0);

  TBL = (const xb_vecN_2xf32 *)cnst_tbl;
  X   = (const xb_vecN_2xf32 *)x;
  Y   = (const xb_vecN_2xf32 *)y;
  Z   = (      xb_vecN_2xf32 *)z;
  X_va = BBE_LAN_2XF32_PP(X);
  Y_va = BBE_LAN_2XF32_PP(Y);
  Z_va = BBE_ZALIGN();

  for ( n=0; n<(N+BBE_SIMD_WIDTH/2-1)/(BBE_SIMD_WIDTH/2); n++ )
  {
    /* a=x.re, b=x.im */
    BBE_LAVN_2XF32_XP( t1, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_LAVN_2XF32_XP( t2, X_va, X, (uint8_t*)x + N*sizeof(complex_float) - (uint8_t*)X );
    BBE_DSELN_2XF32I(b, a, t2, t1, BBE_DSELI_DEINTERLEAVE_2);
    /* c=y.re, d=y.im */
    BBE_LAVN_2XF32_XP( t1, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_LAVN_2XF32_XP( t2, Y_va, Y, (uint8_t*)y + N*sizeof(complex_float) - (uint8_t*)Y );
    BBE_DSELN_2XF32I(d, c, t2, t1, BBE_DSELI_DEINTERLEAVE_2);
    /* Find maximum absolute values of x and y */
    aa = BBE_ABSN_2XF32(a);    ab = BBE_ABSN_2XF32(b);
    ac = BBE_ABSN_2XF32(c);    ad = BBE_ABSN_2XF32(d);
    xmax=aa;
    ymax=ac;
    BBE_MAXN_2XF32T(xmax,xmax,ab,BBE_OEQN_2XF32(xmax,xmax));
    BBE_MAXN_2XF32T(ymax,ymax,ad,BBE_OEQN_2XF32(ymax,ymax));

    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  e1 = BBE_REPNX16(tmp, 1);  e2 = BBE_REPNX16(tmp, 0);
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  plusInf = BBE_MOVN_2XF32_FROMNX16(tmp);
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  zeroi = BBE_REPNX16(tmp, 0);  zerof = BBE_MOVN_2XF32_FROMNX16(zeroi);
    BBE_LPNX16_IP(tmp, TBL,  2*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  pow2_23 = BBE_MOVN_2XF32_FROMNX16(tmp);
    BBE_LPNX16_XP(tmp, TBL, -8*(int)sizeof(int16_t));  tmp = BBE_REPNX16C(tmp, 0);  maxDenorm = BBE_MOVN_2XF32_FROMNX16(tmp);
    /* Get the exponent of xmax and ymax
     * and calculate sy=3*126-ey-(MAX(126, MAX(ey, ex))) */
    ex = BBE_MOVNX16_FROMN_2XF32(xmax);
    ey = BBE_MOVNX16_FROMN_2XF32(ymax);
    ex = BBE_SRLINX16(ex, 7);
    ey = BBE_SRLINX16(ey, 7);
    sy2 = BBE_MAXNX16(ex,   ey);
    sy2 = BBE_MAXNX16(e1, sy2);
    sy = BBE_SUBNX16(e2, ey);
    sy = BBE_SUBNX16(sy, sy2);/* -132...252 */

    /* Implementaion of 
     * cc = ldexpf(c, sy)            
     * dd = ldexpf(d, sy) */
    sy2 = BBE_SRLINX16(sy,  1);  sy = BBE_SUBNX16(sy, sy2);
    sy = BBE_SELNX16I(sy2, sy, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
    sy = BBE_SLLINX16(sy, 7);
    BBE_DSELNX16I(sy2, sy, sy, zeroi, BBE_DSELI_INTERLEAVE_1);

    t1 = BBE_MOVN_2XF32_FROMNX16(sy2);
    t2 = BBE_MOVN_2XF32_FROMNX16(sy);

    cc = c*t1;  cc = cc*t2;
    dd = d*t1;  dd = dd*t2;
    /*--------------------*/
    /* v=a*cc+b*dd; */
    t1 = BBE_MULN_2XF32(a,  cc);
    t2 = BBE_MULN_2XF32(b,  dd);
    v  = BBE_ADDN_2XF32(t1, t2);
    /* w=b*cc-a*dd; */
    t1 = BBE_MULN_2XF32(b,  cc);
    t2 = BBE_MULN_2XF32(a,  dd);
    w  = BBE_SUBN_2XF32(t1, t2);
    /* den=c*cc+d*dd; */
    t1  = BBE_MULN_2XF32(c, cc);
    t2  = BBE_MULN_2XF32(d, dd);
    den = BBE_ADDN_2XF32(t1, t2);
    /* prepare nominator and denominator for division */
    b_yinf = BBE_ANDNOTBN_2(BBE_OEQN_2XF32(ymax,plusInf), BBE_UEQN_2XF32(xmax,plusInf));
    b_dendenorm = BBE_OLEN_2XF32(den, maxDenorm);
    t1 = BBE_MULN_2XF32(v  , pow2_23); v   = BBE_MOVN_2XF32T(t1, v  , b_dendenorm);
    t1 = BBE_MULN_2XF32(w  , pow2_23); w   = BBE_MOVN_2XF32T(t1, w  , b_dendenorm);
    t1 = BBE_MULN_2XF32(den, pow2_23); den = BBE_MOVN_2XF32T(t1, den, b_dendenorm);
    BBE_CONSTN_2XF32T(v  , 0, b_yinf);
    BBE_CONSTN_2XF32T(w  , 0, b_yinf);
    BBE_CONSTN_2XF32T(den, 1, b_yinf);

    /* divide v, w by den */
    b_denzero = BBE_OEQN_2XF32(den, zerof);
    t1 = BBE_RECIPN_2XF32(den);
    t1 = BBE_MOVN_2XF32T(plusInf, t1, b_denzero);
    v  = BBE_MULN_2XF32(t1, v);
    w  = BBE_MULN_2XF32(t1, w);

    /* z.re=v; z.im=w */
    BBE_DSELN_2XF32I(t2, t1, w, v, BBE_DSELI_INTERLEAVE_2);
    BBE_SAVN_2XF32_XP( t1, Z_va, Z, (uint8_t*)z + N*sizeof(complex_float) - (uint8_t*)Z );
    BBE_SAVN_2XF32_XP( t2, Z_va, Z, (uint8_t*)z + N*sizeof(complex_float) - (uint8_t*)Z );
  }
  BBE_SAN_2XF32POS_FP( Z_va, Z );
#endif
} /* vdividecxf() */

#endif/* !HAVE_VFPU */

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
  NatureDSP_Baseband library. Math functions
    Raise To a Power
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Tables */
#include "pow2f_tbl.h"

#define MAX_ALLOCA_SF32   ( (MAX_ALLOCA_SZ/2)/sizeof(float32_t) )
/*-------------------------------------------------------------------------
Raise To a Power

Description: These functions compute the value of the first argument x 
raised to the power of the second argument y.

Data format: IEEE-754 Std. single precision floating-point.

Special cases:
      x   |   y    | Result |  Extra Conditions    
  --------+--------+--------+---------------------
    +/-0  | y      | +/-inf | odd y<0
    +/-0  | y      | +inf   | even y<0
    +/-0  | y      | +/-0   | odd y>0
    +/-0  | y      | 0      | even y>0
    +/-1  | +/-inf | 1      | 
    1     | y      | 1      | any y including NaN
    x     | +/-0   | 1      | any x including NaN
    x     | y      | NaN    | finite x<0 and finite non-integer y (see note 2)
    x     | -inf   | +inf   | |x|<1
    x     | -inf   | 0      | |x|>1
    x     | +inf   | 0      | |x|<1
    x     | +inf   | +inf   | |x|>1
    -inf  | y      | -0     | y an odd integer <0
    -inf  | y      | 0      | y<0 and not an odd integer
    -inf  | y      | -inf   | y an odd integer >0
    -inf  | y      | +inf   | y>0 and not an odd integer
    +inf  | y      | 0      | y<0
    +inf  | y      | +inf   | y>0

Accuracy: 2 ULP under condition that |y|<=100

Notes:
1. Raise to a power functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-point
   exceptions.
2. If x<0 is finite, y is finite and not an integer value, then the respective
   result z is set to NaN, errno is assigned the value EDOM, and the "invalid"
   floating-point exception is raised.

Parameters:
Input:
x[N]    Input data
y[N]    Power values
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vpowf,( float32_t * restrict z, 
       const float32_t * restrict x, 
       const float32_t * restrict y, 
       int N ))
#else
void vpowf ( float32_t * restrict z, 
       const float32_t * restrict x, 
       const float32_t * restrict y, 
       int N )
{
  int i, n, _N;
  __fenv_t fenv;
  xb_vecN_2xf32 xin, y0, one, ef_reg;
  float32_t ALIGN(32) scratch0[MAX_ALLOCA_SF32];
  float32_t ALIGN(32) scratch1[MAX_ALLOCA_SF32];
  const xb_vecN_2xf32  * restrict pX;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
        xb_vecN_2xf32  * restrict pZ = (      xb_vecN_2xf32  *)z;
  const xb_vecN_2xf32  * restrict SCR0_rd;
        xb_vecN_2xf32  * restrict SCR0_wr;
  const xb_vecN_2xf32  * restrict SCR1_rd;
        xb_vecN_2xf32  * restrict SCR1_wr;
  const xb_vecN_2xf32  * restrict TBL_LOG2;
  const xb_vecN_2xf32  * restrict TBL_POW2;
  const int32_t        * restrict TBL;
  xb_vecN_2xf32 v_edom, v_erange, v_fe_inv, v_fe_divz, v_fe_ovfl;
  int er_edom, er_erange, fe_inv, fe_divz, fe_ovfl;
  vboolN_2 b_spec_nan, b_xyinf, b_zinf, b_xyfin_zinf, b_snan;
  vboolN_2 b_edom, b_erange, b_fe_inv, b_fe_divz, b_fe_ovfl;
  xb_vecN_2x32Uv flags_x, flags_y, maskx, v_snan, flags_xy;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) return;
  /* Table of different constants used in computations */
  static const int32_t c_tbl[] =
  {
    -126,
    -150,
    (int32_t)0x007FFFFF,/* max denormalized floating-point number / mantissa mask */
    (int32_t)0x4B800000,/* 2^24 */
    (int32_t)0x3F3504F3,/* sqrt(0.5) */
    (int32_t)0x3F000000,/*  0.5 */
    (int32_t)0xBF000000,/* -0.5 */
    -252,
    254
  };
  TBL_LOG2 = (const xb_vecN_2xf32  *)log2f_coef;
  TBL_POW2 = (const xb_vecN_2xf32  *)pow2f_coef;
  
  er_edom = er_erange = fe_inv = fe_divz = fe_ovfl = 0;
  b_edom = BBE_XORBN_2(b_edom, b_edom);
  b_erange = BBE_XORBN_2(b_erange, b_erange);
  b_fe_inv = BBE_XORBN_2(b_fe_inv, b_fe_inv);
  b_fe_divz = BBE_XORBN_2(b_fe_divz, b_fe_divz);
  b_fe_ovfl = BBE_XORBN_2(b_fe_ovfl, b_fe_ovfl);
  /* Clear exception enable flags and exception status flags. */
  __feholdexcept(&fenv);
  /*------------------------------*/
  /* Compute x^y as 2^(y*log2(x)) */
  /*------------------------------*/
  for (i = 0; i < (int)((N + MAX_ALLOCA_SF32 - 1) / MAX_ALLOCA_SF32); i++)
  {
    _N = XT_MIN(MAX_ALLOCA_SF32, N - i*MAX_ALLOCA_SF32);
    /*
    * Compute log2(x) in form of y0+y1
    */ 
    pX = (xb_vecN_2xf32 *)(x + i*MAX_ALLOCA_SF32);
    SCR0_wr = (xb_vecN_2xf32 *)scratch0;
    SCR1_wr = (xb_vecN_2xf32 *)scratch1;
    TBL = c_tbl;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH-1)); n++)
    {
      xb_vecN_2x32v e_reg, temp0, temp1, x_i32;
      xb_vecNx16 x_i16, h;
      xb_vecN_2xf32 t0, t1;
      vboolN_2 b_denorm, b_small;
      BBE_LVN_2XF32_IP(xin, pX, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

      /* Take mantissa and exponent of input */
      xin = BBE_ABSN_2XF32(xin);
      temp0 = BBE_LSN_2X32_I((const xb_int32v *)TBL, 0* sizeof(float32_t));
      h = BBE_REPNX16(BBE_MOVNX16_FROMN_2X32(temp0), 0);/* load exp=-126 */
      temp0 = BBE_LSN_2X32_I((const xb_int32v *)TBL, 1 * sizeof(float32_t)); temp0 = BBE_REPN_2X32(temp0, 0);/* load exp=-150 */
      t0 = BBE_LSN_2XF32_I((const xtfloat *)TBL, 2 * sizeof(float32_t));  t0 = BBE_REPN_2XF32(t0, 0);/* load max denormalized floating-point number */
      t1 = BBE_LSN_2XF32_I((const xtfloat *)TBL, 3 * sizeof(float32_t)); t1 = BBE_REPN_2XF32(t1, 0);/* load 2^24 */
      /* process denormalized values */
      b_denorm = BBE_OLEN_2XF32(xin, t0);
      t1 = BBE_MULN_2XF32(xin, t1);
      xin = BBE_MOVN_2XF32T(t1, xin, b_denorm);
      h = (BBE_MOVNX16T(BBE_MOVNX16_FROMN_2X32(temp0), h, BBE_MOVN_FROMN_2(b_denorm)));
      /* extract exponent */
      x_i16 = BBE_MOVNX16_FROMN_2XF32(xin);
      x_i16 = BBE_SRLINX16(x_i16, 7);/* 23 - 16 */
      x_i16 = BBE_SELNX16I(x_i16, x_i16, BBE_SELI_INTERLEAVE_1_ODD);
      h = BBE_ADDNX16(h, x_i16);

      /* extract mantissa */
      temp0 = BBE_MOVN_2X32_FROMN_2XF32(t0);/* load mantissa mask */
      temp1 = BBE_LSN_2X32_I((const xb_int32v *)TBL, 5 * sizeof(float32_t));
      temp1 = BBE_REPN_2X32(temp1, 0);/* load 0.5 */
      x_i32 = BBE_MOVN_2X32_FROMN_2XF32(xin);
      x_i32 = BBE_ANDN_2X32(x_i32, temp0);
      x_i32 = BBE_ORN_2X32(x_i32, temp1);
      xin = BBE_MOVN_2XF32_FROMN_2X32(x_i32);
      /* adjust the mantissa to range [ sqrt(0.5) ; sqrt(2.0) ) */
      t0 = BBE_LSN_2XF32_I((const xtfloat *)TBL, 4 * sizeof(float32_t)); t0 = BBE_REPN_2XF32(t0, 0);/* load sqrt(0.5) */
      b_small = BBE_OLTN_2XF32(xin, t0);
      t0 = BBE_ADDN_2XF32(xin, xin);
      BBE_SUBNX16T(h, h, 1, BBE_MOVN_FROMN_2(b_small));
      xin = BBE_MOVN_2XF32T(t0, xin, b_small);
      one = BBE_CONSTN_2XF32(1);
      xin = BBE_SUBN_2XF32(one, xin);
      {
        xb_vecNx16 tmp;
        tmp = BBE_SRAINX16(h, 15);
        h = BBE_SELNX16I(tmp, h, BBE_SELI_INTERLEAVE_1_EVEN);
        e_reg = BBE_MOVN_2X32_FROMNX16(h);
      }
      ef_reg = BBE_FLOATN_2X32(e_reg, 0);

      /* save results to scratch buffer */
      BBE_SVN_2XF32_IP(xin, SCR0_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(ef_reg, SCR1_wr, 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");

    pY = (const xb_vecN_2xf32 *)(y + i*MAX_ALLOCA_SF32);
    SCR0_rd = (const xb_vecN_2xf32 *)scratch0;
    SCR1_rd = (const xb_vecN_2xf32 *)scratch1;
    SCR0_wr = (xb_vecN_2xf32 *)scratch0;
    SCR1_wr = (xb_vecN_2xf32 *)scratch1;
    // TBL = c_tbl + 5;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
      xb_vecN_2xf32 p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, px0, px1;
      xb_vecN_2xf32 t0, t1, t2, w0, w1;
      BBE_LVN_2XF32_IP(xin, SCR0_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(ef_reg, SCR1_rd, 2 * BBE_SIMD_WIDTH);
      /* evaluate polynomial approximation */
      /* Load table of coefficients */
      px0 = BBE_LVN_2XF32_I(TBL_LOG2, 0);
      px1 = BBE_LVN_2XF32_I(TBL_LOG2, 2*BBE_SIMD_WIDTH);
      p0 = BBE_REPN_2XF32(px0, 0);
      p1 = BBE_REPN_2XF32(px0, 1);
      p2 = BBE_REPN_2XF32(px0, 2);
      p3 = BBE_REPN_2XF32(px0, 3);
      p4 = BBE_REPN_2XF32(px0, 4);
      p5 = BBE_REPN_2XF32(px0, 5);
      p6 = BBE_REPN_2XF32(px0, 6);
      p7 = BBE_REPN_2XF32(px0, 7);
      p8 = BBE_REPN_2XF32(px1, 0);
      p9 = BBE_REPN_2XF32(px1, 1);
      p10 = BBE_REPN_2XF32(px1, 2);
      p11 = BBE_REPN_2XF32(px1, 3);
      p12 = BBE_REPN_2XF32(px1, 4);
      p13 = BBE_REPN_2XF32(px1, 5);
      #if 1
      BBE_MULAN_2XF32(p1, xin, p0);
      BBE_MULAN_2XF32(p2, xin, p1);
      BBE_MULAN_2XF32(p3, xin, p2);
      BBE_MULAN_2XF32(p4, xin, p3);
      BBE_MULAN_2XF32(p5, xin, p4);
      BBE_MULAN_2XF32(p6, xin, p5);
      BBE_MULAN_2XF32(p7, xin, p6);
      BBE_MULAN_2XF32(p8, xin, p7);
      BBE_MULAN_2XF32(p9, xin, p8);
      t2=p9;
      #else
      {
        xb_vecN_2xf32 xin2;
        xin2 = BBE_MULN_2XF32(xin, xin);
        BBE_MULAN_2XF32(p2, xin2, p0);
        BBE_MULAN_2XF32(p4, xin2, p2);
        BBE_MULAN_2XF32(p6, xin2, p4);
        BBE_MULAN_2XF32(p8, xin2, p6);
        BBE_MULAN_2XF32(p3, xin2, p1);
        BBE_MULAN_2XF32(p5, xin2, p3);
        BBE_MULAN_2XF32(p7, xin2, p5);
        BBE_MULAN_2XF32(p9, xin2, p7);
        BBE_MULAN_2XF32(p9, xin, p8);
        t2 = p9;
      }

      #endif
      /* next coefficients are computed in extended precision */
      t0 = BBE_MULN_2XF32(xin, t2);
      t1 = t0; BBE_MULSN_2XF32(t1, xin, t2);
      w0 = BBE_ADDN_2XF32(t0, p10);
      w1 = BBE_SUBN_2XF32(w0, p10);
      w1 = BBE_SUBN_2XF32(t0, w1);
      w1 = BBE_SUBN_2XF32(w1, t1);
      t0 = w0; t1 = w1;
      w0 = BBE_MULN_2XF32(xin, t0);
      w1 = w0; BBE_MULSN_2XF32(w1, xin, t0);
      t0 = w0;
      BBE_MULSN_2XF32(w1, t1, xin); t1 = w1;
      w0 = BBE_ADDN_2XF32(t0, p11);
      w1 = BBE_SUBN_2XF32(w0, p11);
      w1 = BBE_SUBN_2XF32(t0, w1);
      w1 = BBE_SUBN_2XF32(w1, t1);
      t0 = w0; t1 = w1;
      xin = BBE_NEGN_2XF32(xin);
      w0 = BBE_MULN_2XF32(t0, xin);
      w1 = w0; BBE_MULSN_2XF32(w1, t0, xin);
      t0 = w0;
      BBE_MULSN_2XF32(w1, t1, xin); t1 = w1;
      /* multiply by log2(e) */
      w0 = BBE_MULN_2XF32(t0, p12);
      w1 = w0;
      BBE_MULSN_2XF32(w1, t0, p12);
      BBE_MULAN_2XF32(w1, t1, p12);
      BBE_MULSN_2XF32(w1, t0, p13);
      t0 = w0; t1 = w1;

      /* add exponent */
  
      w0 = BBE_ADDN_2XF32(t0, ef_reg);
      w1 = BBE_SUBN_2XF32(w0, ef_reg);
      w1 = BBE_SUBN_2XF32(t0, w1);
      t1 = BBE_SUBN_2XF32(w1, t1);
      t0 = w0;
      /* save results to scratch buffer */
      BBE_SVN_2XF32_IP(t0, SCR0_wr, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(t1, SCR1_wr, 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");
    /*
    * Compute 2^(y*log2(x))
    */
    pY = (const xb_vecN_2xf32 *)(y + i*MAX_ALLOCA_SF32);
    SCR0_rd = (const xb_vecN_2xf32 *)scratch0;
    SCR1_rd = (const xb_vecN_2xf32 *)scratch1;
    SCR0_wr = (xb_vecN_2xf32 *)scratch0;
    SCR1_wr = (xb_vecN_2xf32 *)scratch1;
   // TBL = c_tbl + 5;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
      xb_vecN_2xf32 yin, xy, dxy;
      xb_vecN_2xf32 t0, t1;
      xb_vecN_2xf32 zero;
      /* Load computed log2(x), add normalization factor and combine with y*/
      BBE_LVN_2XF32_IP(t0, SCR0_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(t1, SCR1_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(yin, pY, 2 * BBE_SIMD_WIDTH);
      zero = BBE_ZERON_2XF32();
      /* compute y*log2(x) and separate result into integer and fractional parts */
      xy = BBE_FIROUNDN_2XF32(BBE_MULN_2XF32(yin, t0));
      dxy = BBE_NEGN_2XF32(xy);
      BBE_MULAN_2XF32(dxy, yin, t0);
      BBE_MULAN_2XF32(dxy, yin, t1);
      dxy = BBE_MINN_2XF32(dxy, 1.0f);
      dxy = BBE_MAXN_2XF32(dxy, -1.0f);

      xy = BBE_MINN_2XF32(xy, 129.f);
      xy = BBE_MAXN_2XF32(xy, -151.f);
      BBE_SVN_2XF32_IP(dxy, SCR0_wr, 2 * BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(xy, SCR1_wr, 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");
    SCR0_rd = (const xb_vecN_2xf32 *)scratch0;
    SCR1_rd = (const xb_vecN_2xf32 *)scratch1;
    SCR0_wr = (xb_vecN_2xf32 *)scratch0;
    SCR1_wr = (xb_vecN_2xf32 *)scratch1;
    // TBL = c_tbl + 5;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
      xb_vecN_2xf32 xy, dxy, c0, c1;
      xb_vecN_2xf32 px, p0, p1, p2, p3, p4 ,p5 , p6, zout;
      BBE_LVN_2XF32_IP(dxy, SCR0_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(xy, SCR1_rd, 2 * BBE_SIMD_WIDTH);
      /* compute 2^fract */
      px = BBE_LVN_2XF32_I(TBL_POW2, 0);
      p0 = BBE_REPN_2XF32(px, 0);
      p1 = BBE_REPN_2XF32(px, 1);
      p2 = BBE_REPN_2XF32(px, 2);
      p3 = BBE_REPN_2XF32(px, 3);
      p4 = BBE_REPN_2XF32(px, 4);
      p5 = BBE_REPN_2XF32(px, 5);
      p6 = BBE_REPN_2XF32(px, 6);
      /* NOTE: do not change the order of computations and way of polynomial decomposition ! */
      BBE_MULAN_2XF32(p1, dxy, p0);
      BBE_MULAN_2XF32(p2, dxy, p1);
      BBE_MULAN_2XF32(p3, dxy, p2);
      BBE_MULAN_2XF32(p4, dxy, p3);
      BBE_MULAN_2XF32(p5, dxy, p4);
      BBE_MULAN_2XF32(p6, dxy, p5);
      zout = p6;
      /* apply integer part */

      {
        xb_vecNx16 tmp, e1, e2;
        tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(xy, 0));
        tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254));
        e1 = BBE_SRLINX16(tmp, 1);
        e2 = BBE_SUBNX16(tmp, e1);
        e1 = BBE_SLLINX16(e1, 7);
        e2 = BBE_SLLINX16(e2, 7);
        e1 = BBE_SELNX16I(e1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
        e2 = BBE_SELNX16I(e2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
        c0 = BBE_MOVN_2XF32_FROMNX16(e2);
        c1 = BBE_MOVN_2XF32_FROMNX16(e1);
      }
      zout = BBE_MULN_2XF32(zout, c1);
      zout = BBE_MULN_2XF32(zout, c0);
      BBE_SVN_2XF32_IP(zout, SCR0_wr, 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");
    /*
    * Consider special cases, take decision on output values, and update the
    * error states.
    */


    pX = (const xb_vecN_2xf32*)(x + i*MAX_ALLOCA_SF32);
    pY = (const xb_vecN_2xf32*)(y + i*MAX_ALLOCA_SF32);
    SCR0_rd = (const xb_vecN_2xf32*)scratch0;

    SCR1_wr = ( xb_vecN_2xf32*)scratch1;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
      xb_vecN_2xf32 xin, yin, zout, zero, xabs;
      xb_vecN_2x32v x_i32, y_i32, temp0, temp1, zero_i32;
      xb_vecN_2xf32 t0, spec, half;
      vboolN_2 b_zero, b_one, b_Inf, b_NaN1, b_NaN2, b_notspec, b_xeqz;
      vboolN_2 b_sx, b_sy, b_yint, b_yeqz, b_tmp0, b_xinf, b_xeq1, b_yinf, b_yodd;
      BBE_LVN_2XF32_IP(xin, pX, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(yin, pY, 2 * BBE_SIMD_WIDTH);
      half = BBE_CONSTN_2XF32(3);
      zero = BBE_CONSTN_2XF32(0);
      zero_i32 = BBE_ZERON_2X32();
      x_i32 = BBE_MOVN_2X32_FROMN_2XF32(xin);
      y_i32 = BBE_MOVN_2X32_FROMN_2XF32(yin);

      /* Take sign of x and y */

      maskx = BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(0x80000000));
      temp0 = BBE_ANDN_2X32(maskx, BBE_MOVN_2X32_FROMN_2XF32(xin));
      b_sx = BBE_EQN_2X32(temp0, maskx);
      temp0 = BBE_ANDN_2X32(maskx, BBE_MOVN_2X32_FROMN_2XF32(yin));
      b_sy = BBE_EQN_2X32(temp0, maskx);

      xabs = BBE_ABSN_2XF32(xin);

      /* check if y is integer */
      t0 = BBE_FITRUNCN_2XF32(yin);
      b_yint = BBE_OEQN_2XF32(t0, yin);

      /* check if y is odd */
      temp0 = BBE_TRUNCN_2XF32(yin, 0);
      b_tmp0 = BBE_NEQN_2X32(temp0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(MAX_INT32)));
      temp1 = BBE_ANDN_2X32(temp0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(1)));
      temp0 = BBE_MOVN_2X32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2X32(temp1), BBE_MOVNX16_FROMN_2X32(temp0), BBE_MOVN_FROMN_2(b_tmp0)));
      b_yodd = BBE_EQN_2X32(temp0, BBE_MOVN_2X32_FROMNX16(BBE_MOVVA16C(1)));
      b_yodd = BBE_ANDBN_2(b_yodd, b_yint);

      /* process special numbers */
      b_NaN1 = BBE_ANDNOTBN_2(b_sx, b_yint);                        /* x<0 && y is not an integer --> z=NaN                */
      b_NaN2 = BBE_UNN_2XF32(xabs, yin);                            /* isnan(x)||isnan(y) --> z=NaN                        */
      b_yeqz = BBE_OEQN_2XF32(zero, yin);                           /* |y|==0                                              */
      b_xeq1 = BBE_OEQN_2XF32(one, xabs);                           /* |x|==1                                              */
      b_yinf = BBE_OEQN_2XF32(BBE_ABSN_2XF32(yin), plusInff.f);  /* |y| == Inf                                          */
      b_one = BBE_NOTBN_2(b_sx);
      b_one = BBE_ANDBN_2(b_xeq1, BBE_ORBN_2(b_yinf, b_one));  /* |x|==1 && ( |y|==Inf || x>0 )                       */
      b_one = BBE_ORBN_2(b_one, b_yeqz);                            /* ( |x|==1 && ( |y|==Inf || x>0 ) ) || y==0 --> z=1.0 */

      b_xeqz = BBE_OEQN_2XF32(xabs, zero);                          /* x==0                    */
      b_xinf = BBE_OEQN_2XF32(xabs, plusInff.f);                    /* x==INF                  */
      b_zero = BBE_ANDNOTBN_2(b_xeqz, b_sy);                       /* x==0   && y>0 --> z=0.0 */
      b_Inf = BBE_ANDNOTBN_2(b_xinf, b_sy);                         /* x==INF && y>0 --> z=INF */

      b_zero = BBE_ORBN_2(b_zero, BBE_ANDBN_2(b_xinf, b_sy));        /* x==0   && y<0 --> z=INF */
      b_Inf = BBE_ORBN_2(b_Inf, BBE_ANDBN_2(b_xeqz, b_sy));         /* x==INF && y<0 --> z=0.0 */

      b_sx = BBE_ANDBN_2(b_sx, b_yodd);

      /* Save special numbers and mask for special numbers */
      spec = BBE_MOVN_2XF32T(qNaNf.f, half, b_NaN1);
      BBE_CONSTN_2XF32T(spec, 0, b_zero);
      spec = BBE_MOVN_2XF32T(plusInff.f, spec, b_Inf);
      spec = BBE_MOVN_2XF32T(qNaNf.f, spec, b_NaN2);
      BBE_CONSTN_2XF32T(spec, 1, b_one);

      b_notspec = BBE_OEQN_2XF32(spec, half);

      /* load precomputed values of 2^(y*log2(x)) from the scratch */
      BBE_LVN_2XF32_IP(t0, SCR0_rd, 2 * BBE_SIMD_WIDTH);

      /* Replace result with special numbers if needed */
      zout = BBE_MOVN_2XF32T(t0, spec, b_notspec);

      /* Restore sign */
      BBE_NEGN_2XF32T(zout, zout, b_sx);
      BBE_SVN_2XF32_IP(zout, pZ, 2 * BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(spec, SCR1_wr, 2 * BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder");
    pX = (const xb_vecN_2xf32*)(x + i*MAX_ALLOCA_SF32);
    pY = (const xb_vecN_2xf32*)(y + i*MAX_ALLOCA_SF32);
    SCR0_rd = (const xb_vecN_2xf32*)(z + i*MAX_ALLOCA_SF32);
    //SCR0_rd = (const xb_vecN_2xf32*)scratch0;
    SCR1_rd = (const xb_vecN_2xf32*)(const xb_vecN_2xf32*)scratch1;
    for (n = 0; n < (_N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
      xb_vecN_2xf32 xin, yin, zout, zero, xabs;
      xb_vecN_2xf32 spec, half;
      vboolN_2 b_NaN2, b_notspec, b_xeqz, b_xeqz_yltz;
      vboolN_2 b_xinf, b_yinf;
      BBE_LVN_2XF32_IP(xin, pX, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(yin, pY, 2 * BBE_SIMD_WIDTH);
      half = BBE_CONSTN_2XF32(3);
      zero = BBE_CONSTN_2XF32(0);

      /* Take sign of x and y */
      flags_x = BBE_CLSFYN_2XF32(xin);
      flags_y = BBE_CLSFYN_2XF32(yin);

      xabs = BBE_ABSN_2XF32(xin);


      /* process special numbers */
      b_NaN2 = BBE_UNN_2XF32(xabs, yin);                            /* isnan(x)||isnan(y) --> z=NaN                        */
      b_yinf = BBE_OEQN_2XF32(BBE_ABSN_2XF32(yin), plusInff.f);  /* |y| == Inf                                          */

      b_xeqz = BBE_OEQN_2XF32(xabs, zero);                          /* x==0                    */
      b_xinf = BBE_OEQN_2XF32(xabs, plusInff.f);                    /* x==INF                  */

      /* load precomputed values of 2^(y*log2(x)) from the scratch */
      BBE_LVN_2XF32_IP(zout, SCR0_rd, 2 * BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(spec, SCR1_rd, 2 * BBE_SIMD_WIDTH);
      b_notspec = BBE_OEQN_2XF32(spec, half);

      /*
      * Perform additional analysis of input data for Error Handling.
      */

      /* check x and y for signalling NaN. */

      flags_xy = BBE_ORN_2X32U(flags_x, flags_y);
      maskx = BBE_MOVN_2X32U_FROMNX16(BBE_MOVVA16C(0x40));
      v_snan = BBE_ANDN_2X32U(flags_xy, maskx); /* sNaN */
      b_snan = BBE_NEQN_2X32U(v_snan, BBE_ZERON_2X32U());

      b_spec_nan = BBE_UNN_2XF32(spec, spec);

      /* EDOM conditions:
      * A) x<0 && x is finite && y is not an integer && y is finite
      * B) x is NaN and y!=0
      * C) y is NaN and x!=1
      * We use that (A or B or C) if z==NaN. */
      b_edom = BBE_ORBN_2(b_edom, b_spec_nan);

      b_xeqz_yltz = BBE_ANDBN_2(b_xeqz, BBE_OLTN_2XF32(yin, zero));
      /* Check if x, y are finite and z is infinite. */
      b_xyinf = BBE_ORBN_2(b_xinf, b_yinf);
      b_zinf = BBE_OEQN_2XF32(zout, plusInff.f);
      b_xyfin_zinf = BBE_ANDNOTBN_2(b_zinf, b_xyinf);

      /* (x==0)&&(y<0) || x is finite && y is finite && z is infinite -> set ERANGE
      * First term is shared with FE_DIVBYZERO condition! */
      b_erange = BBE_ORBN_2(b_erange, b_xyfin_zinf);

      /* x<0 && x is finite && y is not an integer && y is finite --> raise "invalid" exception */
      b_fe_inv = BBE_ORBN_2(b_fe_inv, BBE_ANDNOTBN_2(b_spec_nan, b_NaN2));
      /* "invalid" exception should be also raised if either input is a signalling NaN. */
      b_fe_inv = BBE_ORBN_2(b_fe_inv, b_snan);
      /* x==0 && y<0 --> raise "divide-by-zero" exception */
      b_fe_divz = BBE_ORBN_2(b_fe_divz, b_xeqz_yltz);

      /* x!=0 && x is finite && y is finite && z is infinite --> raise "overflow" exception */
      b_fe_ovfl = BBE_ORBN_2(b_fe_ovfl, BBE_ANDBN_2(b_xyfin_zinf, b_notspec));


    }
  }

   {
     /* FE_DIVBYZERO is a subset of ERANGE, thus we use the former to trace the
     * condition (x==0)&&(y<0).  */
     b_erange = BBE_ORBN_2(b_erange, b_fe_divz);
     v_edom = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_edom);
     v_erange = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(1), BBE_ZERON_2XF32(), b_erange);
     v_fe_inv = BBE_MOVN_2XF32T(FE_INVALID, BBE_ZERON_2XF32(), b_fe_inv);
     v_fe_ovfl = BBE_MOVN_2XF32T(FE_OVERFLOW, BBE_ZERON_2XF32(), b_fe_ovfl);
     v_fe_divz = BBE_MOVN_2XF32T(FE_DIVBYZERO, BBE_ZERON_2XF32(), b_fe_divz);

     er_edom = BBE_RMAXNUMN_2XF32(v_edom);
     er_erange = BBE_RMAXNUMN_2XF32(v_erange);
     fe_inv = BBE_RMAXNUMN_2XF32(v_fe_inv);
     fe_ovfl = BBE_RMAXNUMN_2XF32(v_fe_ovfl);
     fe_divz = BBE_RMAXNUMN_2XF32(v_fe_divz);


     if (0 != er_erange) { __Pragma("frequency_hint never"); errno = ERANGE; };
     /* EDOM takes precedence over ERANGE! */
     if (0 != er_edom) { __Pragma("frequency_hint never"); errno = EDOM; };

     /* Restore exception enable flags and status flags, suppress undesired status flags. */
     __fesetenv(&fenv);
     __feraiseexcept(fe_inv | fe_ovfl | fe_divz);
  }
} /* vpowf() */
#endif

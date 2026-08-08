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
    Polar to Cartesian Conversion
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
DISCARD_FUN(complex_float, spolarf, ( float32_t r, float32_t a ))
#else
#include "sinf_tbl.h"
#include "inv2pif_tbl.h"

/*-------------------------------------------------------------------------
Polar to Cartesian Conversion 

Description: These functions convert input pairs of magnitude and phase
values into Cartesian coordiantes of respective points on the complex plane. 

Representation:
vpolar           16-bit signed fixed-point format
                 Input phase data phase[] are Q15 angular values
                 normalized by pi.
                 Input magnitude data x[] are of abn arbitrary fixed-
                 point format Qx.
                 Fixed-point format for resulting coordinates on the 
                 complex plane z[] is Qx-sh.
vpolarf,spolarf  IEEE-754 Std. single precision floating-point format
                 Input phase data phase[] are in radians.

Accuracy:
9.2e-5 (3 in Q15) - for vpolar
2 ULP - vpolarf, spolarf
3 ULP - vfastpolarf

Note for vpolarf,spolarf:
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the conversion result is (0,0).

Input domain for 'fast' version vfastpolarf()
|a|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
r[N]       Magnitude data
a[N]       Phase data
sh         Right bit shift amount, 0..15 (vpolar)
N          Length of vectors
Output:
z[N]       Points on the complex plane

Restrictions:
z,x,phase  Aligned on 32-byte boundary
z,x,phase  Must not overlap
N          Multiple of 16 (vpolar) or 8 (vpolarf,vfastpolarf)
-------------------------------------------------------------------------*/

complex_float spolarf ( float32_t r, float32_t a )
{
  float32_t phase, phase2, ys, yc, t, jf;
  int ss, sc, ji;
  float32_t y_re, y_im;
  union
  {
    complex_float ALIGN(8) z;
    struct { float32_t re, im; }s;
  }z;
  int nondom;

  phase = a;
  /* find sine and cosine of the phase */
  phase2 = (float32_t)XT_ABS_S(phase);
  nondom = phase2 > sinf_maxval.f;
  /* argument reduction */
  jf = XT_MUL_S(phase, inv2pif.f);
  jf = XT_FIROUND_S(jf);
  ji = XT_TRUNC_S(jf, 0);
  XT_MSUB_S(phase, pi2fc[0].f, jf);
  XT_MSUB_S(phase, pi2fc[1].f, jf);
  XT_MSUB_S(phase, pi2fc[2].f, jf);

  /* adjust signs */
  ss = (((ji) >> 1) & 1);
  sc = ((ji + 1) >> 1) & 1;
  /* compute sine/cosine via minmax polynomial  */
  phase2 = phase*phase;
  ys = polysinf_tbl[0].f;
  t = polysinf_tbl[1].f; XT_MADD_S(t, ys, phase2); ys = t;
  t = polysinf_tbl[2].f; XT_MADD_S(t, ys, phase2); ys = t;
  ys = ys*phase2;
  t = phase; XT_MADD_S(t, ys, phase); ys = t;
  yc = polycosf_tbl[0].f;
  t = polycosf_tbl[1].f; XT_MADD_S(t, yc, phase2); yc = t;
  t = polycosf_tbl[2].f; XT_MADD_S(t, yc, phase2); yc = t;

  t = XT_CONST_S(1);
  XT_MADD_S(t, yc, phase2); yc = t;
  /* select real/imag */
  y_re = ys; XT_MOVEQZ_S(y_re, yc, (ji & 1));
  y_im = yc; XT_MOVEQZ_S(y_im, ys, (ji & 1));
  /* apply the sign */
  XT_MOVNEZ_S(y_re, -y_re, sc);
  XT_MOVNEZ_S(y_im, -y_im, ss);

  XT_MOVNEZ_S(y_re, 0.f, nondom);
  XT_MOVNEZ_S(y_im, 0.f, nondom);
  y_re *= r;
  y_im *= r;
  /* combine y_re,y_im */
  z.s.re = y_re;
  z.s.im = y_im;
  return z.z;
} /* spolarf() */

#endif/* !HAVE_VFPU */

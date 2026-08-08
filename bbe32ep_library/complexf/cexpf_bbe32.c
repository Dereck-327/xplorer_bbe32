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
    Complex Exponential
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
DISCARD_FUN(complex_float, sexpf, ( float32_t phase ))
#else
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"

/*-------------------------------------------------------------------------
Complex Exponential

Description: These functions compute complex number lying on the unit 
circle from input phase data.

Representation:
cbexp         16-bit signed fixed-point format
              Input data phase[] are Q15 angular values normalized by pi.
              Output data z[] format is Q(15-sh), where sh is the shift
              control argument
cbexpf,sexpf  IEEE-754 Std. single precision floating-point format
              Input phase data phase[] are in radians. See the Note.
              
Accuracy:
3.7e-4 (12 in Q15) for cbexp
2 ULP for cbexpf, sexpf
3 ULP for cbfastexpf

Note for (non-fast versions):
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the respective result is 0+0j.

Input domain for 'fast' version cfastbexpf()
|phase|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase[N]  Phase values
sh        (cbexp only) left bit shift amount, [0..15]
N         Length of vectors
Output:
z[N]      Results

Restrictions:
z,phase   Aligned on 32-byte boundary
z,phase   Must not overlap
N         Multiple of 16 (cbexp) or 8 (cbexpf,cbfastexpf)
-------------------------------------------------------------------------*/

complex_float sexpf ( float32_t phase )
{
  float32_t phase2, ys, yc, t, jf;
  int ss, sc, ji;
  float32_t y_re, y_im;
  int nondom;
  union
  {
    complex_float ALIGN(8) z;
    struct { float32_t re, im; }s;
  }z;
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
  /* combine y_re,y_im */
  z.s.re = y_re;
  z.s.im = y_im;
  return z.z;
} /* sexpf() */

#endif/* !HAVE_VFPU */

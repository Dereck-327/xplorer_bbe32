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
    Complex Exponential with Linearly Evolving Phase
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
DISCARD_FUN(float32_t, cexplinf, ( complex_float * restrict z, float32_t phase0, float32_t delta ))
#else
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"

/*-------------------------------------------------------------------------
Complex Exponential with Linearly Evolving Phase

Description: These functions compute a series of complex numbers lying on
the unit circle, with the phase angle being incremented for each successive
complex number.

Representation:
cbexplin,           16-bit signed fixed-point format
cbexplin_fast       Initial phase and phase increment (phase0 and delta) are 
                    Q15 angular values normalized by pi.
                    Output data z[] format is Q15.
cbexplinf,          IEEE-754 Std. single precision floating-point format
cexplinf            Initial phase and phase increment (phase0 and delta) are
                    in radians. See Note 3.

Accuracy:
cbexplin            Absolute error for re/im components does not exceed
                    3.7e-4 (12 in Q15)
cbexplin_fast       Absolute error for re/im components for n-th result is
                    ~(floor(n/8)+1)*3.7e-4, n=0..N-1
cbexplinf,cexplinf  2 ULP 
cbfastexplinf       3 ULP

Notes:
1. cbexplin() function computes each complex value from a full-precision
   linearly evolving angle, thus the function is able to generate long data
   sequences with no risk of phase errors accumulation.
2. cbexplin_fast() function is approximately twice faster but less accurate 
   because it performs successive rotations by small angle on each iteration.
   Avoid using this function for long sequences or tiny phase increments.
3. For the non-fast floating-point functions (cbexplinf, cexplinf) , the 
   phase (initial and updated) should belong to the range [-102940.0, 
   102940.0], otherwise the respective resultis 0+0j.
4. Scalar function cexplinf uses the initial phase angle (phase0) to compute 
   a single complex number, and returns the once incremented phase value. 

Input domain for 'fast' version (cbfastexplinf)
|phase|<804.2477
for all initial and updated phases 
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase0  Initial phase angle (phase of the first complex number in a series)
delta   Phase increment step
N       Length of output vector
Output:
z[N]    Complex numbers
Returned Value:
        Updated phase, phase0+N*delta

Restrictions:
z       Aligned on 32-byte boundary
z       Must not overlap
N       Multiple of 16 (cbexplln, cbexplln_fast) or 8 (cbexplinf,cbfastexplinf)
-------------------------------------------------------------------------*/
float32_t cexplinf ( complex_float * restrict z, float32_t phase0, float32_t delta )
{
  float32_t ph, phase2, ys, yc, t, jf;
  int ss, sc, ji;
  float32_t y_re, y_im;
  int nondom;
  union
  {
    complex_float ALIGN(8) z;
    struct { float32_t re, im; }s;
  }cplx;

  NASSERT(z);
  ph = phase0;
  /* find sine and cosine of the phase */
  phase2 = (float32_t)XT_ABS_S(ph);
  nondom = phase2 > sinf_maxval.f;
  /* argument reduction */
  jf = XT_MUL_S(ph, inv2pif.f);
  jf = XT_FIROUND_S(jf);
  ji = XT_TRUNC_S(jf, 0);
  XT_MSUB_S(ph, pi2fc[0].f, jf);
  XT_MSUB_S(ph, pi2fc[1].f, jf);
  XT_MSUB_S(ph, pi2fc[2].f, jf);

  /* adjust signs */
  ss = (((ji) >> 1) & 1);
  sc = ((ji + 1) >> 1) & 1;
  /* compute sine/cosine via minmax polynomial  */
  phase2 = ph*ph;
  ys = polysinf_tbl[0].f;
  t = polysinf_tbl[1].f; XT_MADD_S(t, ys, phase2); ys = t;
  t = polysinf_tbl[2].f; XT_MADD_S(t, ys, phase2); ys = t;
  ys = ys*phase2;
  t = ph; XT_MADD_S(t, ys, ph); ys = t;
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
  cplx.s.re = y_re;
  cplx.s.im = y_im;
  z[0] = cplx.z;

  phase0 += delta;
  return phase0;
} /* cexplinf() */

#endif/* !HAVE_VFPU */

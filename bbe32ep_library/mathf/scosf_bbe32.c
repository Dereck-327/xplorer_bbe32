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
    Sine/Cosine
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
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Sine/Cosine 

Description: These functions compute sine or cosine of input data

Representation:
vsine,vcos,    16-bit signed fixed-point format Q15 for input/output data
ssine,scos     It is assumed that input angular values are normalized by pi.
               That is, Fixed-point functions actually compute sin(pi*x) or
               cos(pi*x)
vsinef,vcosf,  IEEE-754 Std. single precision floating-point format for
ssinef,scosf   input/output data. Input data are treated as angular values
               specified in radians. Floating-point functions limit the
               rangw of allowable input values, see note 2.

Accuracy:
2 LSB - vsine(),vcos(), ssine(), scos(),
2 ULP - vsinef(), vcosf(), ssinef(), scosf(),
3 ULP - vfastsinef(),vfastcosf()

Notes for non-fast versions:
1. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Floating-point functions require that input value belongs to the 
   closed range [-102940.0,102940.0], otherwise the respective result
   is NaN.

Input domain for 'fast' versions vfastsinef(),vfastcosf()
|x|<804.2477
The output value is not defined outside of this range or accuracy is degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vsine,vcos) or 8 (vsinef,vcosf,vfastsinef,vfastcosf)
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t scosf ( float32_t x )
{
  vbool1 b_nan, b_outl, b_sn, b_neg;
  xb_int32v ji;
  xtfloat jf, x2, t, ys, yc, y;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  /* Take the absolute value. */
  x = XT_ABS_S(x);

  /* Check if the input value is a NaN, or does not belong to the allowed range, */
  b_nan = XT_UN_S(x, x);
  b_outl = XT_ULT_S(sinf_maxval.f, x);
  if (BBE_MOVAB1(b_outl))
  {
    __Pragma("frequency_hint never");
    if (BBE_MOVAB1(b_nan))
    {
      errno = EDOM;
    }

    return (qNaNf.f);
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  /* Argument reduction. */
  jf = XT_MUL_S(x, inv2pif.f);
  jf = XT_FIROUND_S(jf);
  ji = XT_TRUNC_S(jf, 0);
  XT_MSUB_S(x, pi2fc[0].f, jf);
  XT_MSUB_S(x, pi2fc[1].f, jf);
  XT_MSUB_S(x, pi2fc[2].f, jf);
  x2 = XT_MUL_S(x, x);

  /* Compute sine via minmax polynomial  */
  t = polysinf_tbl[0].f;                         ys = t;
  t = polysinf_tbl[1].f; XT_MADD_S(t, ys, x2); ys = t;
  t = polysinf_tbl[2].f; XT_MADD_S(t, ys, x2); ys = t;
  t = XT_MUL_S(ys, x2);
  ys = x; XT_MADD_S(ys, t, x);;

  /* Compute cosine via minmax polynomial  */
  t = polycosf_tbl[0].f;                         yc = t;
  t = polycosf_tbl[1].f; XT_MADD_S(t, yc, x2); yc = t;
  t = polycosf_tbl[2].f; XT_MADD_S(t, yc, x2); yc = t;
  t = XT_CONST_S(1);     XT_MADD_S(t, yc, x2); yc = t;

  /* Determine the sign and sin/cos selector from the input range. */
  {
    xb_vecNx16 h0, h1, h2;
    vboolN s;
    vboolN bsc, bneg;
    h0 = BBE_MOVNX16_FROM32(ji);

    h1 = BBE_SLLINX16(h0, 14);
    {
      xb_vecNx16 tt;
      tt =0x4000;
      h1 = BBE_ADDNX16(h1, tt);
    }

    h2 = BBE_SLLINX16(h0, 15);

    h0 = BBE_SELNX16I(h2, h1, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    s = BBE_LTNX16(h0, BBE_ZERONX16());
    BBE_EXTRACTB(bsc, bneg, s);
      b_sn = BBE_MOV1_FROMN(bsc);
      b_neg = BBE_MOV1_FROMN(bneg);
  }
  /* Select sine or cosine. */
  y = yc; XT_MOVT_S(y, ys, b_sn);
  /* Adjust the sign. */
  ys = XT_NEG_S(y); XT_MOVT_S(y, ys, b_neg);
  BBE_MOVSCFV(SCF);

  return (y);

} /* scosf() */
#else
DISCARD_FUN(float32_t,scosf,( float32_t x ));
#endif

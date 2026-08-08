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
    Tangent
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
/* Angular argument reduction constants. */
#include "sinf_tbl.h"
/* tan/cotan approximation polynomial coeffs. */
#include "tanf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Tangent 

Description: These functions compute tangent of input data

Representation:
vtan,stan    Signed fixed-point format
             Input data are 16-bit Q15 angular values normalized by pi,
             i.e. fixed-point functions actually compute tan(pi*x).
             Output data are 32-bit Q16.15 values.
vtanf,stanf  IEEE-754 Std. single precision floating-point format for
             input/output data. Input data are treated as angular values
             specified in radians. Floating-point functions limit the
             rangw of allowable input values, see note 3.

Accuracy:
For the fixed-point functions, accuracy depends on the input value x,
as shown in the table below:
   Range of |x|    | Absolute error  | Relative error
-------------------+-----------------+----------------
 [-pi/4; pi/4]     |     1 (Q15)     |
 [pi/4; 7pi/16]    |    15 (Q15)     |    4.6e-4
 [7pi/16; 31pi/64] |   242 (Q15)     |    1.5e-3
-------------------+-----------------+----------------
2 ULP for vtanf(),stanf()
3 ULP for vfasttanf()

Notes for non-fast versions:
1. Fixed-point function result is not defined if input value x is
   +/-pi/2 (+/-8192 in Q15 normalized by pi).
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating-point functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfasttanf():
|x|<804.2477
The output value is not defined outside of this range or accuracy is 
degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vtan) or 8 (vtanf,vfasttanf)
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t stanf ( float32_t x )
{
  xtfloat xin, ztan, zout;
  xb_int32v ji, xin_i;
  xtfloat jf, p, p2, p3, p4, g;
  xtfloat pi2fc0, pi2fc1, pi2fc2;
  xtfloat cf0, cf1, cf2, cf3, cf4, cf5, cf6;
  vbool1 b_cot, b_outl, b_sx;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  xin = x;

  /*
  * Take absolute value of the argument and remember the actual sign. Also
  * check for outliers.
  */

  xin_i = BBE_MOV32_FROMF32(xin);
  p = XT_ABS_S(xin);
  /* NaN or out-of-range? */
  b_outl = XT_ULT_S(tanf_maxval, p);

  if (vbool1_rtor_xtbool(b_outl))
  {
    __Pragma("frequency_hint never");
    if (BBE_MOVAB1(XT_UN_S(xin, xin)))
    {
      errno = EDOM;
    }

    return (qNaNf.f);
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */


  /*
  * Argument reduction.
  */

  /* Determine the pi/2-wide range the input value belongs to. */
  jf = XT_MUL_S(p, inv2pif.f);

  jf = XT_FIROUND_S(jf);
  ji = XT_TRUNC_S(jf, 0);

  /* First 72 bits of pi/2 as a sum of 3 single precision terms. */
  pi2fc0 = pi2fc[0].f;
  pi2fc1 = pi2fc[1].f;
  pi2fc2 = pi2fc[2].f;

  /* Reduce p to [-pi/4,pi/4]. */
  XT_MSUB_S(p, jf, pi2fc0);
  XT_MSUB_S(p, jf, pi2fc1);
  XT_MSUB_S(p, jf, pi2fc2);

  /*
  * Compute the polynomial approximation g(p^2) = tan(p)/p-1. We use a combination
  * of Horner's rule and Estrin's method to evaluate the polynomial.
  */

  cf0 = polytanf_tbl[0].f;
  cf1 = polytanf_tbl[1].f;
  cf2 = polytanf_tbl[2].f;
  cf3 = polytanf_tbl[3].f;
  cf4 = polytanf_tbl[4].f;
  cf5 = polytanf_tbl[5].f;
  cf6 = polytanf_tbl[6].f;

  p2 = XT_MUL_S(p, p);
  p3 = XT_MUL_S(p2, p);
  p4 = XT_MUL_S(p2, p2);

  XT_MADD_S(cf2, cf1, p2); cf1 = cf2;
  XT_MADD_S(cf4, cf3, p2); cf2 = cf4;
  XT_MADD_S(cf6, cf5, p2); cf3 = cf6;

  g = cf0;
  XT_MADD_S(cf1, g, p4); g = cf1;
  XT_MADD_S(cf2, g, p4); g = cf2;
  XT_MADD_S(cf3, g, p4); g = cf3;
  ztan = p;

  /* Free term of the polynomial is zero, thus we obtain the 3rd power. */
  XT_MADD_S(ztan, g, p3);
  zout = ztan;
  /*
  * Conditionally reciprocate the cotan approximation.
  */
  /* Compute the cotangent for odd-numbered segments. */
  {
    xb_vecNx16 h0, h1;
    vboolN s;
    vboolN bcot;
    h0 = BBE_MOVNX16_FROM32(ji);

    h1 = BBE_SLLINX16(h0, 15);
    h0 = BBE_SELNX16I(h1, h1, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    s = BBE_LTNX16(h0, BBE_ZERONX16());
    BBE_EXTRACTB(bcot, bcot, s);
    b_cot = BBE_MOV1_FROMN(bcot);
  }
  if (BBE_MOVAB1(b_cot))
  {
    __Pragma("frequency_hint frequent");
    zout = XT_RECIP_S(ztan);
  }
  /* Compute the sign adjustment term. */
   {
     xb_int32Uv hx;
     hx = XT_CLSFY_S(xin);
     hx = BBE_OPERATOR_AND32U(hx, 0x1);
     b_sx = BBE_OPERATOR_NEQ32U(hx, 0);
   }


  b_sx = BBE_OPERATOR_XORB1(b_sx, b_cot);

  /* Adjust the sign. */
  jf = XT_NEG_S(zout); XT_MOVT_S(zout, jf, b_sx);


  BBE_MOVSCFV(SCF);

  return (zout);

} /* stanf() */
#else
DISCARD_FUN(float32_t,stanf,( float32_t x ))
#endif

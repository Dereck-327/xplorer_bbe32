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
    Modify the Exponent of a Floating-Point Number
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
/*-------------------------------------------------------------------------
Modify the Exponent of a Floating-Point Number

Description: These functions multiply input values by 2^n, where n is an
exponent adjustment term. If the result overflows, functions return 
HUGE_VALF with the proper sign. If, on the contrary, the result underflows,
functions return zero with proper sign.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: exact

Note:
Exponent modification functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-point
exceptions.

Parameters:
Input:
x[N]    Input data
n[N]    Exponent adjustment terms
N       Length of input/output data vectors
Output:
y[N]    Results

Restrictions:
y,x,n   Aligned on 32-byte boundary
y,x,n   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t sldexpf ( float32_t x, int32_t n )
{
  int EXP;
  int e0, e1, e2;
  float32_t s0, s1, s2;
  float32_t y0, y1, y2;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  /* split input 2^n by 3 factors each in range 2^-127...2^127 */
  EXP = n;
  EXP = XT_MAX(-278, EXP);
  EXP = XT_MIN(277, EXP);

  e2 = XT_MAX((-126), XT_MIN(127, EXP));
  e1 = XT_MAX((-126), XT_MIN(127, (EXP - e2)));
  e0 = (EXP - e2) - e1;

  //e0 = (EXP >> 2);
  //e1 = (EXP - e0);
  //e2 = (e1 >> 1);
  //e1 = (e1 - e2);
  s0 = XT_WFR(((uint32_t)(e0 + 127)) << 23);
  s1 = XT_WFR(((uint32_t)(e1 + 127)) << 23);
  s2 = XT_WFR(((uint32_t)(e2 + 127)) << 23);

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  /* multiply the input x by those factors */
  y0 = XT_MUL_S(x, s0);
  y1 = XT_MUL_S(y0, s1);
  y2 = XT_MUL_S(y1, s2);

  /* process errors caused by wrong input data or overflow */
  {
    vbool1 b_xinf, b_yinf, b_edom, b_erange;
    BBE_MOVSCFV(SCF);
    b_xinf = XT_OEQ_S(XT_ABS_S(x), plusInff.f);
    b_yinf = XT_OEQ_S(XT_ABS_S(y2), plusInff.f);
    b_erange = BBE_ANDNOTB1(b_yinf, b_xinf);
    b_edom = XT_UN_S(x, x);
    if ((xtbool)b_edom) errno = EDOM;
    else if ((xtbool)b_erange)
    {
      __Pragma("frequency_hint never");
      errno = ERANGE;
      __feraiseexcept(FE_OVERFLOW);
    }
  }
  return y2;
} /* sldexpf() */
#else
DISCARD_FUN(float32_t,sldexpf,( float32_t x, int32_t n ))
#endif

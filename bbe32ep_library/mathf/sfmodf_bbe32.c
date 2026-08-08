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
    Modulus
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
/*-------------------------------------------------------------------------
Modulus

Description: these functions compute the floating-point remainder that
results from dividing the first argument by the second argument. The result
is less than the second argument and has the same sign as the first argument.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
    x    |   y    | Result | Extra Conditions
  -------+--------+--------+-------------------
   +/-0  |   y    |  +/0   | y!=0
    x    | +/-Inf |   x    | x != +/-Inf
  +/-Inf |   y    |  NaN   | for any y
    x    | +/-0   |  NaN   | for any x
    x    |   y    |  NaN   | |x/y|>=2^24


Accuracy: 2 ULP

Notes:
1. Modulus functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. Modulus functions limit the range of allowable input values, as follows:
   A) If |x/y|>=2^24, then the respective result z is set to NaN
   B) If x==+/-Inf and/or y==+/-0, functions set output value z to NaN, raise
      the "invalid" floating-point exception, and assign EDOM to errno.

Parameters:
Input:
x[N]    Input values
y[N]    Modulus values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t sfmodf ( float32_t x, float32_t y )
{
  static const union ufloat32uint32 _2minus24 = { 0x33800000 };/* 2^-24 */
  xb_int32v sgn;
  float32_t x0, x24, y0, q0, q1, z0, z1, zero, one;
  vbool1 b_qbig, b_yeqz, b_xinf, b_yinf, b_nan;
  vbool1 b_edom, b_fe_inv;
  int err, excepts;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  b_nan = XT_UN_S(x, y);

  excepts = 0;
  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  //BBE_MOVSCFV(XTENSA_SCF_SET_EXCEPTION_ENABLE(SCF, 0)); /* Clear FP flags */

  zero = XT_CONST_S(0);
  one = XT_CONST_S(1);
  sgn = BBE_OPERATOR_AND32(BBE_MOV32_FROMF32(x), MIN_INT32);
  x0 = XT_ABS_S(x);
  y0 = XT_ABS_S(y);
  /* check the ranges of input
  if |x/y|>=2^24, x is finite and y nonzero, return qNaN */
  b_xinf = XT_UEQ_S(x0, plusInff.f);
  b_yinf = XT_OEQ_S(y0, plusInff.f);
  b_yeqz = XT_UEQ_S(y0, zero);

  x24 = XT_MUL_S(x0, _2minus24.f);
  b_qbig = XT_OLE_S(y0, x24);

  q1 = q0 = XT_DIV_S(x0, y0);
  q0 = XT_FIFLOOR_S(q0);
  q1 = XT_SUB_S(q1, one);
  z0 = x0; XT_MSUBN_S(z0, q0, y0);
  z1 = x0; XT_MSUBN_S(z1, q1, y0);
  XT_MOVT_S(z0, z1, XT_OLT_S(z0, zero));

  XT_MOVT_S(z0, x0, BBE_ANDNOTB1(b_yinf, b_xinf));
  z0 = BBE_MOVF32_FROM32(BBE_OPERATOR_OR32(BBE_MOV32_FROMF32(z0), sgn));
  XT_MOVT_S(z0, qNaNf.f, b_qbig);

  /* b_edom = ( x==Inf || x==NaN || y==0 || y==NaN ) */
  b_edom = BBE_OPERATOR_ORB1(b_yeqz, b_xinf);
  /* b_fe_inv = b_edom && x!=NaN && y!=NaN */
  b_fe_inv = BBE_ANDNOTB1(b_edom, b_nan);
  err = errno; XT_MOVT(err, EDOM, b_edom); errno = err;
  XT_MOVT(excepts, FE_INVALID, b_fe_inv);
  BBE_MOVSCFV(SCF);
  __feraiseexcept(excepts);

  return (z0);

} /* sfmodf() */
#else
DISCARD_FUN(float32_t,sfmodf,( float32_t x, float32_t y ))
#endif

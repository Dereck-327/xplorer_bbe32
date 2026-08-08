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
    Reciprocal Square Root
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
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Reciprocal Square Root

Description: Evaluate the reciprocal square root of input value x and store
result to y: y = 1/x^0.5.

Representation:
vrsqrt,srsqrt  32-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-
               defined, provided that it is even.
               Fixed-point format Qy for output data is Qy = 31-Qx/2.
               For example, if Qx == 30 then Qy = 31-30/2 = 16
vrsrtf         IEEE-754 Std. single precision floating-point format

Accuracy:
vrsqrt,srsqrt  1.1e-4 (worst case relative error)
vrsrtf,srsqrtf 2 ULP
vfastrsrtf     3 ULP

Notes (for non-fast versions)::
1. Fixed-point functions return MIN_INT32 (0x80000000) for a negative or
   zero input value.
2. Floating-point reciprocal square root conforms to IEEE-754 Std rSqrt 
   operation in respect to signaling error conditions by means of floating-
   point exceptions.
3. Floating-point reciprocal square root limits the range of allowable
   input values, as follows:
   A) If x<0, functions raise the "invalid" floating-point exception,
      assign EDOM to errno and set output value y to NaN.
   B) If x==+/-0, functions set output value y to +/-HUGE_VALF, raise the
      "divide by zero" floating-point exception, and assign ERANGE to errno.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vrsqrt) or 8 (vrsqrtf,vfastrsqrtf)
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t srsqrtf ( float32_t x )
{
  float32_t z, z_err, spec;
  float32_t zero, one, half, t;
  vbool1 b_zero, b_inf, b_edom, b_nan;

  int err;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  zero = XT_CONST_S(0);
  one = XT_CONST_S(1);
  half = XT_CONST_S(3);

  /* Check input for special cases */
  b_edom = XT_ULT_S(x, zero);
  b_nan = XT_UN_S(x, x);

  /* 1-st rsqrt approximation */
  spec = z = XT_RSQRT0_S(x);

  /* Instruction XT_RSQRT0_S raises all necessary exceptions *
  * for invalid or out-of-range input data.                 */
  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  /* compute approximation error */
  z_err = one;
  t = XT_MUL_S(x, z);
  XT_MSUB_S(z_err, t, z);
  /* 2-nd rsqrt approximation */
  t = XT_MUL_S(half, z);
  XT_MADDN_S(z, t, z_err);
  /* compute approximation error */
  z_err = one;
  t = XT_MUL_S(x, z);
  XT_MSUB_S(z_err, t, z);
  /* 3-rd rsqrt approximaiton */
  t = XT_MUL_S(half, z);
  XT_MADDN_S(z, t, z_err);

  /* Check input for special cases */
  b_inf = XT_OEQ_S(spec, zero);

  b_zero = XT_OEQ_S(x, zero);
  /* Put special values to output: *
  * z=+/-Inf if x=+/-0;           *
  * z=0      if x=Inf.            */
  XT_MOVT_S(z, spec, BBE_OPERATOR_ORB1(b_zero, b_inf));

  BBE_MOVSCFV(SCF); /* Restore floating-point status register to block unnecessary exceptions */
  err = errno;
  /* set errno to EDOM if (x<0 or x==NaN) */
  XT_MOVT(err, EDOM, b_edom);
  /* set errno to ERANGE if (x==0) */
  XT_MOVT(err, ERANGE, b_zero);
  errno = err;

  return z;
} /* srsqrtf() */
#else
DISCARD_FUN(float32_t,srsqrtf,( float32_t x ))
#endif

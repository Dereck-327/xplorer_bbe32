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
    Square Root
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
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Square Root

Description: These functions compute the positive square root of input
vector elements.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
vsqrtf, ssqrtf  2 ULP 
vfastsqrtf      3 ULP

Notes:
1. Square root functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. For a negative input value, functions raise the "invalid" floating-point
   exception, assign EDOM to errno and set the respective output value z to NaN.
3. Negative zero (-0) is considered as a valid input, the result is also -0.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t ssqrtf ( float32_t x )
{
  vbool1 b_edom;
  int err;

  /* set errno to EDOM if (x<0 or x==NaN) */
  b_edom = XT_ULT_S(x, 0.0f);
  err = errno;
  XT_MOVT(err, EDOM, vbool1_rtor_xtbool(b_edom));
  errno = err;
  return XT_SQRT_S(x);
} /* ssqrtf() */
#else
DISCARD_FUN(float32_t,ssqrtf,( float32_t x ))
#endif

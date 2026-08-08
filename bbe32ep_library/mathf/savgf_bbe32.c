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
    Average of Two Values
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Average of Two Values

Description: These functions compute the average of two arguments.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP

Input domain for 'fast' version vfastavgf():
|x+y|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]   Input data
y[N]   Input data
N      Length of input/output data vectors
Output:
z[N]   Results

Restrictions:
z,x,y  Aligned on 32-byte boundary
z,x,y  Must not overlap
N      Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t savgf ( float32_t x, float32_t y )
{
  float32_t half;
  half = (float32_t)XT_CONST_S(3);
  x = XT_MUL_S(half, x);
  XT_MADD_S(x, half, y);
  return x;
} /* savgf() */
#else
DISCARD_FUN(float32_t,savgf,( float32_t x, float32_t y ))
#endif

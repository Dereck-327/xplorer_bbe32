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
    Integer To Floating Value Conversion
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
Integer To Floating Value Conversion

Description: These functions convert integer input values to floating 
values and scale them by 2^t.

Data format: Signed 32-bit integer on input, 
             IEEE-754 Std. single precision floating-point on output.

Parameters:
Input:
x[N]  Input integers
N     Length of input/output data vectors
t     Scale factor
Output:
y[N]  Conversion results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
t     Must belong to [-126,126]
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t sint2float ( int32_t x, int t )
{
  float32_t y, sc;
  uint32_t s;

  ASSERT(t >= -126 && t <= 126);

  /* Convert the input value to floating-point */
  y = ((float32_t)x);
  /* Scale the value */
  s = ((uint32_t)(t + 127)) << 23;
  sc = XT_WFR(s);/* sc=2^t */
  y = y*sc;
  return y;
} /* sint2float() */
#else
DISCARD_FUN(float32_t,sint2float,( int32_t x, int t ))
#endif

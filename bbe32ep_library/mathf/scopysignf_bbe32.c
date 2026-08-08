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
    Copy Sign
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
Copy Sign

Description: These functions copy the sign bit of the second argument to the
first argument and return the result.

Data format: IEEE-754 Std. single precision floating-point.

Note:
In terms of IEEE-754 Std, this is a quiet-computational operation which
treats floating-point numbers and NaNs alike, and does not raise any 
floating-point exceptions.

Parameters:
Input:
x[N]    Input data
y[N]    Sign data
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t scopysignf ( float32_t x , float32_t y)
{
  xb_vecN_2x32v ax, ay;

  /*
  * take the sign of x and y;
  * use fixed-point comparison to distinguish +0/-0
  */
  ax = BBE_MOV32_FROMF32(x);
  ay = BBE_MOV32_FROMF32(y);
  ay = BBE_OPERATOR_AND32(ay, 0x80000000);
  ax = BBE_OPERATOR_AND32(ax, 0x7fffffff);
  ax = BBE_OPERATOR_OR32(ax, ay);
  x = BBE_MOVF32_FROM32(ax);
  return x;
  
} /* scopysignf() */
#else
DISCARD_FUN(float32_t,scopysignf,( float32_t x , float32_t y))
#endif

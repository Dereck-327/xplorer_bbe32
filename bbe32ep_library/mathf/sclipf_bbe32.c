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
    Clipping
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
Clipping

Description: These functions limit the absolute value of inputs. If 
magnitude of input value is less than the second argument then it is left
unchanged. Otherwise it is replaced by absolute value of the second argument
or its negation, depending on the sign of the input value.

Data format: IEEE-754 Std. single precision floating-point.

Parameters:
Input:
x[N]  Input data
y     Limiting value
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t sclipf ( float32_t x, float32_t y )
{
  float32_t xabs, yabs, z;
  xb_vecN_2x32v az, ax;
  vbool1 b_nan;
  b_nan = XT_UN_S(y, y);
  xabs = XT_ABS_S(x);
  yabs = XT_ABS_S(y);
  z = XT_MINNUM_S(xabs, yabs);
  ax = BBE_MOV32_FROMF32(x);
  az = BBE_MOV32_FROMF32(z);
  ax = BBE_OPERATOR_AND32(ax, 0x80000000);
  az = BBE_OPERATOR_OR32(ax, az);
  z = BBE_MOVF32_FROM32(az);
  XT_MOVT_S(z,x,b_nan);
  return (z);

} /* sclipf() */
#else
DISCARD_FUN(float32_t,sclipf,( float32_t x, float32_t y ))
#endif

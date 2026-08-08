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
  NatureDSP_Baseband library. Vector Mathematics
    Reciprocal of Q15 Numbers
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Reciprocal

Description: Evaluate the reciprocal of input value x and store result
to y: y = 1/x.

Representation:
vrecip16,srecip16  16-bit signed fixed-point format
vrecipf            IEEE-754 Std. single precision floating-point format

Fixed-point routines compute reciprocals for Q15 input data, and return the
fractional and exponential parts of the result. Since the reciprocal of 
a 16-bit Q15 is at least 1.0 in magnitude, functions return fractional part
frac in Q(15-exp) format, where exp is the exponential part of the respective
result. Full result can be restored in 32-bit Q15 format by sign extending
the fractional part to 32 bits and shifting it to the left by exp bit positions.
Scalar fixed-point function returns packed 32-bit result, where exponential part
resides in 16 MSBs and fractional part is located in 16 LSBs.

Special cases:
      x    |  Result |  Extra Conditions    
  ---------|---------|---------------------
   +/-Inf  | +/-0    | vrecipf,srecipf
   +/-0    | +/-Inf  |  
  ---------|---------|---------------------
     0     |   not   | vrecip16,srecip16
           | defined |

Input domain for vfastrecipf:
|x|>2.94e-39, |x|<Inf
The output value is not defined outside of this range.

Accuracy:
vrecip16,srecip16  1 LSB of the fractional part
vrecipf,srecipf    1 ULP 
vfastrecipf        2 ULP 

Parameters:
Input:
x[N]           Input data vector
N              Length of vectors
Output:
vrecip16,srecip16
frac[N]        Fractional part of reciprocals, Q(15-exp)
exp[N]         Exponent of reciprocals (1...16)
vrecip:
y[N]           Reciprocals

Restrictions:
y,x,fract,exp  Aligned on 32-byte boundary
y,x,fract,exp  Must not overlap
N              Multiple of 16 (vrecip16), 8 (vrecipf,vfastrecipf)
-------------------------------------------------------------------------*/

#if !(HAVE_RECIP && 1)

DISCARD_FUN(uint32_t, srecip16, (int16_t x))

#else

uint32_t srecip16(int16_t x)
{
    xb_vecNx16 z0, z1, z2, x0, y0, e0, b;
    xb_vecNx40 xw, a;
    vsaN nsa, shft;
    vboolN    b0;
    uint16_t   y, e;

    xb_vecNx16 SpecX, _1;
    vboolN SpecXMask;

    // init correction variables to process special X value (see "if (xi == 16384)" line in generic code - file recip16.c)
    _1 = BBE_MOVVINT16(1);
    //SpecX = BBE_MOVVINT16(0x4000);
    SpecX = 0x4000;

    z0 = 0;
    z1 = BBE_MOVVINT16(1);
    z2 = BBE_MOVVINX16(BBE_MOVVI_INT16_MAXINT);
    shft = BBE_MOVVSA32(23);

    x0 = BBE_MOVVA16(x);
    b0 = BBE_EQNX16(x0, z0);
    nsa = BBE_NSANX16(x0);
    y0 = BBE_SLLNX16(x0, nsa);
    xw = BBE_UNPKSNX16(y0);
    xw = BBE_SLLINX40(xw, 24);

    SpecXMask = BBE_EQNX16(y0, SpecX); // calculate mask of special values

    BBE_RECIPLUNX40_0(a, y0, b, xw);
    BBE_MULUSANX16(a, b, y0);
    y0 = BBE_PACKVNX40(a, shft);
    y0 = BBE_MOVNX16T(z2, y0, b0);
    y0 = BBE_MOVNX16T(SpecX, y0, SpecXMask); // correction of fractional part

    e0 = BBE_MOVVVS(nsa);
    e0 = BBE_ADDNX16(e0, z1);
    BBE_ADDNX16T(e0, e0, _1, SpecXMask); // correction of exponential part

    y = BBE_MOVAV16(y0);
    e = BBE_MOVAV16(e0);
    return (uint32_t)((e << 16) | y);
} /* srecip16() */

#endif

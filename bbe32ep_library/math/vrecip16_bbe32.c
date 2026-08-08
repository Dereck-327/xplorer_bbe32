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

#if !(HAVE_VSAMATH && HAVE_RECIP && 1)

DISCARD_FUN(void, vrecip16, (int16_t * restrict  frac,
                             int16_t * restrict  exp,
                             const int16_t * x,
                             int N))

#else

void vrecip16 ( int16_t * restrict fract, 
                int16_t * restrict exp, 
          const int16_t * restrict x, 
          int N)
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pF = (xb_vecNx16 *)fract;
    xb_vecNx16 * restrict pE = (xb_vecNx16 *)exp;

    int k;
    xb_vecNx16 z0, z2, x0, y0, e0, b, _1;
    xb_vecNx40 xw, a;
    vsaN shft, nsa24;
    vsaN vsa24;
    vboolN    b0;
    xb_vecNx16 SpecY;
    vboolN SpecXMask;

    NASSERT_ALIGN(fract, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(exp, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N%BBE_SIMD_WIDTH == 0);

    z0 = 0;
    z2 = BBE_MOVVINX16(BBE_MOVVI_INT16_MAXINT);
    shft = BBE_MOVVSA32(23);
    vsa24 = BBE_MOVVSA32(24);

    // init correction variables to process special x[i] value (see "if (xi == 16384)" line in generic code)
    _1 = BBE_MOVVINT16(1);
    SpecY = 0x4000;

    for (k = 0; k<(N >> LOG2_BBE_SIMD_WIDTH); ++k)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
        b0 = BBE_EQNX16(x0, z0);

        nsa24 = BBE_NSANX16(x0);
        xw = BBE_UNPKSNX16(x0);
        nsa24 = BBE_ADDSAVSN(24, nsa24);
        xw = BBE_SLSNX40(xw, nsa24);

        //// calculate mask of special values
        x0 = BBE_PACKVNX40(xw, vsa24);
        SpecXMask = BBE_EQNX16(x0, SpecY);

        BBE_RECIPLUNX40_0(a, y0, b, xw);
        BBE_RECIPLUNX40_1(a, y0, b, xw);

        BBE_MULUSANX16(a, b, y0);
        y0 = BBE_PACKVNX40(a, shft);
        y0 = BBE_MOVNX16T(z2, y0, b0);
        y0 = BBE_MOVNX16T(SpecY, y0, SpecXMask); // correction of fractional part

        nsa24 = BBE_ADDSAVSN(-23, nsa24);
        e0 = BBE_MOVVVS(nsa24);
        BBE_ADDNX16T(e0, e0, _1, SpecXMask); // correction of exponential part

        BBE_SVNX16_IP(y0, pF, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(e0, pE, 2 * BBE_SIMD_WIDTH);
    }
} /* vrecip16() */

#endif

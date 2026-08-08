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
    Division of Q15 Numbers
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
Division

Representation:
vdivide,sdivide    16-bit signed fixed-point format (Q15)
vdividef,sdividef  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of the 
division result. Fixed-point format for the fractional part is 16-bit 
Q(15-exp), where exp denotes the respective exponential value. Full division 
result can be restored in 48-bit Q31 format by sign extending the fractional 
part to 64 bits and shifting it to the left by 16+exp bit positions
Scalar fixed-point function returns packed 32-bit result, where exponential part
resides in 16 MSBs and fractional part is located in 16 LSBs.

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number  (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number 
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                       (fixed-point functions)
          |         | defined |

Accuracy:
vdivide,sdivide   1 LSB of the fractional part
vdividef,dividef  1 ULP

Parameters:
Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdivide,sdivide
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[N]    Exponential part of quotients, -14..16
vdividef,sdividef
z[N]      Quotients

Restrictions:
z,x,y,fract,exp   Aligned on 32-byte boundary
z,x,y,fract,exp   Must not overlap
N                 Multiple of 16 (vdivide) or 8 (vdividef)
-------------------------------------------------------------------------*/

#if !HAVE_DIV

DISCARD_FUN(uint32_t, sdivide, (int16_t x, int16_t y))

#else

uint32_t sdivide (int16_t x,int16_t y)
{
    uint16_t   r, e;
    xb_vecNx16 z0, z1, x0, y0, e0, e1, r0;
    xb_vecNx40 xw;
    vsaN nsax, nsay;
    z0 = BBE_MOVPINT16(16);
    z1 = BBE_MOVVINT16(1);

    x0 = BBE_MOVVA16(x);
    y0 = BBE_MOVVA16(y);
    nsax = BBE_NSANX16(x0);
    nsay = BBE_NSANX16(y0);
    x0 = BBE_SLANX16(x0, nsax);
    y0 = BBE_SLANX16(y0, nsay);
    xw = BBE_MULNX16(x0, z0);
    r0 = BBE_QUONX32(xw, y0);
    e0 = BBE_MOVVVS(nsax);
    e1 = BBE_MOVVVS(nsay);
    e0 = BBE_SUBNX16(e1, e0);
    e0 = BBE_ADDNX16(e0, z1);
    r = BBE_MOVAV16(r0);
    e = BBE_MOVAV16(e0);
    return (uint32_t)((e << 16) | r);
} /* sdivide() */

#endif

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
#include "NatureDSP_Baseband_complex.h"

/*-------------------------------------------------------------------------
Division of Complex Numbers

Representation:
vdividecx         16-bit signed fixed-point format (Q15)
vdividecxf,cdivf  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of real
and imaginary components of the division result. Fixed-point format for the
fractional part is 16-bit Q(15-exp), where exp denotes the respective 
exponential value. Full division result can be restored in 48-bit Q31 format
by sign extending the fractional part to 64 bits and shifting it to the left
by 16+exp bit positions

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number          (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number (see *)
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                               (fixed-point functions)
          |         | defined |

* - quotent of infinite number to finite number may have NaN in its real or 
    imaginary part depending on signs of input arguments

Input domain for 'fast' version vfastdividecxf is limited by usage of direct 
formula without normalization of inputs, so input domain is:
|real(x)*real(y)|<Inf
|imag(x)*imag(y)|<Inf
|real(y)*real(y)|<Inf
|imag(y)*imag(y)|<Inf
1.1755e-038 < |real(x)*real(y) + imag(x)*imag(y)| < Inf
1.1755e-038 < |real(y)*real(y) + imag(y)*imag(y)| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
vdividecx         1.5 LSB of the fractional part
vdividecxf,cdivf  2 ULP
vfastdividecxf    2 ULP

Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdividecx:
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[2*N]  Exponential part of quotients, -14..16. Exponential values for
          real and imaginary components are interleaved, with exponential
          parts of real components going first (at even indices)
vdividecxf,cdivf
z[N]      Quotients
Restrictions:
x,y,z,fract,exp   Aligned on 32-byte boundary
x,y,z,fract,exp   Must not overlap
N                 Multiple of 8 (vdividecx) or 4 (vdividecxf, vfastdividecxf)
---------------------------------------------------------------------------*/

#if !HAVE_DIV

DISCARD_FUN(void, vdividecx, (complex_fract16 * restrict fract, 
                 int16_t         * restrict exp, 
           const complex_fract16 * restrict x, 
           const complex_fract16 * restrict y,
           int N))

#else

void vdividecx(complex_fract16 * restrict fract, 
                 int16_t         * restrict exp, 
           const complex_fract16 * restrict x, 
           const complex_fract16 * restrict y,
           int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pY = (const xb_vecNx16 *)y;
    xb_vecNx16 * restrict pF = (xb_vecNx16 *)fract;
    xb_vecNx16 * restrict pE = (xb_vecNx16 *)exp;
    xb_vecNx16 x0, y0, r0, e0, vvww, den16, invd;
    xb_vecNx40 den, vw, _4000L, invd40;
    vsaN       evw, ed, RndShift, NoShift, shift;

    int n;

    NASSERT_ALIGN(fract, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(exp, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

    _4000L = BBE_MOVWA40(0, 0x40000000);

    RndShift = BBE_MOVVSA32(-1); // result shift will be (15 - (-1))
    NoShift = BBE_MOVVSA32(15); // result shift will be (15 - 15)
    shift = BBE_MOVVSA32(24);

    for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(y0, pY, 2 * BBE_SIMD_WIDTH);

        den = BBE_MAGINX16C(y0, y0); // |y0|^2

        vw = BBE_MULNX16J(x0, y0); // complex multiplication of x0 to conjugate of y0

        // normalize nominator and denominator 
        evw = BBE_NSANX40(vw);
        vw = BBE_SLANX40(vw, evw);
        vvww = BBE_PACKHNX40(vw);
        // calc 1/den with saturetion
        ed = BBE_NSANX40(den);
        ed = BBE_SUBSVSN(ed, RndShift);    // correction "nsa" to "count of leading zeros"
        den = BBE_SLANX40(den, ed);
        den16 = BBE_PACKVNX40(den, shift);
        invd = BBE_QUONX32U(_4000L, den16); // unsigned division
        // saturation
        invd40 = BBE_UNPKUNX16(invd);
        invd = BBE_PACKNVNX40(invd40, NoShift);

        // multiplication complex value vw to (1/den)
        vw = BBE_MULNX16(vvww, invd);
        vw = BBE_SLLINX40(vw, 1); // multipte result to 2 

        r0 = BBE_PACKNVNX40(vw, RndShift); // rounding and saturation

        e0 = BBE_SUBSVSN(ed, evw);

        BBE_SVNX16_IP(r0, pF, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(e0, pE, 2 * BBE_SIMD_WIDTH);
    }
} /* vdividecx() */

#endif

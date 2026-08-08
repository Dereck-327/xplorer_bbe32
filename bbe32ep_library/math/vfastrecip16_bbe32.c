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
    Reciprocal for Pseudo-Floating Point Format
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
Reciprocal for Pseudo-Floating Point Format

Description: calculates y=1/x; input and output values are represented
by 7-bit signed exponent and 16-bit signed mantissa, Q(15+exp).
If the magnitude of result is so large that the result cannot be
represented in the pseudo-floating point format, then the function
returns 32767/-32768,-64 (+/-1.8e19) depending on the sign bit of the
input mantissa.

Relative accuracy: 7.7e-3 (~7 bits of precision) in worst case

Parameters:
Input:
xmant[N]                Mantissa of input values, -32768..32767
xexp[N]                 Exponent of input values, -64..63
Output:
ymant[N]                Mantissa of output values, -32768..32767
yexp[N]                 Exponent of output values, -64..63

Restrictions:
xmant,ymant,xexp,yexp   Aligned on 32-byte boundary
xmant,ymant,xexp,yexp   Must not overlap
N                       Multiple of 16
-------------------------------------------------------------------------*/

void vfastrecip16 ( int16_t * restrict ymant,
                    int16_t * restrict yexp,
              const int16_t * restrict xmant,
              const int16_t * restrict xexp,
              int N )
{
    xb_vecNx16 * restrict YM;
    xb_vecNx16 * restrict YE;
    const xb_vecNx16 * XM;
    const xb_vecNx16 * XE;

    xb_vecNx16 xm, xe;
    xb_vecNx16 ym, ye;
    vsaN   vsa_xn, vsa_ye;

    int n;

    NASSERT_ALIGN(xmant, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(xexp, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(ymant, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(yexp, (2 * BBE_SIMD_WIDTH));
    NASSERT((N % BBE_SIMD_WIDTH) == 0);

    //
    // Process data.
    //

    YM = (xb_vecNx16*)ymant;
    YE = (xb_vecNx16*)yexp;
    XM = (const xb_vecNx16*)xmant;
    XE = (const xb_vecNx16*)xexp;

#ifdef COMPILER_XTENSA
#pragma ymemory( XM )
#pragma ymemory( XE )
#endif

    for (n = 0; n<N / BBE_SIMD_WIDTH; n++)
    {
        // Load input data, Q(15+xe)
        BBE_LVNX16_IP(xm, XM, +2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(xe, XE, +2 * BBE_SIMD_WIDTH);

        vsa_xn = BBE_NSANX16(xm);
        // Mantissa normalization; Q(15+xn+xe) <- Q(15+xe) + xn
        xm = BBE_SLANX16(xm, vsa_xn);

#if (XCHAL_HAVE_BBEN_ADVPRECISION==1)
        {
            vsaN   vsa_xe;
            vsa_xe = BBE_MOVVSV(xe, 0);
            vsa_xn = BBE_ADDSVSN(vsa_xe, vsa_xn);
        }
#else
        {
            xb_vecNx16 xe1;
            xe1 = BBE_ADDNX16(xe, BBE_MOVVVS(vsa_xn));
            xe1 = BBE_MAXNX16(xe1, -64);
            xe1 = BBE_MINNX16(xe1, 63);
            vsa_xn = BBE_MOVVSV(xe1, 0);
        }
#endif

        // Take reciprocals; Q(15-1-xn-xe) <- Q29/Q(15+xn+xe)
        BBE_FPRECIPNX16_0(ym, vsa_ye, xm, vsa_xn);
        BBE_FPRECIPNX16_1(ym, vsa_ye, xm, vsa_xn);

        // ye <=  -1 - xn
        ye = BBE_MOVVVS(vsa_ye);

        // Save results, Q(15+ye)
        BBE_SVNX16_IP(ym, YM, +2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(ye, YE, +2 * BBE_SIMD_WIDTH);
    }
} /* vfastrecip16() */

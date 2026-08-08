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
    Inverse Square Root for Pseudo-Floating Point Format
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
Inverse Square Root for Pseudo-Floating Point Format

Description: calculates y=1/sqrt(x); input and output values are represented
by 16-bit signed mantissa and 7-bit signed exponent, Q(15+exp). Returns 
-32768,0 (-1.0) for input values with negative mantissa. For zero mantissa
the result is 32767,-64 (~1.8e19).

Relative accuracy: 2.0e-3 (~9 bits of precision) in worst case

Parameters:
Input:
xmant[N]                Mantissa of input values, -32768..32767
xexp[N]                 Exponent of input values, -64..62
Output:
ymant[N]                Mantissa of output values, -32768,16384..32767
yexp[N]                 Exponent of output values, -64,-40..31

Restrictions:
xmant,ymant,xexp,yexp   Aligned on 32-byte boundary
xmant,ymant,xexp,yexp   Must not overlap
N                       Multiple of 16
-------------------------------------------------------------------------*/

void vfastrsqrt ( int16_t * restrict ymant,
                  int16_t * restrict yexp,
            const int16_t * restrict xmant,
            const int16_t * restrict xexp,
            int N )
{
    xb_vecNx16 * restrict YM;
    xb_vecNx16 * restrict YE;
    const xb_vecNx16 *          XM;
    const xb_vecNx16 *          XE;

    xb_vecNx16 xm, xe, eadd;
    xb_vecNx16 ym, ye;
    xb_vecNx16 zm;
    xb_vecNx16 c0;
    xb_vecNx16 _inv1, _1;

    vsaN   vsa_xn, vsa_xe, vsa_ye;
    vboolN vb0, vb1;

    int n;

    NASSERT_ALIGN(xmant, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(xexp, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(ymant, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(yexp, (2 * BBE_SIMD_WIDTH));

    NASSERT((N % BBE_SIMD_WIDTH) == 0);

    //
    // Setup special values to handle non-positive mantissa on input.
    //

    c0 = 0;
    zm = BBE_MOVVINX16(BBE_MOVVI_INT16_MININT);
    //OLD: _inv1 = BBE_MOVVINT16(0xfffe);
    _inv1 = 0xfffe;
    _1 = BBE_MOVVINT16(1);

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
        // Check if input mantissa is zero
        vb1 = BBE_EQNX16(xm, c0);
        // Force zero exponent for zero mantissa values.
        BBE_LVNX16F_IP(xe, XE, +2 * BBE_SIMD_WIDTH, vb1);

        vsa_xn = BBE_NSANX16(xm);
        // Mantissa normalization; Q(15+xn+xe) <- Q(15+xe) + xn
        xm = BBE_SLANX16(xm, vsa_xn);

        // NOTE: If exp<=-62 or exp>=62 sometimes result of BBE_FPRSQRTNX16_0()/BBE_FPRSQRTNX16_1() is incorrect (overflow of VSA register).
        // =>    process exponent itself: exp == EvenExp + ZeroOrOne, calculate 1/sqrt(x) with small exponent ZeroOrOne and make correction (ResExp - EvenExp/2)
        eadd = BBE_ANDNX16(xe, _inv1);
        xe = BBE_ANDNX16(xe, _1);

#if XCHAL_HAVE_BBEN_ADVPRECISION
        vsa_xe = BBE_MOVVSV(xe, 0);
        vsa_xe = BBE_ADDSVSN(vsa_xe, vsa_xn); // Calculate full exponent after mantissa normalization (exclude BIG EVEN VALUE eadd)
#else
        {
            xb_vecNx16 xe1;
            xe1 = BBE_ADDNX16(xe, BBE_MOVVVS(vsa_xn));
            xe1 = BBE_MAXNX16(xe1, -64);
            xe1 = BBE_MINNX16(xe1, 63);
            vsa_xe = BBE_MOVVSV(xe1, 0);
        }
#endif
        // Take reciprocal square roots; Q(15-1-((xe+xn)&~1)/2)
        BBE_FPRSQRTNX16_0(ym, vsa_ye, xm, vsa_xe);
        BBE_FPRSQRTNX16_1(ym, vsa_ye, xm, vsa_xe);

#ifdef BBE_SUBSRA1SVSN
        vsa_ye = BBE_SUBSRA1SVSN(vsa_ye, BBE_MOVVSV(eadd, 0));
        ye = BBE_MOVVVS(vsa_ye);
#else
        ye = vsa_ye;
        ye = BBE_SUBNX16(ye, BBE_SRAINX16(eadd, 1));
        vsa_ye = ye;
#endif

        vb0 = BBE_LTNX16(xm, c0); // Check if input mantissa is negative

        // Force a special negative number on output for non-positive mantissa on input.
        ym = BBE_MOVNX16T(zm, ym, vb0);
        ye = BBE_MOVNX16T(c0, ye, vb0);

        // Save results, Q(15+ye)
        BBE_SVNX16_IP(ym, YM, +2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(ye, YE, +2 * BBE_SIMD_WIDTH);
    }
} /* vfastrsqrt() */

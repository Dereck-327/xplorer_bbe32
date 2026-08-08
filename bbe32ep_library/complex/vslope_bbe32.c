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
    Correction of Vector Slope
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"


/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_complex.h"
#include "cbexpTbl.h"

/*-------------------------------------------------------------------------
Correction of Vector Slope

Description: Function adjusts the phase of complex values from input vector
by a sequence of linearly evolving angle values.

Representation: 16-bit signed fixed-point

Accuracy: absolute error for re/im component does not exceed 4.0e-4 (13 in Q15).

Note:
Phase adjustment of k-th complex input value is performed through the 
multiplication by a conjugated rotating factor lying on the unit circle. This
factor is calculated from the full-precision phase angle: phase0+k*delta, thus
the function is able to process long data sequences with no risk of phase 
errors accumulation.

Parameters:
Input:
phase0   Initial phase divided by pi, Q15
delta    Phase increment divided by pi, Q15
N        Size of input/output array
Input/Output:
x[N]     Input/output complex values, Q15; real and imaginary parts are 
         Interleaved with the real part stored at even indices
Returned Value:
         Updated phase, phase0+N*delta

Restrictions:
x        Aligned on 32-byte boundary
N        Multiple of 16
-------------------------------------------------------------------------*/

int16_t vslope (complex_fract16 * restrict x, int16_t phase0, int16_t delta, int N)
{
    int n;
    const xb_vecNx16 * restrict pT = (const xb_vecNx16 *)cbexpTbl;
    xb_vecNx16 * restrict pX = (xb_vecNx16 *)x;
    xb_vecNx16 tsin0l, tsin0h, tsin1, tcos1;
    xb_vecNx16 z, yc0, ys0, yc1, ys1, yc2, ys2, addx;
    xb_vecNx16 _7, _4000, _6434;
    xb_vecNx16 v0, v1, u0, u1;
    vselN ix0, ix01, ix1;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH) == 0);

    _7 = BBE_MOVVA16(7);
    _4000 = BBE_MOVVA16(0x4000);
    yc2 = BBE_MOVVA16(32763);
    _6434 = BBE_MOVVA16(6434);

    // init argument for first iteration
    z = BBE_SEQNX16();      // sequence: 0,1,2, ... 7
    v0 = BBE_MOVVA16(delta);
    z = BBE_MULNX16PACKL(z, v0); // sequence: d, 2*d, ... 7*d
    v0 = BBE_MOVVA16(phase0);
    z = BBE_ADDNX16(z, v0); // sequence: a0+d, a0+2*d, ... a0+7*d

    addx = BBE_MOVVA16(delta << LOG2_BBE_SIMD_WIDTH);

    BBE_LVNX16_XP(tsin0l, pT, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(tsin0h, pT, 2 * BBE_SIMD_WIDTH);
    tcos1 = BBE_SELNX16I(tsin0h, tsin0h, BBE_SELI_EXTRACT_HI_HALVES);
    tsin1 = BBE_SELNX16I(tsin0h, tsin0h, BBE_SELI_EXTRACT_LO_HALVES);

    for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        // first approximation 
        ix0 = BBE_MOVVSELNX16(z, 12);

        v0 = BBE_SRAINX16(z, 9);
        v0 = BBE_ANDNX16(v0, _7);
        ix1 = BBE_MOVVSELNX16(v0, 0);

        v0 = BBE_ADDNX16(z, _4000);
        ix01 = BBE_MOVVSELNX16(v0, 12);

        yc0 = BBE_SHFLNX16(tsin0l, ix01);
        ys0 = BBE_SHFLNX16(tsin0l, ix0);

        // second approximation 
        
        yc1 = BBE_SHFLNX16(tcos1, ix1);
        ys1 = BBE_SHFLNX16(tsin1, ix1);

        v0 = BBE_SLLINX16(z, 4);
        v0 = BBE_POLYNX16_OFF(v0, 13, 0); // offset in interval

        // last appoximation: ~ x*pi 
        // "yc2": see initialization of loop
        ys2 = BBE_MULNX16PACKQ(v0, _6434);

        // 2 complex multiplies 
        u0 = BBE_SELNX16I(ys0, yc0, BBE_SELI_INTERLEAVE_1_LO);
        u1 = BBE_SELNX16I(ys0, yc0, BBE_SELI_INTERLEAVE_1_HI);
        v0 = BBE_SELNX16I(ys1, yc1, BBE_SELI_INTERLEAVE_1_LO);
        v1 = BBE_SELNX16I(ys1, yc1, BBE_SELI_INTERLEAVE_1_HI);
        yc1 = BBE_MULNX16CPACKQ(u0, v0); // complex multiplcation of low parts
        ys1 = BBE_MULNX16CPACKQ(u1, v1); // complex multiplcation of hi parts

        u0 = BBE_SELNX16I(ys2, yc2, BBE_SELI_INTERLEAVE_1_LO);
        u1 = BBE_SELNX16I(ys2, yc2, BBE_SELI_INTERLEAVE_1_HI);
        v0 = BBE_MULNX16CPACKQ(u0, yc1); // complex multiplcation of low parts
        v1 = BBE_MULNX16CPACKQ(u1, ys1); // complex multiplcation of hi parts

        // load block of complex numbers from array
        u0 = BBE_LVNX16_I(pX, 0 * (2 * BBE_SIMD_WIDTH));
        u1 = BBE_LVNX16_I(pX, 1 * (2 * BBE_SIMD_WIDTH));

        // third multiplication (u0, u1 to conjugate of v0, v1)
        yc1 = BBE_MULNX16JPACKQ(u0, v0); // complex multiplcation of low parts
        ys1 = BBE_MULNX16JPACKQ(u1, v1); // complex multiplcation of hi parts

        // save complex coordinates
        BBE_SVNX16_IP(yc1, pX, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(ys1, pX, 2 * BBE_SIMD_WIDTH);

        // update argument
        z = BBE_ADDNX16(z, addx);
    }

    return(BBE_MOVAV16(z));
} /* vslope() */

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
    Complex Vector Exponent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"


/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_complex.h"
#include "sineTbl.h"

/*-------------------------------------------------------------------------
Complex Exponential with Linearly Evolving Phase

Description: These functions compute a series of complex numbers lying on
the unit circle, with the phase angle being incremented for each successive
complex number.

Representation:
cbexplin,           16-bit signed fixed-point format
cbexplin_fast       Initial phase and phase increment (phase0 and delta) are 
                    Q15 angular values normalized by pi.
                    Output data z[] format is Q15.
cbexplinf,          IEEE-754 Std. single precision floating-point format
cexplinf            Initial phase and phase increment (phase0 and delta) are
                    in radians. See Note 3.

Accuracy:
cbexplin            Absolute error for re/im components does not exceed
                    3.7e-4 (12 in Q15)
cbexplin_fast       Absolute error for re/im components for n-th result is
                    ~(floor(n/8)+1)*3.7e-4, n=0..N-1
cbexplinf,cexplinf  2 ULP 
cbfastexplinf       3 ULP

Notes:
1. cbexplin() function computes each complex value from a full-precision
   linearly evolving angle, thus the function is able to generate long data
   sequences with no risk of phase errors accumulation.
2. cbexplin_fast() function is approximately twice faster but less accurate 
   because it performs successive rotations by small angle on each iteration.
   Avoid using this function for long sequences or tiny phase increments.
3. For the non-fast floating-point functions (cbexplinf, cexplinf) , the 
   phase (initial and updated) should belong to the range [-102940.0, 
   102940.0], otherwise the respective resultis 0+0j.
4. Scalar function cexplinf uses the initial phase angle (phase0) to compute 
   a single complex number, and returns the once incremented phase value. 

Input domain for 'fast' version (cbfastexplinf)
|phase|<804.2477
for all initial and updated phases 
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase0  Initial phase angle (phase of the first complex number in a series)
delta   Phase increment step
N       Length of output vector
Output:
z[N]    Complex numbers
Returned Value:
        Updated phase, phase0+N*delta

Restrictions:
z       Aligned on 32-byte boundary
z       Must not overlap
N       Multiple of 16 (cbexplln, cbexplln_fast) or 8 (cbexplinf,cbfastexplinf)
-------------------------------------------------------------------------*/
int16_t cbexplin_fast(complex_fract16 * restrict z, int16_t   phase0, int16_t   delta, int N)
{
    int n;
    xb_vecNx16 ph, t, addph;
    xb_vecNx40 tw;
    const xb_vecNx16 * restrict T = (const xb_vecNx16 *)sineTbl;
    xb_vecNx16 * restrict pZ = (xb_vecNx16*)z;
    xb_vecNx16 t1, t3, t4;
    xb_vecNx16 c, s, cs0, re, im;
    xb_vecNx16 c4000;

    NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH) == 0);
    if (N <= 0) return (phase0);

    t1 = BBE_LVNX16_I(T, 2 * BBE_SIMD_WIDTH);
    t3 = BBE_LVNX16_I(T, 6 * BBE_SIMD_WIDTH);
    t4 = BBE_LVNX16_I(T, 8 * BBE_SIMD_WIDTH);

    // init argument for first iteration
    ph = BBE_SEQNX16();      // sequence: 0,1,2, ... 7
    t = BBE_MOVVA16(delta);
    tw = BBE_MULNX16(ph, t); // sequence: d, 2*d, ... 7*d
    ph = BBE_PACKLNX40(tw);
    t = BBE_MOVVA16(phase0);
    ph = BBE_ADDNX16(ph, t); // sequence: a0+d, a0+2*d, ... a0+7*d

    addph = BBE_MOVVA16(delta << LOG2_BBE_SIMD_WIDTH);

    // calculate sin(ph) for first block
    {
        vselN sel;
        xb_vecNx16 ofs, tmp_t, x0, y0, SrcX;

        SrcX = ph; //store sign

        x0 = BBE_ABSNX16(ph);
        x0 = BBE_SLLINX16(x0, 1);

        // POLI_SU
        ofs = BBE_POLYNX16_OFF(x0, 12, 1); // offset in interval
        sel = BBE_MOVVSELNX16(x0, 12);     // interval number (index in sine table)

        tmp_t = BBE_SHFLNX16(t1, sel);
        t = BBE_MULNX16PACKQ(tmp_t, ofs);

        tmp_t = BBE_SHFLNX16(t3, sel);
        y0 = BBE_ADDNX16(tmp_t, t);

        t = BBE_MULNX16PACKQ(y0, ofs);

        tmp_t = BBE_SHFLNX16(t4, sel);
        y0 = BBE_ADDSNX16(tmp_t, t);

        s = BBE_MULSGNNX16(SrcX, y0); // correct sign
    }

    // calculate cos(ph) for first block
    {
        vselN sel;
        xb_vecNx16 ofs, tmp_t, x0, y0, SrcX;

        c4000 = BBE_MOVVA16(0x4000);
        ph = BBE_ADDNX16(ph, c4000); // cos(x) == sin(x+16384)

        SrcX = ph; //store sign

        x0 = BBE_ABSNX16(ph);
        x0 = BBE_SLLINX16(x0, 1);

        // POLI_SU
        ofs = BBE_POLYNX16_OFF(x0, 12, 1); // offset in interval
        sel = BBE_MOVVSELNX16(x0, 12);     // interval number (index in sine table)

        tmp_t = BBE_SHFLNX16(t1, sel);
        t = BBE_MULNX16PACKQ(tmp_t, ofs);

        tmp_t = BBE_SHFLNX16(t3, sel);
        y0 = BBE_ADDNX16(tmp_t, t);

        t = BBE_MULNX16PACKQ(y0, ofs);

        tmp_t = BBE_SHFLNX16(t4, sel);
        y0 = BBE_ADDSNX16(tmp_t, t);

        c = BBE_MULSGNNX16(SrcX, y0); // correct sign
    }

    // calculate sin(Delta), cos(Delta)
    {
        vselN sel;
        xb_vecNx16 ofs, tmp_t, x0, y0, SrcX, addph2;

        addph2 = BBE_ADDNX16(addph, c4000);

        // build interleaved vector of pairs (Delta+16384, Delta) for calculate pairs (cos(Delta), sin(Delta))
        ph = BBE_SELNX16I(addph, addph2, BBE_SELI_INTERLEAVE_1_LO);

        SrcX = ph; //store sign

        x0 = BBE_ABSNX16(ph);
        x0 = BBE_SLLINX16(x0, 1);

        // POLI_SU
        ofs = BBE_POLYNX16_OFF(x0, 12, 1); // offset in interval
        sel = BBE_MOVVSELNX16(x0, 12);     // interval number (index in sine table)

        tmp_t = BBE_SHFLNX16(t1, sel);
        t = BBE_MULNX16PACKQ(tmp_t, ofs);

        tmp_t = BBE_SHFLNX16(t3, sel);
        y0 = BBE_ADDNX16(tmp_t, t);

        t = BBE_MULNX16PACKQ(y0, ofs);

        tmp_t = BBE_SHFLNX16(t4, sel);
        y0 = BBE_ADDSNX16(tmp_t, t);

        cs0 = BBE_MULSGNNX16(SrcX, y0); // correct sign
    }

    // save first block
    re = BBE_SELNX16I(s, c, BBE_SELI_INTERLEAVE_1_LO);
    im = BBE_SELNX16I(s, c, BBE_SELI_INTERLEAVE_1_HI);
    

    // process other blocks
    for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH ); n++)
    {
        BBE_SVNX16_IP(re, pZ, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(im, pZ, 2 * BBE_SIMD_WIDTH);
        re = BBE_MULNX16CPACKQ(re, cs0);
        im = BBE_MULNX16CPACKQ(im, cs0);
        //BBE_SVNX16_IP(re, pZ, +2 * BBE_SIMD_WIDTH);
        //BBE_SVNX16_IP(im, pZ, +2 * BBE_SIMD_WIDTH);
    }

    return phase0 + N*delta;
} /* cbexplin_fast() */

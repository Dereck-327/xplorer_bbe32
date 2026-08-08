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
#include "cbexpTbl.h"

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
int16_t cbexplin(complex_fract16 * restrict z, int16_t   phase0, int16_t   delta, int N)
{
  int n;
  const xb_vecNx16 * restrict pT = (const xb_vecNx16 *)cbexpTbl;
  xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
  xb_vecNx16 tsin0l, tsin0h, tsin1, tcos1;
  xb_vecNx16 x, yc0, ys0, yc1, ys1, yc2, ys2, addx;
  xb_vecNx16 _4000, _6434;
  xb_vecNx16 v0, v1, u0, u1;
  xb_vecNx40 w0;
  vselN ix0, ix01, ix1;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH) == 0);

  //tsin0l = BBE_LVNX16_I(pT, 0 * (2 * BBE_SIMD_WIDTH));
  //tsin0h = BBE_LVNX16_I(pT, 1 * (2 * BBE_SIMD_WIDTH));
  //tsin1  = BBE_LVNX16_I(pT, 2 * (2 * BBE_SIMD_WIDTH));
  //tcos1  = BBE_LVNX16_I(pT, 3 * (2 * BBE_SIMD_WIDTH));

  _4000 = BBE_MOVVA16(0x4000);
  yc2 = BBE_MOVVA16(32763);
  _6434 = BBE_MOVVA16(6434);

  // init argument for first iteration
  x = BBE_SEQNX16();      // sequence: 0,1,2, ... 7
  v0 = BBE_MOVVA16(delta);
  w0 = BBE_MULNX16(x, v0); // sequence: d, 2*d, ... 7*d
  x = BBE_PACKLNX40(w0);
  v0 = BBE_MOVVA16(phase0);
  x = BBE_ADDNX16(x, v0); // sequence: a0+d, a0+2*d, ... a0+7*d

  addx = BBE_MOVVA16(delta << LOG2_BBE_SIMD_WIDTH);

  tsin0l = BBE_LVNX16_I(pT, 0 * (2 * BBE_SIMD_WIDTH));
  tsin0h = BBE_LVNX16_I(pT, 1 * (2 * BBE_SIMD_WIDTH));
  //tsin1 = BBE_LVNX16_I(pT, 2 * (2 * BBE_SIMD_WIDTH));
  //tcos1 = BBE_LVNX16_I(pT, 3 * (2 * BBE_SIMD_WIDTH));

  for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    // first approximation 
    ix0 = BBE_MOVVSELNX16(x, 12);

    /* OLD:
    v0 = BBE_SRAINX16(x, 9);
    v0 = BBE_ANDNX16(v0, _7);
    ix1 = BBE_MOVVSELNX16(v0, 0);
    */
    ix1 = BBE_MOVVSELNX16(x, 9);

    v0 = BBE_ADDNX16(x, _4000);
    ix01 = BBE_MOVVSELNX16(v0, 12);

        
    yc0 = BBE_SHFLNX16(tsin0l, ix01);
    ys0 = BBE_SHFLNX16(tsin0l, ix0);

    // second approximation
    tcos1 = BBE_SELNX16I(tsin0h, tsin0h, BBE_SELI_EXTRACT_HI_HALVES);
    tsin1 = BBE_SELNX16I(tsin0h, tsin0h, BBE_SELI_EXTRACT_LO_HALVES);
    yc1 = BBE_SHFLNX16(tcos1, ix1);
    ys1 = BBE_SHFLNX16(tsin1, ix1);

    v0 = BBE_SLLINX16(x, 4);
    v0 = BBE_POLYNX16_OFF(v0, 13, 0); // offset in interval

    // last appoximation: ~ x*pi 
    // "yc2": see initialization of loop
    ys2 = BBE_MULNX16PACKQ(v0, _6434);

    // 2 complex multiplies 
    BBE_DSELNX16I(u1, u0, ys0, yc0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(v1, v0, ys1, yc1, BBE_DSELI_INTERLEAVE_1);
    yc1 = BBE_MULNX16CPACKQ(u0, v0); // complex multiplcation of low parts
    ys1 = BBE_MULNX16CPACKQ(u1, v1); // complex multiplcation of hi parts

    BBE_DSELNX16I(u1, u0, ys2, yc2, BBE_DSELI_INTERLEAVE_1);
    v0 = BBE_MULNX16CPACKQ(u0, yc1); // complex multiplcation of low parts
    v1 = BBE_MULNX16CPACKQ(u1, ys1); // complex multiplcation of hi parts

    // save complex coordinates
    BBE_SVNX16_IP(v0, pZ, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(v1, pZ, 2 * BBE_SIMD_WIDTH);

    // update argument
    x = BBE_ADDNX16(x, addx);
  }

  return(BBE_MOVAV16(x));
} /* cbexplin() */

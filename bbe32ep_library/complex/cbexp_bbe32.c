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
Complex Exponential

Description: These functions compute complex number lying on the unit 
circle from input phase data.

Representation:
cbexp         16-bit signed fixed-point format
              Input data phase[] are Q15 angular values normalized by pi.
              Output data z[] format is Q(15-sh), where sh is the shift
              control argument
cbexpf,sexpf  IEEE-754 Std. single precision floating-point format
              Input phase data phase[] are in radians. See the Note.
              
Accuracy:
3.7e-4 (12 in Q15) for cbexp
2 ULP for cbexpf, sexpf
3 ULP for cbfastexpf

Note for (non-fast versions):
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the respective result is 0+0j.

Input domain for 'fast' version cfastbexpf()
|phase|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase[N]  Phase values
sh        (cbexp only) left bit shift amount, [0..15]
N         Length of vectors
Output:
z[N]      Results

Restrictions:
z,phase   Aligned on 32-byte boundary
z,phase   Must not overlap
N         Multiple of 16 (cbexp) or 8 (cbexpf,cbfastexpf)
-------------------------------------------------------------------------*/

void cbexp ( complex_fract16 * restrict z, const int16_t   * restrict phase, int sh, int N )
{
  int n;
  const xb_vecNx16 * restrict pA = (const xb_vecNx16 *)phase;
  const xb_vecNx16 * restrict pT = (const xb_vecNx16 *)cbexpTbl;
  xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
  vsaN sh_v;
  xb_vecNx16 tsin0l, tsin0h, tsin1, tcos1;
  xb_vecNx16 x, yc0, ys0, yc1, ys1, yc2, ys2;
  xb_vecNx16 _4000, _6434;
  xb_vecNx16 v0, v1, u0, u1;
  xb_vecNx40 w0, w1;
  vselN ix0, ix01, ix1;

  NASSERT(0 <= sh && sh <= 15);
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(phase, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH) == 0);

  tsin0l = BBE_LVNX16_I(pT, 0 * (2 * BBE_SIMD_WIDTH));
  tsin0h = BBE_LVNX16_I(pT, 1 * (2 * BBE_SIMD_WIDTH));
  //tsin1 = BBE_LVNX16_I(pT, 2 * (2 * BBE_SIMD_WIDTH));
  //tcos1 = BBE_LVNX16_I(pT, 3 * (2 * BBE_SIMD_WIDTH));

  _4000 = BBE_MOVVA16(0x4000);
  yc2 = BBE_MOVVA16(32763);
  _6434 = BBE_MOVVA16(6434);

  sh += 15;
  sh_v = BBE_MOVVSA32(sh);

  for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    BBE_LVNX16_IP(x, pA, 2 * BBE_SIMD_WIDTH);

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

    x = BBE_SLLINX16(x, 4);
    x = BBE_POLYNX16_OFF(x, 13, 0); // offset in interval

    // last appoximation: ~ x*pi 
    // "yc2": see initialization of loop
    ys2 = BBE_MULNX16PACKQ(x, _6434);

    // 2 complex multiplies 
    BBE_DSELNX16I(u1, u0, ys0, yc0, BBE_DSELI_INTERLEAVE_1);
    BBE_DSELNX16I(v1, v0, ys1, yc1, BBE_DSELI_INTERLEAVE_1);
    yc1 = BBE_MULNX16CPACKQ(u0, v0); // complex multiplcation of low parts
    ys1 = BBE_MULNX16CPACKQ(u1, v1); // complex multiplcation of hi parts

    BBE_DSELNX16I(u1, u0, ys2, yc2, BBE_DSELI_INTERLEAVE_1);
    w0 = BBE_MULRNX16C(u0, yc1, sh_v); // complex multiplcation of low parts
    w1 = BBE_MULRNX16C(u1, ys1, sh_v); // complex multiplcation of hi parts

    v0 = BBE_PACKVNX40(w0, sh_v);
    v1 = BBE_PACKVNX40(w1, sh_v);

    // save complex coordinates
    BBE_SVNX16_IP(v0, pZ, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(v1, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* cbexp() */

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
    Polar to Cartesian Conversion
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
#include "cbexpTbl.h"

/*-------------------------------------------------------------------------
Polar to Cartesian Conversion 

Description: These functions convert input pairs of magnitude and phase
values into Cartesian coordiantes of respective points on the complex plane. 

Representation:
vpolar           16-bit signed fixed-point format
                 Input phase data phase[] are Q15 angular values
                 normalized by pi.
                 Input magnitude data x[] are of abn arbitrary fixed-
                 point format Qx.
                 Fixed-point format for resulting coordinates on the 
                 complex plane z[] is Qx-sh.
vpolarf,spolarf  IEEE-754 Std. single precision floating-point format
                 Input phase data phase[] are in radians.

Accuracy:
9.2e-5 (3 in Q15) - for vpolar
2 ULP - vpolarf, spolarf
3 ULP - vfastpolarf

Note for vpolarf,spolarf:
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the conversion result is (0,0).

Input domain for 'fast' version vfastpolarf()
|a|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
r[N]       Magnitude data
a[N]       Phase data
sh         Right bit shift amount, 0..15 (vpolar)
N          Length of vectors
Output:
z[N]       Points on the complex plane

Restrictions:
z,x,phase  Aligned on 32-byte boundary
z,x,phase  Must not overlap
N          Multiple of 16 (vpolar) or 8 (vpolarf,vfastpolarf)
-------------------------------------------------------------------------*/

void vpolar ( complex_fract16 * restrict z, 
               const int16_t   * restrict x, 
               const int16_t   * restrict phase, 
               int sh, int N)
{
  const xb_vecNx16* restrict pA;
  const xb_vecNx16* restrict pA1;
  const xb_vecNx16* restrict pR;
  xb_vecNx16* restrict pZ;
  const xb_vecNx16* restrict pTbl;
  int n;
  xb_vecNx16 c0, c2, c3;

  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(phase, (2 * BBE_SIMD_WIDTH));
  NASSERT(0 <= sh || sh <= 15);
  NASSERT(N%BBE_SIMD_WIDTH == 0);
  if (N <= 0) return;

  pTbl = (const xb_vecNx16*)(sineTbl + 16);
  NASSERT_ALIGN(pTbl, (2 * BBE_SIMD_WIDTH));
  pR = (const xb_vecNx16*)(x);
  pA = (const xb_vecNx16*)(phase); pA1 = pA;
  pZ = (xb_vecNx16*)(z);

  sh += 15;

  const xb_vecNx16* restrict pA2 = pA;
  const xb_vecNx16* restrict pA3 = pA;

  BBE_LVNX16_IP(c0, pTbl, 4 * BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(c2, pTbl, 2 * BBE_SIMD_WIDTH);
  BBE_LVNX16_XP(c3, pTbl, -3 * 2 * BBE_SIMD_WIDTH);

  __Pragma("loop_count min=1");
  for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    xb_vecNx16 x, sgnx, yc, ys, r;
    vselN idx, idx2;

    BBE_LVNX16_XP(x, pA, 2 * BBE_SIMD_WIDTH);
    x = BBE_ABSNX16(x);
    idx = BBE_MOVVSELNX16(x, 11);
    x = BBE_ADDNX16(x, x);
    x = BBE_POLYNX16_OFF(x, 12, 1);

    ys = BBE_SHFLNX16(c0, idx);
    ys = BBE_ADDNX16(BBE_MULNX16PACKQ(ys, x), BBE_SHFLNX16(c2, idx));
    ys = BBE_ADDSNX16(BBE_MULNX16PACKQ(ys, x), BBE_SHFLNX16(c3, idx));

    BBE_LVNX16_XP(sgnx, pA1, 2 * BBE_SIMD_WIDTH);
    ys = BBE_MULSGNNX16(sgnx, ys);

    // compute cosine
    BBE_LVNX16_XP(x, pA2, 2 * BBE_SIMD_WIDTH);
    x = BBE_ADDNX16(x, BBE_MOVVA16(0x4000));
    x = BBE_ABSNX16(x);
    idx2 = BBE_MOVVSELNX16(x, 11);
    x = BBE_ADDNX16(x, x);
    x = BBE_POLYNX16_OFF(x, 12, 1);

    yc = BBE_SHFLNX16(c0, idx2);
    yc = BBE_ADDNX16(BBE_MULNX16PACKQ(yc, x), BBE_SHFLNX16(c2, idx2));
    yc = BBE_ADDSNX16(BBE_MULNX16PACKQ(yc, x), BBE_SHFLNX16(c3, idx2));

    BBE_LVNX16_XP(sgnx, pA3, 2 * BBE_SIMD_WIDTH);
    sgnx = BBE_ADDNX16(sgnx, BBE_MOVVA16(0x4000));
    yc = BBE_MULSGNNX16(sgnx, yc);

    BBE_LVNX16_XP(r, pR, 2 * BBE_SIMD_WIDTH);
    ys = BBE_PACKVNX40(BBE_MULRNX16(ys, r, sh), sh);
    yc = BBE_PACKVNX40(BBE_MULRNX16(yc, r, sh), sh);

    // interleave result
    BBE_DSELNX16I(ys, yc, ys, yc, BBE_DSELI_INTERLEAVE_1);

    // save complex coordinates
    BBE_SVNX16_IP(yc, pZ, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(ys, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* pvcartesian() */

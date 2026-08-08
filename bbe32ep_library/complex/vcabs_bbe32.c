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
    Complex Vector Magnitude
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
Magnitude Of Complex Number

Description: These functions multiply complex input number by its conjugate
and take the square root of the result.

Representation:
vcabs          16-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-defined.
               Fixed-point format for output data is Qz = Qx+sh/2, where sh
               is the shift control argument. If a resulting value is too
               large to be represented in Qz format, then it is saturated
               by 32767.
vcabsf,scabsf  IEEE-754 Std. single precision floating-point format


Input domain for 'fast' version vfastcabsf():
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
1 LSB - vcabs
2 ULP - vcabsf,scabsf
3 ULP - vfastcabsf

Parameters:
Input:
x[N]  Complex numbers
sh    (vcabs only) bi-directional shift control, an even integer from
      the range [-2,30]. Positive value corresponds to a left shift
N     Length of vectors
Output:
z[N]  Magnitudes

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vcabs) or 8 (vcabsf,vfastcabsf)
-------------------------------------------------------------------------*/

#if !(HAVE_NSAENX40 && HAVE_RECIP && 1)

DISCARD_FUN(void, vcabs, (int16_t   * restrict z, const complex_fract16 * restrict x, int sh, int N))

#else

void vcabs(int16_t   * restrict z, const complex_fract16 * restrict x, int sh, int N)
{
  const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
        xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
  static const ALIGN(32) int16_t SEL[BBE_SIMD_WIDTH] = { 0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15 };
  const xb_vecNx16 * restrict pS = (const xb_vecNx16 *)SEL;

  int n;
  xb_vecNx16 x0, x1, y0;
  xb_vecNx16 z0;
  xb_vecNx16 b, c;
  xb_vecNx40 A0, a;
  vsaN       shft, nsa;
  vselN      sel0;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % BBE_SIMD_WIDTH == 0);
  NASSERT((-2 <= sh) && (sh <= 30));

  x0 = BBE_LVNX16_I(pS, 0 * 2 * BBE_SIMD_WIDTH);
  sel0 = BBE_MOVVSELNX16(x0, 0);

  shft = BBE_MOVVSA32(sh+1);
  z0 = BBE_MOVVINT16(-5);
  for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); ++n)
  {
    BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(x1, pX, 2 * BBE_SIMD_WIDTH);

    A0 = BBE_MAGINX16C(x1, x0);
    A0 = BBE_SLSNX40(A0, shft);
    //A0 = BBE_SLSINX40(A0, 1);
    nsa = BBE_NSAENX40(A0);
    A0 = BBE_SLLNX40(A0, nsa);

    BBE_RSQRTLUNX40_0(a, b, c, A0);
    BBE_RSQRTLUNX40_1(a, b, c, A0);
    BBE_MULUUSNX16(a, c, b);

    BBE_RECIPLUNX40_0(A0, y0, b, a);
    BBE_RECIPLUNX40_1(A0, y0, b, a);

    BBE_MULUSANX16(A0, b, y0);
    b = BBE_MOVVVS(nsa);
    b = BBE_SRAINX16(b, 1);
    b = BBE_ADDNX16(b, z0);
    nsa = BBE_MOVVSVS(b);
    A0 = BBE_SRSNX40(A0, nsa);
    y0 = BBE_PACKHNX40(A0);
    y0 = BBE_SHFLNX16(y0, sel0);
    BBE_SVNX16_IP(y0, pZ, 2 * BBE_SIMD_WIDTH);
  }
} /* cvmag() */

#endif

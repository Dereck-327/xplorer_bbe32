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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
    2/4 Tap Interpolators
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"

/*-------------------------------------------------------------------------
2/4 Tap Interpolators

Representation: 16-bit signed fixed-point

Functions calculate 
y[k]=taps[0]*x1[k]+taps[1]*x2[k];
or 
y[k]=taps[0]*x1[k]+taps[1]*x2[k]+taps[2]*x3[k]+taps[3]*x4[k];
for all elements of input vectors
where taps[] is real vector and x<1..4>[] are complex vectors

Input:
N             size of complex vectors
x1[2*N]       Input complex vectors of length N, Q15
x2[2*N] 
x3[2*N] 
x4[2*N] 
taps[2 or 4]  Taps, Q15
Output:
y[2*N]        Output complex vector of length N, Q15

Return value:
none

Restrictions:
x1,x2,x3,x4,y       Must be aligned on 32-byte boundary
x1,x2,x3,x4,y,taps  Must not overlap
N                   Must be a multiple of 8
-------------------------------------------------------------------------*/

void interp_2tap ( int16_t * restrict y,
             const int16_t * restrict x1,
             const int16_t * restrict x2,
             const int16_t * restrict taps,
             int N )
{
  const xb_vecNx16    * restrict px1 = (xb_vecNx16   *)x1;
  const xb_vecNx16    * restrict px2 = (xb_vecNx16   *)x2;
  xb_vecNx16    * restrict py = (xb_vecNx16   *)y;

  int n;
  xb_vecNx16    t0, t1, vx0, vx1, t;
  xb_vecNx40    A10;

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x1);
  NASSERT_ALIGN32(x2);
  NASSERT((N & 7) == 0);

  t0 = BBE_MOVVA16(taps[0]);
  t1 = BBE_MOVVA16(taps[1]);
  for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); ++n)
  {
    BBE_LVNX16_IP(vx0, px1, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(vx1, px2, 2 * BBE_SIMD_WIDTH);
    A10 = BBE_MULNX16(vx0, t0);
    BBE_MULANX16(A10, vx1, t1);
    t = BBE_PACKQNX40(A10);
    BBE_SVNX16_IP(t, py, 2 * BBE_SIMD_WIDTH);
  }
} /* interp_2tap() */

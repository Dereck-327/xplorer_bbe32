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
    Complex Vector Product
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Dot Product of Complex Vectors

Description: These routines take two complex vectors and calculate their
dot product.

Representation:
cvdot   Signed fixed-point format
        Input vectors x and y are comprised of 32-bit complex data with
        16-bit signed real and imaginary components of arbitrary formats
        Qx and Qy. Complex dot product is computed in a pair of 40-bit 
        accumulators, which are then rounded, shifted to the right by rsh
        bit positions and saturated to form a 64-bit complex result with
        Qx+Qy-rsh fractional bits in 32-bit real and imaginary components.
cvdotf  IEEE-754 Std. single precision floating-point format for real/imaginary
        components of 64-bit input data and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (cvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cvdot) or 4 (cvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/

complex_fract32 cvdot  ( const complex_fract16 * restrict x,
                         const complex_fract16 * restrict y,
                         int rsh,
                         int N )
{
  complex_fract32 ALIGN(64) z; 
  const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pY = (const xb_vecNx16 *)y;
  int k;
  xb_vecNx16 x0, y0;
  xb_vecNx40 A0;
  xb_c40   A;
  vsaN VsaShift;

  if (N <= 0)
  {
      y0 = 0;
      BBE_SV4X16_I(y0, (&z), 0);
      return z;
  }
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N>0 && N % (BBE_SIMD_WIDTH / 2) == 0);
  NASSERT(rsh >= 0);

  VsaShift = BBE_MOVVSA32(rsh - 3);
  A0 = 0;
  for (k = 0; k<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); k++)
  {
      BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_XP(y0, pY, 2 * BBE_SIMD_WIDTH);

      BBE_MULANX16C(A0, x0, y0);
  }
  {
      vsaN vsa = BBE_MOVVSA32(rsh);
      A0 = BBE_RNDADJNX40(A0, VsaShift); // rounding code 
      A = BBE_RADDNX40C(A0);
      A0 = BBE_MOVNX40_FROMC40(A);
      A0 = BBE_SRANX40(A0, vsa);
      y0 = BBE_MOVSVWL(A0);
      BBE_SV4X16_I(y0, (&z), 0);
      return z;
  }
} /* cvdot() */

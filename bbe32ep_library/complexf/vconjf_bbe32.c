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
  NatureDSP_Baseband library. Complex Math functions
    Conjugate Of Complex Number
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, vconjf, ( complex_float * restrict z, const complex_float * restrict x, int N ))
#else

/*-------------------------------------------------------------------------
Conjugate Of Complex Number

Description: These functions negate the imaginary part of complex input 
numbers.

Representation: IEEE-754 Std. single precision floating-point format

Accuracy: Exact

Parameters:
Input:
x[N]  Complex numbers
N     Length of vectors
Output:
z[N]  Conjugated input numbers

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 4
-------------------------------------------------------------------------*/

void vconjf ( complex_float * restrict z, const complex_float * restrict x, int N )
{
  int n;
  const xb_vecN_2xf32 *restrict px;
        xb_vecN_2xf32 *restrict pz;
  xb_vecN_2xf32 in, out;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 4) == 0);

  px = (const xb_vecN_2xf32 *)x;
  pz = (      xb_vecN_2xf32 *)z;

  for (n = 0; n < (N>>(LOG2_BBE_SIMD_WIDTH-2)); n++)
  {
    BBE_LVN_2XF32_IP(in, px, 2*BBE_SIMD_WIDTH);
    out = BBE_CONJN_2XF32(in);
    BBE_SVN_2XF32_IP(out, pz, 2*BBE_SIMD_WIDTH);
  }
} /* vconjf() */

#endif/* !HAVE_VFPU */

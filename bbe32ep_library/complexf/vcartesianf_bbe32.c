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
    Cartesian To Polar Conversion
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
DISCARD_FUN(void, vcartesianf, ( float32_t     * restrict z,
                                 float32_t     * restrict phase,
                           const complex_float * restrict x,
                           int N ))
#else

/*-------------------------------------------------------------------------
Cartesian To Polar Conversion

Description: These functions convert Cartesian coordinates of points on
the complex plane to the polar system (magnitude and phase).

Representation: IEEE-754 Std. single precision floating-point format

Input domain for 'fast' version vfastcartesianf():
|real(x)|>1.1755e-038
|imag(x)|>1.1755e-038
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy: 
2 ULP for vcartesianf(),scartesianf()
3 ULP for vfastcartesianf()

Parameters:
Input:
x[N]       Cartesian coordinates
N          Length of vectors
Output:
z[N]       Magnitude data
phase[N]   Phase data

Restrictions:
z,phase,x  Aligned on 32-byte boundary
z,phase,x  Must not overlap
N          Multiple of 8
-------------------------------------------------------------------------*/

void vcartesianf ( float32_t     * restrict z,
                   float32_t     * restrict phase,
             const complex_float * restrict x,
             int N )
{
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(phase, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  /* Find the argument of complex numbers */
  vargf(phase, x, N);
  /* Find the magnitude of complex numbers */
  vcabsf(z, x, N);
} /* vcartesianf() */

#endif/* !HAVE_VFPU */

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
    Inverse Square Root
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"

/*-------------------------------------------------------------------------
Reciprocal Square Root

Description: Evaluate the reciprocal square root of input value x and store
result to y: y = 1/x^0.5.

Representation:
vrsqrt,srsqrt  32-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-
               defined, provided that it is even.
               Fixed-point format Qy for output data is Qy = 31-Qx/2.
               For example, if Qx == 30 then Qy = 31-30/2 = 16
vrsrtf         IEEE-754 Std. single precision floating-point format

Accuracy:
vrsqrt,srsqrt  1.1e-4 (worst case relative error)
vrsrtf,srsqrtf 2 ULP
vfastrsrtf     3 ULP

Notes (for non-fast versions)::
1. Fixed-point functions return MIN_INT32 (0x80000000) for a negative or
   zero input value.
2. Floating-point reciprocal square root conforms to IEEE-754 Std rSqrt 
   operation in respect to signaling error conditions by means of floating-
   point exceptions.
3. Floating-point reciprocal square root limits the range of allowable
   input values, as follows:
   A) If x<0, functions raise the "invalid" floating-point exception,
      assign EDOM to errno and set output value y to NaN.
   B) If x==+/-0, functions set output value y to +/-HUGE_VALF, raise the
      "divide by zero" floating-point exception, and assign ERANGE to errno.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vrsqrt) or 8 (vrsqrtf,vfastrsqrtf)
-------------------------------------------------------------------------*/

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)

DISCARD_FUN(int32_t, srsqrt, (int32_t x))

#else

int32_t srsqrt(int32_t x)
{
    xb_vecNx40   a_vec, x_vec, z0, z1;
    xb_vecNx16   tmp, exp;
    xb_vecNx16   /*x0, x1,*/ b_vec, cn_vec, _2 = BBE_MOVVINT16(2 + 15 + 9);
    vsaN c_vec;
    int32_t      y;
    vboolN    b0;
    z0 = 0;
    z1 = BBE_MOVWA32(0x80000000);

    x_vec = BBE_MOVWA32(x);
    b0 = BBE_LTNX40(x_vec, z0);
    x_vec = BBE_SLLINX40(x_vec, 1);
    c_vec = BBE_NSAENX40(x_vec);
    x_vec = BBE_SLLNX40(x_vec, c_vec);
    BBE_RSQRTLUNX40_0(a_vec, b_vec, cn_vec, x_vec);
    BBE_RSQRTLUNX40_1(a_vec, b_vec, cn_vec, x_vec);
    BBE_MULUUSNX16(a_vec, cn_vec, b_vec);
    tmp = BBE_MOVVVS(c_vec);
    tmp = BBE_SRAINX16(tmp, 1);
    exp = BBE_SUBNX16(tmp, _2);
    c_vec = BBE_MOVVSV(exp, 0);
    a_vec = BBE_SLSNX40(a_vec, c_vec);
    a_vec = BBE_MOVNX40T(z1, a_vec, b0);
    y = BBE_MOVAW32(a_vec);
    return y;
} /* srsqrt() */

#endif 

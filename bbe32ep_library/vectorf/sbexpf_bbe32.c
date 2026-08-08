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
  NatureDSP_Baseband library. Vector Operations
    Common Block Exponent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Vector Operations. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Common Block Exponent

Description: These functions compute base-2 exponent adjustment term needed
to normalize data in the input vector. Exact meaning of normalization depends
on the data format (see below).

Representation: 
vbexp,vbexp_fast        16-bit signed fixed-point format
                        For each input value functions count the number of
                        redundant sign bits (as if the value was loaded in
                        a 32-bit register) and return the minimum result over
                        the input data vector.
sbexp                   32-bit signed fixed-point format
                        Count the number of redundant sign bits and return 
                        the result.
vbexpf,sbexpf           IEEE-754 Std. single precision floating-point format
                        For each finite input value x, functions estimate the
                        integer E(x), such that 0.5 <= |x|*2^E(x) < 1 and
                        -128 <= E(x) <= 148. The minimum value of E(x) over
                        input data vector is the result.

Special cases:
   x    |  Result |    Extra Conditions    
--------+---------|---------------------------
0       |    0    |
+/-Inf  | -129    | floating-point functions
NaN     |    0    |
--------|---------|---------------------------
0       |   31    |
-32768  |   16    | fixed-point functions
32767   |   16    |

Parameters:
Input:
x[N]    Input data
N       Length of data vector
Returned Value: exponent adjustment term, or zero if N<=0
Restrictions:
vbexp(), vbexpf():
  No restrictions
vbexp_fast():
  x     Aligned on 32-byte boundary
  N     Multiple of 16
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(int,sbexpf,( float32_t x ))
#else
int sbexpf ( float32_t x )
{
  /*
    ufloat32uint32 ux;
  float32_t xx;
  int e;
  xx = fabsf(x);
  if (isnan(xx) || x == 0.f) return 0;
  ux.f = xx;
  e = 0;
  if ((ux.u >> 23) == 0)  // multiply denorm numbers by 2^24 
  {
    xx = xx*16777216.f;
    e -= 24;
  }
  ux.f = xx;
  e += ((int)(ux.u >> 23) - 126);
  return -e; 
  */

  unsigned e;
  xtfloat f0;
  int32_t u0, exp0;
  xtbool b0, b1;
  xtbool zero, bnan;
  f0 = x;
  f0 = XT_ABS_S(f0);
  b0 = XT_UN_S(f0, f0);
  b1 = XT_OEQ_S(f0, 0.f);
  zero = (b0);
  bnan = (b1);
  zero = XT_ORB(zero, bnan);
  exp0 = 0;
  u0 = XT_RFR(f0);
  u0 = ((uint32_t)u0) >> 23;
  if (u0 == 0)
  {
    f0 = XT_MUL_S(f0, 16777216.f);
    exp0 -= 24;
  }
  u0 = XT_RFR(f0);
  u0 = ((uint32_t)u0) >> 23;
  exp0 += ((int)(u0)-126);
  e = exp0;
  XT_MOVT(e, 0, zero);
  return -(int32_t)e;

} /* sbexpf() */
#endif

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
  NatureDSP_Baseband library. Math functions
    Logarithms
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include <errno.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* tables */
#include "logf_tbl.h"
#include "sqrt2f_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"

/*-------------------------------------------------------------------------
Logarithms

Description: These function compute base-2, base-10 or natural logarithm of
input data.

Representation:
vlog2,vlogn,vlog10,     Signed fixed-point format
slog2,slogn,slog10      Input data are 32-bit Q16.15, results are 16-bit Q4.11.
                        Here are a few examples for the base-2 logarithm:

                          Function | Input Data Q16.15 <real> | Result Q4.11 <real>
                        -----------+--------------------------+--------------------
                           slog2   | 65536 <2.0>              | 2048 <1.0>
                           slog2   | 2147483647 <65535.99997> | 32767 <15.9995>
                           slog2   | 1 <3.052e-5>             | -30720 <-15.0>
                        -----------+--------------------------+--------------------
vlog2f,vlognf,vlog10f,  IEEE-754 Std. single precision floating-point format
slog2f,slognf,slog10f

Accuracy:
1 LSB for the fixed-point functions
2 ULP for the floating-point functions

Notes:
1. Fixed-point Functions return -32768 for a negative or zero input.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating point functions limit the range of allowable input values:
   A) If x<0, the result is set to NaN, errno is assigned the value EDOM, and
      "invalid" floating-point exception is raised
   B) If x==0, the result is set to minus infinity, errno is assigned the value 
      ERANGE, and "divide-by-zero" floating-point exception is raised

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vlog2,vlogn,vlog10) or 8 (vlog2f,vlognf,vlog10f)
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t slognf ( float32_t x )
{
  float32_t one, two, zero, inf;
  float32_t y, y1, y2, y3, x2, x4;
  vbool1 b_xinf, b_subn, b_small;
  unsigned ux;
  int e;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */
  zero = XT_CONST_S(0);
  one = XT_CONST_S(1);
  two = XT_CONST_S(2);
  inf = plusInff.f;

  /* Check for out-of-domain values: x==NaN, x==0, x<0 */
  if (BBE_MOVAB1(XT_ULE_S(x, zero)))
  {
    vbool1 bzero;
    int err;
    __Pragma("frequency_hint never");
    bzero = XT_OEQ_S(x, zero);
    err = EDOM; XT_MOVNEZ(err, ERANGE, BBE_MOVAB1(bzero));
    errno = err;
    /* generate FE_DIVBYZERO and FE_INVALID with right output (-inf or qNaN) */
    x = XT_RSQRT0_S(XT_NEG_S(XT_ABS_S(x)));
    return x;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
 // BBE_MOVSCFV(XTENSA_SCF_SET_EXCEPTION_ENABLE(SCF, 0)); /* Clear FP flags */

  b_xinf = XT_OLE_S(inf, x);
  /* implementation of frexpf */
  b_subn = XT_OLT_S(x, realminf.f);
  y = zero; XT_MADDN_S(y, x, 8388608.f); /* protect from OVFL exception */
  XT_MOVT_S(x, y, b_subn);
  e = 0; XT_MOVNEZ(e, -23, BBE_MOVAB1(b_subn));
  ux = XT_RFR(x);
  e += (ux >> 23) - 126;
  ux &= 0x807FFFFF;
  ux |= 0x3F000000;
  x = XT_WFR(ux);
  /*--------------------------*/

  b_small = XT_OLT_S(x, sqrt0_5f.f);
  XT_MOVT_S(x, x*two, b_small);
  XT_MOVNEZ(e, e - 1, BBE_MOVAB1(b_small));

  x = XT_SUB_S(x, one);
  x2 = XT_MUL_S(x, x);
  x4 = XT_MUL_S(x2, x2);
  /* Combination of Estrin`s and Horner schemes */
  y1 = logf_tbl[1].f; XT_MSUB_S(y1, x, logf_tbl[0].f);
  y2 = logf_tbl[3].f; XT_MSUB_S(y2, x, logf_tbl[2].f);
  y3 = logf_tbl[5].f; XT_MSUB_S(y3, x, logf_tbl[4].f);
  y = logf_tbl[7].f; XT_MSUB_S(y, x, logf_tbl[6].f);

  XT_MADD_S(y2, y1, x2);
  XT_MADD_S(y, y3, x2);
  XT_MADD_S(y, y2, x4);
  XT_MADD_S(x, x2, y); y = x;
  XT_MADD_S(y, (float32_t)e, ln2.f);
  XT_MOVT_S(y, inf, b_xinf); /* set +Inf for x=+Inf */

  BBE_MOVSCFV(SCF);

  return (y);

} /* slognf() */
#else
DISCARD_FUN(float32_t,slognf,( float32_t x ))
#endif

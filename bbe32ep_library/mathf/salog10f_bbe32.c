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
    Antilogarithm and Exponential
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
/* Tables */
#include "expf_tbl.h"
#include "alog10f_tbl.h"
#include "pow2f_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"
/*-------------------------------------------------------------------------
Antilogarithm and Exponential

Description: These functions compute base-10 or natural antilogarithm of
input data.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 2 ULP.

Note:
Antilogarithm functions conform to ANSI C requirements on standard math 
library functions in respect to treatment of errno and floating-point
exceptions.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t salog10f ( float32_t x )
{
  float32_t y, dy;
  float32_t p0;
  float32_t y1;
  xb_int32v e1, e2;
  /* Is a NaN or is less than zero; is equal to zero. */
  vbool1 b_max, b_nan, b_inf;

  int err;
  xb_int32v SCF; /* Floating-point Status and Control Register values. */

  /* Check input for values that are out of domain/range */
  b_nan = XT_UN_S(x, x);                  /* x==NaN               */
  b_max = XT_ULE_S(alog10fminmax[1].f, x);/* x>=38.5318 or x==NaN */
  b_inf = XT_UEQ_S(x, plusInff.f);        /* x==+Inf or x==NaN    */

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */
  if ((xtbool)(XT_OLE_S(x, alog10fminmax[0].f)))
  {
    BBE_MOVSCFV(SCF);
    return 0.f;
  }
  err = errno;
  XT_MOVT(err, EDOM, b_nan);

  /*
  * Multiply by 1/log10(2)
  */
  p0 = XT_MUL_S(x, log2_10[0].f);
  p0 = XT_FIROUND_S(p0);
  dy = XT_NEG_S(p0);
  XT_MADD_S(dy, x, log2_10[0].f);
  XT_MADD_S(dy, x, log2_10[1].f);
  XT_MOVT_S(dy, 0.f, b_inf);

  /* compute 2^x */
  /* Load and replicate polynomial coefficients */

  y = pow2f_coef[0].f; y1 = pow2f_coef[1].f;
  XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[2].f;
  XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[3].f;
  XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[4].f;
  XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[5].f;
  XT_MADD_S(y1, dy, y); y = y1; y1 = pow2f_coef[6].f;
  XT_MADD_S(y1, dy, y); y = y1; 
  /* resulted scaling */

  p0 = XT_MAX_S(XT_MIN_S(p0, 129.f), -151.f);
  /* Apply exponential part to the result */
  {
    xb_vecNx16 tmp, v1, v2;
    tmp = BBE_MOVNX16_FROMN_2X32(BBE_TRUNCN_2XF32(BBE_MOVN_2XF32_FROMF32(p0), 0));
    tmp = BBE_ADDNX16(tmp, BBE_MOVVA16(254));
    v1 = BBE_SRLINX16(tmp, 1);
    v2 = BBE_SUBNX16(tmp, v1);
    v1 = BBE_SLLINX16(v1, 7);
    v2 = BBE_SLLINX16(v2, 7);
    v1 = BBE_SELNX16I(v1, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
    v2 = BBE_SELNX16I(v2, BBE_ZERONX16(), BBE_SELI_INTERLEAVE_1_EVEN);
    e1 = BBE_MOV32_FROMN_2X32(BBE_MOVN_2X32_FROMNX16(v2));
    e2 = BBE_MOV32_FROMN_2X32(BBE_MOVN_2X32_FROMNX16(v1));
  }

  y = XT_MUL_S(y, BBE_MOVF32_FROM32(e2));
  y = XT_MUL_S(y, BBE_MOVF32_FROM32(e1));
  XT_MOVT_S(y, qNaNf.f, b_nan);
  BBE_MOVSCFV(SCF);
  /* set errno to ERANGE if x>=38.5318 but not +INF */
  if ((xtbool)BBE_ANDNOTB1(b_max, b_inf))
  {
    __Pragma("frequency_hint never");
    err = ERANGE;
    __feraiseexcept(FE_OVERFLOW);
  }
  errno = err;
  return y;

} /* salog10f() */
#else
DISCARD_FUN(float32_t,salog10f,( float32_t x ))
#endif

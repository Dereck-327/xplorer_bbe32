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
    Logarithms
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"
#include "log2Tbl.h"

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

int16_t slog2 ( int32_t x )
{
    int16_t y;
    const xb_vecNx16 * restrict TBL;
    xb_vecNx16 z0, z1, z3;
    xb_vecNx16 y0;
    xb_vecNx16 t0, t1, t2, t;
    xb_vecNx16 v0, e0, a0, b0;
    xb_vecNx40 xw, z2;
    vsaN       nsa;
    vselN      sel;
    vboolN     mask;
    TBL = (const xb_vecNx16*)log2Tbl;

    t0 = BBE_LVNX16_I(TBL, 0);
    t1 = BBE_LVNX16_I(TBL, 2 * BBE_SIMD_WIDTH);
    t2 = BBE_LVNX16_I(TBL, 4 * BBE_SIMD_WIDTH);

    z0 = BBE_MOVVINT16(24);
    z1 = BBE_MOVPINT16(2);
    z2 = 0;
    z3 = BBE_MOVVINX16(BBE_MOVVI_Q15_M1);

    // Q24.15
    xw = BBE_MOVWA32(x);
    nsa = BBE_NSANX40(xw);
    mask = BBE_LTNX40(z2, xw);
    // Q4.12
    xw = BBE_SLLNX40(xw, nsa);
    xw = BBE_SRAINX40(xw, 22);
    v0 = BBE_PACKLNX40(xw);
    e0 = BBE_MOVVVS(nsa);
    // 24 - nsa
    e0 = BBE_SUBNX16(z0, e0);
    sel = BBE_MOVVSELNX16(v0, 12);//Q0 
    // ---- POLY.SU --- //
    v0 = BBE_POLYNX16_OFF(v0, 12, 0);//Q12

    t = BBE_MULNX16PACKQ(t2, v0);
    t = BBE_ADDNX16(t1, t);
    b0 = BBE_MULNX16PACKQ(t, v0);
    a0 = BBE_ADDNX16(t0, b0);

    // Q11 <- Q11 + Q0*Q11  
    b0 = BBE_MULNX16PACKL(e0, z1);
    a0 = BBE_ADDNX16(a0, b0);
    a0 = BBE_SHFLNX16(a0, sel);

    y0 = BBE_MOVNX16T(a0, z3, mask);
    y = BBE_MOVAV16(y0);
    return y;
} /* slog2() */

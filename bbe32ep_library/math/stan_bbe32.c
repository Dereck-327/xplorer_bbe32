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
    Tangent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"
#include "tanTbl.h"

/*-------------------------------------------------------------------------
Tangent 

Description: These functions compute tangent of input data

Representation:
vtan,stan    Signed fixed-point format
             Input data are 16-bit Q15 angular values normalized by pi,
             i.e. fixed-point functions actually compute tan(pi*x).
             Output data are 32-bit Q16.15 values.
vtanf,stanf  IEEE-754 Std. single precision floating-point format for
             input/output data. Input data are treated as angular values
             specified in radians. Floating-point functions limit the
             rangw of allowable input values, see note 3.

Accuracy:
For the fixed-point functions, accuracy depends on the input value x,
as shown in the table below:
   Range of |x|    | Absolute error  | Relative error
-------------------+-----------------+----------------
 [-pi/4; pi/4]     |     1 (Q15)     |
 [pi/4; 7pi/16]    |    15 (Q15)     |    4.6e-4
 [7pi/16; 31pi/64] |   242 (Q15)     |    1.5e-3
-------------------+-----------------+----------------
2 ULP for vtanf(),stanf()
3 ULP for vfasttanf()

Notes for non-fast versions:
1. Fixed-point function result is not defined if input value x is
   +/-pi/2 (+/-8192 in Q15 normalized by pi).
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating-point functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfasttanf():
|x|<804.2477
The output value is not defined outside of this range or accuracy is 
degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vtan) or 8 (vtanf,vfasttanf)
-------------------------------------------------------------------------*/

#if !(HAVE_VSAMATH && HAVE_RECIP && 1)

DISCARD_FUN(int32_t, stan, (int16_t x))

#else

int32_t stan(int16_t x)
{
    const void * restrict TBL;
    xb_vecNx16 x0, xa, ofs;
    xb_vecNx16 z0, z3;
    xb_vecNx16 p0, p1, p2, p3, t;
    xb_vecNx16 y0, y1, b;
    xb_vecNx16 t0, t1, t2, t3;
    xb_vecNx40 yw, a, r0, r1, z4;
    vboolN     b0, b1, b2;
    vselN      idx;
    vsaN       nsa;

    TBL = (const xb_vecNx16*)tanTbl;
    // Order 0 derivatives, signed Q15
    BBE_LVNX16_XP(t0, TBL, 2 * BBE_SIMD_WIDTH);
    // Order 1 derivatives, signed Q18
    BBE_LVNX16_IP(t1, TBL, 2 * BBE_SIMD_WIDTH);
    // Order 2 derivatives, signed Q21
    BBE_LVNX16_IP(t2, TBL, 2 * BBE_SIMD_WIDTH);
    // Order 3 derivatives, signed Q24
    BBE_LVNX16_XP(t3, TBL, -3 * 2 * BBE_SIMD_WIDTH);

    z0 = 0;
    z3 = BBE_MOVPINT16(16);
    z4 = BBE_MOVWA32(0x7fffffff);

    x0 = BBE_MOVVA16(x);
    x0 = BBE_SLLINX16(x0, 1);
    xa = BBE_ABSNX16(x0);
    b0 = BBE_LTNX16(x0, z0);
    b1 = BBE_LTUNX16(z3, xa);
    xa = BBE_SLLINX16(x0, 1);
    xa = BBE_ABSSNX16(xa);

    /*POLI_SU*/
    xa = BBE_SLLINX16(xa, 1);
    ofs = BBE_POLYNX16_OFF(xa, 12, 1); // offset in interval
    idx = BBE_MOVVSELNX16(xa, 12);     // interval number (index in sine table)

    //////////////////////////
    p0 = BBE_SHFLNX16(t0, idx);
    p1 = BBE_SHFLNX16(t1, idx);
    p2 = BBE_SHFLNX16(t2, idx);
    p3 = BBE_SHFLNX16(t3, idx);

    y0 = BBE_MULNX16PACKQ(p0, ofs);
    t = BBE_ADDSNX16(y0, p1); // NOTE : with saturation
    y0 = BBE_MULNX16PACKQ(t, ofs);
    y0 = BBE_ADDSNX16(y0, p2);// NOTE : with saturation
    t = BBE_MULNX16PACKQ(y0, ofs);
    y0 = BBE_ADDSNX16(p3, t);
    b2 = BBE_EQNX16(y0, z0);
    r0 = BBE_UNPKSNX16(y0);
    /* r0 = tan(x), if x in [-pi/4;pi/4]*/
    /*------------------------------*/
    nsa = BBE_NSANX16(y0);
    y1 = BBE_ABSNX16(y0);
    // Q(15+nsa) <- Q15 + nsa;
    y1 = BBE_SLANX16(y1, nsa);
    yw = BBE_UNPKSNX16(y1);
    yw = BBE_SLLINX40(yw, 24);

    BBE_RECIPLUNX40_0(a, y1, b, yw);
    BBE_RECIPLUNX40_1(a, y1, b, yw);

    BBE_MULUSANX16(a, b, y1);
    nsa = BBE_SUBSAVSN(22, nsa);
    a = BBE_SRANX40(a, nsa);
    a = BBE_MOVNX40T(z4, a, b2);
    r0 = BBE_MOVNX40T(a, r0, b1);
    r1 = BBE_NEGNX40(r0);
    r0 = BBE_MOVNX40T(r1, r0, b0);
    return BBE_MOVAW32(r0);
} /* stan() */

#endif

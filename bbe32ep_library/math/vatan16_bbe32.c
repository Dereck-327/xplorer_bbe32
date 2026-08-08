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
    Arctangent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"
#include "atanTbl.h"

/*-------------------------------------------------------------------------
Arctangent 

Description: These functions compute the principal value of arctangent.

Representation:
vatan16,satan16  16-bit signed fixed-point format
                 Input data are Q15. Functions compute atan(x)/(pi/4)
                 and output results in Q15 format.
vatanf,satanf    IEEE-754 Std. single precision floating-point format
vfastatanf       Functions compute atan(x) and output results in radians

Special cases:
    Input | Result 
   -------+--------
    +inf  |  pi/2  (floating-point functions)
    -inf  | -pi/2  

Accuracy:
1 LSB for fixed point fixed point functions
1 ULP for vatanf(), satanf()
2 ULP for vfastatanf()

Notes:
1. These functions are much faster than full-quadrant arctangent atan2,
   so they are preferable when the full phase is not required.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.

Input domain for 'fast' version vfastatanf():
|x|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vatan16) or 8 (vatanf,vfastatanf)
-------------------------------------------------------------------------*/

void vatan16 ( int16_t* restrict z, 
         const int16_t* restrict x, 
         int N )
{
    int n;
    const xb_vecNx16 * restrict TBL = (const xb_vecNx16*)atanTbl;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pX1 = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;
    xb_vecNx16 t0, t1, t2;
    xb_vecNx16 x0, y0, ofs, SignX;
    xb_vecNx16 p0, p1, p2, t;
    xb_vecNx40 w0, z0;
    vselN      sel;

    if (N <= 0) return;
    NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    // Order 0 derivatives, signed Q15
    BBE_LVNX16_XP(t0, TBL, 2 * BBE_SIMD_WIDTH);
    // Order 1 derivatives, signed Q19
    BBE_LVNX16_XP(t1, TBL, 2 * BBE_SIMD_WIDTH);
    // Order 2 derivatives, signed Q23
    BBE_LVNX16_XP(t2, TBL, -2 * 2 * BBE_SIMD_WIDTH);
    z0 = BBE_MOVWA32(0);

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);

        //SignX = x0; //store sign

        x0 = BBE_ABSSNX16(x0);
        x0 = BBE_ADDNX16(x0, x0); // x0 = BBE_SLLINX16(x0, 1);

        /*POLI_SU*/
        // Mid-range offsets
        sel = BBE_MOVVSELNX16(x0, 12);
        ofs = BBE_POLYNX16_OFF(x0, 12, 0);

        
        // Select Taylor expansion coefficients
        p0 = BBE_SHFLNX16(t0, sel);// Q15
        p1 = BBE_SHFLNX16(t1, sel);// Q19
        p2 = BBE_SHFLNX16(t2, sel);// Q23

        // Perform Taylor series summation
        w0 = BBE_MULNX16(p0, ofs);
        w0 = BBE_ADDNX40(w0, z0);
        t = BBE_PACKQNX40(w0);
        t = BBE_ADDNX16(p1, t);
        t = BBE_MULNX16PACKQ(t, ofs);
        t = BBE_ADDSNX16(p2, t);

        BBE_LVNX16_XP(SignX, pX1, 2 * BBE_SIMD_WIDTH);
        y0 = BBE_MULSGNNX16(SignX, t); // correct sign

        BBE_SVNX16_XP(y0, pZ, 2 * BBE_SIMD_WIDTH);
    }
} /* vatan16() */

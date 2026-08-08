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
    Full Arctangent
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
Full Arctangent 

Description: These functionx calculate four quadrant arctangent of complex 
argument x and return arg(x)/pi in Q15 format. Q15 is also assumed for
argument. If real and imaginary components are zero, the result is also zero.

Absolute phase accuracy: 1 LSB 

Parameters:
Input:
x[2*N]  Complex input data. Real and imaginary parts are interleaved with
        real parts stored at even indices.
N       Length of vectors
Output:
z[N]    Result

Restrictions:
z,x     Aligned on 32-byte boundary
z,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/

#if !(HAVE_DIV && 1)

DISCARD_FUN(void, vatan2_16, (int16_t * restrict z,
    const int16_t * restrict x,
    int N))

#else

void vatan2_16(int16_t* restrict z,
           const int16_t* restrict x, 
           int N )
{
    xb_vecNx16 x0, y0, x1, y1, t, z0;
    xb_vecNx16 c0, c2, c4000;
    xb_vecNx40 qw, c8000;
    vboolN rneg, ineg, ygtx, bZero, bZero2;

    const xb_vecNx16 * restrict TBL = (const xb_vecNx16*)atanTbl;
    xb_vecNx16 t0, t1, t2;
    xb_vecNx16 ofs;
    xb_vecNx16 p0, p1, p2;
    vselN      sel;

    int n;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pZ = (xb_vecNx16 *)z;

    NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    c0 = BBE_MOVVA16(0);
    c2 = BBE_MOVVA16(2);
    c4000 = BBE_MOVVA16(0x4000);
    c8000 = BBE_MOVWA40(0, 0x8000);

    BBE_LVNX16_XP(t0, TBL, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(t1, TBL, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(t2, TBL, -2 * 2 * BBE_SIMD_WIDTH);
    //__Pragma("loop_count min=1");
    for (n = 0; n <(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_XP(x1, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y1, pX, 2 * BBE_SIMD_WIDTH);

        // separate real and imaginary parts
        x0 = BBE_SELNX16I(y1, x1, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        y0 = BBE_SELNX16I(y1, x1, BBE_SELI_EXTRACT_1_OF_2_OFF_1);

        bZero = BBE_EQNX16(x0, c0);
        bZero2 = BBE_EQNX16(y0, c0);
        bZero = BBE_ANDB(bZero, bZero2); // flag: re == 0 && im == 0

        // Take signs, and reduce the task to 0 <= y < x, so that 0 <= y/x <= 1 and 0 <= atan(y/x) <= pi/4
        rneg = BBE_LTNX16(x0, c0);
        ineg = BBE_LTNX16(y0, c0);

        x0 = BBE_ABSSNX16(x0);
        y0 = BBE_ABSSNX16(y0);

        // exchange X and Y if need
        ygtx = BBE_LTNX16(x0, y0);
        x1 = BBE_MOVNX16T(y0, x0, ygtx);
        y1 = BBE_MOVNX16T(x0, y0, ygtx);

        // Use 32 by 16 signed division. Saturation comes into effect when x==y. Q15 <- ( Q15 + 15 )/Q15 
        qw = BBE_UNPKQNX16(y1);
        t = BBE_QUONX32(qw, x1);

        // --- calculate z = atan(t); ---
        //t = BBE_ABSSNX16(t);
        t = BBE_ADDNX16(t, t);  // t = BBE_SLLINX16(t, 1);

        //POLI_SU
        // Mid-range offsets
        sel = BBE_MOVVSELNX16(t, 12);
        ofs = BBE_POLYNX16_OFF(t, 12, 0);

        // Select Taylor expansion coefficients
        
        p0 = BBE_SHFLNX16(t0, sel);// Q15
        p1 = BBE_SHFLNX16(t1, sel);// Q19
        p2 = BBE_SHFLNX16(t2, sel);// Q23

        // Perform Taylor series summation
        t = BBE_MULNX16PACKQ(p0, ofs);
        t = BBE_ADDNX16(p1, t);
        t = BBE_MULNX16PACKQ(t, ofs);
        z0 = BBE_ADDSNX16(p2, t);
        // ------------------------------

        //z = (S_add_ss(z, 2) >> 2);
        z0 = BBE_ADDSNX16(z0, c2);
        z0 = BBE_SRAINX16(z0, 2);

        // if |im|>|re| then z <- 1/2 - z 
        t = BBE_SUBNX16(c4000, z0);
        z0 = BBE_MOVNX16T(t, z0, ygtx);

        // Restore the input quadrant. 
        qw = BBE_UNPKUNX16(z0);
        qw = BBE_SUBNX40(c8000, qw);
        t = BBE_PACKSNX40(qw);
        z0 = BBE_MOVNX16T(t, z0, rneg);
        BBE_NEGNX16T(z0, z0, ineg);

        // check special argument
        z0 = BBE_MOVNX16T(c0, z0, bZero);

        BBE_SVNX16_IP(z0, pZ, 2 * BBE_SIMD_WIDTH);
    }
} /* vatan2_16() */

#endif

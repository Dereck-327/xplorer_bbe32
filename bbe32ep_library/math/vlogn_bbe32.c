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
#include "lognTbl.h"

#ifndef BBE_ADDSAVSN
#define ADDSAVSN(z,x,y)      \
{                            \
    xb_vecNx16 xx, yy;       \
    xx = x;                  \
    yy = y;                  \
    xx = BBE_ADDNX16(xx, yy); \
    z = xx;                  \
}
#else
#define ADDSAVSN(z,x,y) z=BBE_ADDSAVSN(x,y);
#endif

#ifndef BBE_SUBSAVSN
#define SUBSAVSN(z,x,y)      \
{                            \
    xb_vecNx16 xx, yy;       \
    xx = x;                  \
    yy = y;                  \
    xx = BBE_SUBNX16(xx, yy); \
    z = xx;                  \
}
#else
#define SUBSAVSN(z,x,y) z=BBE_SUBSAVSN(x,y);
#endif

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

void vlogn ( int16_t * restrict y,
       const int32_t * restrict x,
       int N )
{
    int n;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;
    const xb_vecNx16 * restrict TBL;
    xb_vecNx16 z1;
    xb_vecNx16 x0, x1, y0;
    xb_vecNx16 t0, t1, t2;
    xb_vecNx16 p0, p1, p2;
    xb_vecNx16 v0, e0, a0, b0;
    xb_vecNx40 xw, z2;
    vsaN       nsa /*, shft*/;
    vselN      sel;
    vboolN     mask;

    if (N <= 0) return;
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    TBL = (const xb_vecNx16*)logNatTbl;

    BBE_LVNX16_IP(t0, TBL, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(t1, TBL, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(t2, TBL, 2 * BBE_SIMD_WIDTH);

    z1 = BBE_MOVVA16(22713);// ln2, Q15
    z2 = 0;
    //shft = BBE_MOVVSA32(4);
    y0 = BBE_MOVVINX16(BBE_MOVVI_INT16_MININT);

    __Pragma("ymemory(pX)")

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * 2 * BBE_SIMD_WIDTH);
        x1 = BBE_LVNX16_I(pX, -2 * BBE_SIMD_WIDTH);
        // Q24.15
        xw = BBE_MOVSWV(x1, x0);
        nsa = BBE_NSANX40(xw);

        mask = BBE_LTNX40(z2, xw);
        // Q4.12
        {
            vsaN vsa22;
            ADDSAVSN(vsa22, -22, nsa);
            xw = BBE_SLANX40(xw, vsa22);
            v0 = BBE_PACKLNX40(xw);
        }
        SUBSAVSN(nsa, 24, nsa);
        e0 = BBE_MOVVVS(nsa);
        sel = BBE_MOVVSELNX16(v0, 12);//Q0 
        // ---- POLY.SU --- //
        v0 = BBE_POLYNX16_OFF(v0, 12, 0);//Q12

        p2 = BBE_SHFLNX16(t2, sel);//Q11
        p1 = BBE_SHFLNX16(t1, sel);//Q14
        p0 = BBE_SHFLNX16(t0, sel);//Q17

        // Q14 <- Q14 + ( Q17*Q12 - 15 w/ rounding )
        p2 = BBE_MULNX16PACKQ(p2, v0);
        p1 = BBE_ADDNX16(p1, p2);
        // Q11 <- Q11 + ( Q14*Q12 - 15 w/ rounding )
        p1 = BBE_MULNX16PACKQ(p1, v0);
        a0 = BBE_ADDSNX16(p0, p1);

        // Q15 <- Q0*Q15
        //xw = BBE_MULRNX16(e0, z1, shft);
        // Q11 <- Q15-4
        //b0 = BBE_PACKVNX40(xw, shft);
        e0 = BBE_SLLINX16(e0, 6);
        b0 = BBE_MULNX16PACKP(e0, z1);

        a0 = BBE_ADDNX16(a0, b0);
        a0 = BBE_MOVNX16T(a0, y0, mask);

        BBE_SVNX16_IP(a0, pY, 2 * BBE_SIMD_WIDTH);
    }
} /* vlogn() */

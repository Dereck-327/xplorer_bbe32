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
    Angle Unwrapping
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

#ifdef BBE_PACKNVNX40
#define USE_ROUNDING 1 
#else
#define USE_ROUNDING 0  
#endif

/*-------------------------------------------------------------------------
Angle Unwrapping

Description: These Functions correct phase angles to produce smooth phase 
plots. They adjust phase angles from the input vector by adding a multiple
of 2pi when the absolute difference between two consecutive elements of the
input vector is greater than or equal to the jump tolerance. Tolerance values
smaller than pi are equivalent to the tolerance value of pi.

Representation:
aunwrap   16-bit signed fixed-point format
          It is assumed that all real-world angles are normalized by pi when 
          converted to a fixed-point format. For example, -pi represented
          in Q14 format would be -1*2^14 == -16384.
          Fixed-point format of input phase angles x[N] is Q15. 
          To accommodate the output data format to possible 2pi cycles, the
          number of fractional bits for unwrapped angles y[N] is specified
          through the q argument.
          Fixed-point format for the jump tolerance value tol is Q14.
aunwrapf  IEEE-754 Std. single precision floating-point format for input/output
          data

Parameters:
Input:
x[N]      Original phase angles
tol       Jump tolerance. Typical value is pi (16384 in Q14).
q         Fixed-point position for output data, 0..15 (aunwrap)
N         Size of input/output arrays
Output:
y[N]      Unwrapped phase angles

Restrictions:
x,y       Aligned on 32-byte boundary
x,y       Must not overlap
N         Multiple of 16 (aunwrap) or 8 (aunwrapf)
-------------------------------------------------------------------------*/

void aunwrap ( int16_t* restrict y, 
         const int16_t * restrict x, 
         int16_t tol,
         int q, int N )
{
    const xb_vecNx16 * restrict X;
    const xb_vecNx16 * restrict X1;
    xb_vecNx16 * restrict Y;

    xb_vecNx16 t, x2, x1, x0, y0, t0, dx, b0, _1 = BBE_MOVVINT16(1);
    xb_vecNx40 DX0, DX, B0;
    xb_vecNx16 vtol = BBE_MOVVA16(tol), _0 = 0;
    vboolN     lte;
#if USE_ROUNDING
    vsaN       shft = BBE_MOVVSA32(q);
#else
    vsaN       shft = BBE_MOVVSA32(15 - q);
#endif
    vsaN       sh1 = BBE_MOVVSA32(1), sh16 = BBE_MOVVSA32(16);
    int n;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % BBE_SIMD_WIDTH == 0);
    NASSERT(q >= 0 && q <= 15);

    X = (const xb_vecNx16*)x;
    X1 = X;
    Y = (xb_vecNx16*)y;

    x0 = 0;
    b0 = 0;

    __Pragma("loop_count min=1");
    for (n = 0; n<N >> LOG2_BBE_SIMD_WIDTH; n++)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x1, X, 0);
        x2 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_LEFT_1);
        // dx0=(int32_t)x[n]-x[n-1];
        DX0 = BBE_UNPKSNX16(x1);
        BBE_MULSNX16(DX0, x2, _1);
        x2 = BBE_PACKLNX40(DX0);
        DX = BBE_UNPKSNX16(x2);
        // dx =((int16_t)dx0)-dx0;
        DX = BBE_SUBNX40(DX, DX0);
        DX0 = BBE_ABSNX40(DX0);
        //if (dx0<tol) dx=0; 
        t = BBE_PACKVNX40(DX0, sh1);
        lte = BBE_LENX16(vtol, t);
        t = BBE_PACKVNX40(DX, sh16);
        dx = BBE_MOVNX16T(t, _0, lte);

        t0 = BBE_SELNX16I(dx, _0, BBE_SELI_ROTATE_LEFT_1);
        dx = BBE_ADDNX16(dx, t0);
        t0 = BBE_SELNX16I(dx, _0, BBE_SELI_ROTATE_LEFT_2);
        dx = BBE_ADDNX16(dx, t0);
        t0 = BBE_SELNX16I(dx, _0, BBE_SELI_ROTATE_LEFT_4);
        dx = BBE_ADDNX16(dx, t0);
        t0 = BBE_SELNX16I(dx, _0, BBE_SELI_ROTATE_LEFT_8);
        dx = BBE_ADDNX16(dx, t0);

        dx = BBE_ADDNX16(dx, b0);
        b0 = BBE_REPNX16(dx, (BBE_SIMD_WIDTH - 1));
        B0 = BBE_UNPKSNX16(dx);
        B0 = BBE_SLLINX40(B0, 16);
        BBE_LVNX16_XP(x1, X1, 2 * BBE_SIMD_WIDTH);
        BBE_MULANX16(B0, x1, _1);
#if USE_ROUNDING
        y0 = BBE_PACKNVNX40(B0, shft);
#else
        B0 = BBE_RNDADJNX40(B0, shft);
        y0 = BBE_PACKVNX40(B0, shft);
#endif
        BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(x0, X, 2 * BBE_SIMD_WIDTH); //x0 = x1;
    }
} /* aunwrap() */

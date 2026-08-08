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
    Real Vector Product
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Dot Product of Real Vectors

Description: These routines take two real vectors and calculate their dot 
product.

Representation:
rvdot   Signed fixed-point format
        Input vectors x and y are 16-bit signed data of arbitrary formats
        Qx and Qy. Dot product is computed in a 40-bit accumulator, which
        is then rounded, shifted to the right by rsh bit positions and
        saturated to form a 32-bit result with Qx+Qy-rsh fractional bits.
rvdotf  IEEE-754 Std. single precision floating-point format for
        input vectors and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (rvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (rvdot) or 8 (rvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/

int32_t rvdot ( const int16_t * restrict x,
                const int16_t * restrict y,
                int rsh, 
                int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pY = (const xb_vecNx16 *)y;
    int k;
    xb_vecNx16 x0, y0;
    xb_vecNx40 A0;
    xb_int40   A;
    //vsaN VsaShift;

    if (N <= 0) return(0);
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N%BBE_SIMD_WIDTH == 0);
    NASSERT(rsh >= 0);

    //VsaShift = BBE_MOVVSA32(rsh - 3);

    A0 = 0;
    for (k = 0; k<(N >> LOG2_BBE_SIMD_WIDTH); k++)
    {
        BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(y0, pY, 2 * BBE_SIMD_WIDTH);

        BBE_MULANX16(A0, x0, y0);
    }

    //A0 = BBE_RNDADJNX40(A0, VsaShift); // rounding code 
    A = BBE_RADDNX40(A0);
    k = (A >> rsh);
    return k;

} /* rvdot() */

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
    Sum of Squares of a Vector
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
Sum of Squares of a Vector

Description: These routines compute the power of a vector.

Representation:
vpower   Signed fixed-point format
         Input vector elements are 16-bit signed data of arbitrary format
         Qx. Sum of squared values is computed in a 40-bit accumulator,
         which is then rounded, shifted to the right by rsh bit positions and
         saturated to form a 32-bit result with 2*Qx-rsh fractional bits.
vpowerf  IEEE-754 Std. single precision floating-point format for the
         input vector and the result

Parameters:
Input:
x[N]     Input vector
rsh      Right shift amount (vpower)
N        Length of input vector 
Returned Value:
Sum of squares over the input vector

Restrictions:
x        Aligned on 32-byte boundary
N        Multiple of 16 (vpower) or 8 (vpowerf)
rsh>=0   Right shift amount must be non-negative
-------------------------------------------------------------------------*/

int32_t vpower ( const int16_t * restrict x,
                 int rsh,
                 int N )
{
    const xb_vecNx16 * restrict pX0 = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pX1 = (const xb_vecNx16 *)(x + BBE_SIMD_WIDTH);
    int k;
    xb_vecNx16 x0, x1;
    xb_vecNx40 A0;
    xb_int40   A;
    vsaN VsaShift;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N%BBE_SIMD_WIDTH == 0);

    VsaShift = BBE_MOVVSA32(rsh - 3);
    A0 = 0;
    for (k = 0; k<(N >> (LOG2_BBE_SIMD_WIDTH + 1)); k++)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x0, pX0, 4 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(x1, pX1, 4 * BBE_SIMD_WIDTH);
        BBE_MAGIANX16C(A0, x0, x1);
    }
    if (N&(BBE_SIMD_WIDTH * 2 - 1))
    {
        BBE_LVNX16_IP(x0, pX0, 2 * BBE_SIMD_WIDTH);
        BBE_MULANX16(A0, x0, x0);
    }
    A0 = BBE_RNDADJNX40(A0, VsaShift); // rounding code 
    A = BBE_RADDNX40(A0);

    k = (A >> rsh);
    return k;
} /* vpower() */

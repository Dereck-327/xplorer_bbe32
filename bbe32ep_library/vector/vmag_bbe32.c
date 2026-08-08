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
    Square Root of Sum of Squares
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
Square Root of Sum of Squares 

Description: These routines compute the magnitude of a vector

Representation:
vmag    Signed fixed-point format
        Input vector elements are 16-bit signed data of arbitrary format
        Qx. Sum of squared values is computed in a 40-bit accumulator,
        which is then shifted to the right by rsh bit positions. Square
        root approximation of the shifted sum is formatted as a 16-bit
        non-negative value with Qx-rsh/2 fractional bits.
vmagf   IEEE-754 Std. single precision floating-point format for
        input vector and magnitude estimation result

Parameters:
Input:
x[N]    Input data
rsh     Right shift amount (vmag)
N       Length of vector 
Returned Value:
Square root of the sum of squares over the input vector

Restrictions:
x       Aligned on 32-byte boundary
N       Multiple of 16 (vmag) or 8 (vmagf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/

#if !(HAVE_NSAENX40 && HAVE_RECIP && 1)

DISCARD_FUN(int16_t, vmag, (const int16_t * restrict x,
    int rsh,
    int N))

#else

int16_t vmag(const int16_t * restrict x,
               int rsh,
               int N )
{
    xb_vecNx40 w;
    xb_int40 iw;
    const xb_vecNx16* restrict X0 = (const xb_vecNx16*)x;
    const xb_vecNx16* restrict X1 = (const xb_vecNx16*)(x + BBE_SIMD_WIDTH);
    xb_vecNx16 x0, x1;
    int k;
    vsaN vrsh = BBE_MOVVSA32(1 - rsh);
    xb_vecNx16   b, c, _19 = BBE_MOVVINT16(19);
    vsaN nsa;
    vsaN VsaShift;

    if (N <= 0) return(0);
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N%BBE_SIMD_WIDTH == 0);
    NASSERT(rsh >= 0);

    VsaShift = BBE_MOVVSA32(rsh - 3);
    w = 0;
    for (k = 0; k<(N >> (LOG2_BBE_SIMD_WIDTH + 1)); k++)
    {
        // Here BBE_LVNX16_XP provide better schedule than BBE_LVNX16_IP!
        BBE_LVNX16_XP(x0, X0, 4 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(x1, X1, 4 * BBE_SIMD_WIDTH);
        BBE_MAGIANX16C(w, x0, x1);
    }
    if (N&(BBE_SIMD_WIDTH * 2 - 1))
    {
        BBE_LVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
        BBE_MULANX16(w, x0, x0);
    }
    w = BBE_RNDADJNX40(w, VsaShift); // rounding code 
    iw = BBE_RADDNX40(w);
    w = BBE_MOVNX40_FROM40(iw);
    w = BBE_SLSNX40(w, vrsh);
    // take square root
    nsa = BBE_NSAENX40(w);
    w = BBE_SLLNX40(w, nsa);
    BBE_RSQRTLUNX40_0(w, b, c, w);
    BBE_MULUUSNX16(w, c, b);
    BBE_RECIPLUNX40_0(w, c, b, w);
    BBE_MULUSANX16(w, b, c);
    b = BBE_MOVVVS(nsa);
    b = BBE_SRAINX16(b, 1);
    b = BBE_ADDNX16(b, _19);
    nsa = BBE_MOVVSV(b, 0);
    b = BBE_PACKVNX40(w, nsa);
    return BBE_MOVAV16(b);
} /* vmag() */

#endif

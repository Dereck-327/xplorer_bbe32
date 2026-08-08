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
    NatureDSP_Baseband library. FIR filters and Related Functions
    Block Real FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "bkfirf_common.h"

/*-------------------------------------------------------------------------
Block Real FIR Filter

Computes a real FIR filter (direct-form) using IR stored in vector h. The
real data input is stored in vector x. The filter output result is stored
in vector y. The filter calculates N output samples using M coefficients
and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by bkfir[f]_algDelay(M).

Representation:
bkfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
bkfirf  IEEE-754 Std. single precision floating-point format for filter 
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[M]    Filter coefficients; h[0] is to be multiplied by the newest sample
N       Length of sample block
M       Length of filter
x[N]    Input samples
Output:
y[N]    Output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (bkfir) or 8 (bkfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, bkfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,bkfirf_process_16,(float32_t* restrict delay, float32_t* restrict y, const float32_t * restrict x, const float32_t* restrict h, int M, int N))
#else
void bkfirf_process_16(float32_t* restrict delay, float32_t* restrict y, const float32_t * restrict x, const float32_t* restrict h, int M, int N)
{
    int n;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)delay;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)h;

    xb_vecN_2xf32 d0, d1, d2, c0;
    xb_vecN_2xf32 h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5, x6, x7;
    xb_vecN_2xf32 A0, A1, A2, A3;

    NASSERT(N > 0 && !(N & 7));
    NASSERT(M == 16);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delay);
    NASSERT_ALIGN32(h);

    BBE_LVN_2XF32_XP(d1, pD, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_XP(d2, pD, -2 * BBE_SIMD_WIDTH);

    for (n = 0; n < (N / (BBE_SIMD_WIDTH / 2)); n++)
    {
        d0 = d1;
        d1 = d2;
        BBE_LVN_2XF32_IP(d2, pX, 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_XP(c0, pH, 2 * BBE_SIMD_WIDTH);
        h0 = BBE_REPN_2XF32(c0, 0);
        h1 = BBE_REPN_2XF32(c0, 1);
        h2 = BBE_REPN_2XF32(c0, 2);
        h3 = BBE_REPN_2XF32(c0, 3);
        h4 = BBE_REPN_2XF32(c0, 4);
        h5 = BBE_REPN_2XF32(c0, 5);
        h6 = BBE_REPN_2XF32(c0, 6);
        h7 = BBE_REPN_2XF32(c0, 7);

        x0 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_2);
        x1 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_4);
        x2 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_6);
        x3 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
        x4 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_10);
        x5 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_12);
        x6 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_14);
        x7 = d1;

        A0 = BBE_MULN_2XF32(h0, x0);
        A1 = BBE_MULN_2XF32(h1, x1);
        A2 = BBE_MULN_2XF32(h2, x2);
        A3 = BBE_MULN_2XF32(h3, x3);
        BBE_MULAN_2XF32(A0, h4, x4);
        BBE_MULAN_2XF32(A1, h5, x5);
        BBE_MULAN_2XF32(A2, h6, x6);
        BBE_MULAN_2XF32(A3, h7, x7);

        BBE_LVN_2XF32_XP(c0, pH, -2 * BBE_SIMD_WIDTH);
        h0 = BBE_REPN_2XF32(c0, 0);
        h1 = BBE_REPN_2XF32(c0, 1);
        h2 = BBE_REPN_2XF32(c0, 2);
        h3 = BBE_REPN_2XF32(c0, 3);
        h4 = BBE_REPN_2XF32(c0, 4);
        h5 = BBE_REPN_2XF32(c0, 5);
        h6 = BBE_REPN_2XF32(c0, 6);
        h7 = BBE_REPN_2XF32(c0, 7);

        x0 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_2);
        x1 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_4);
        x2 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_6);
        x3 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_8);
        x4 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_10);
        x5 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_12);
        x6 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_14);
        x7 = d2;

        BBE_MULAN_2XF32(A0, h0, x0);
        BBE_MULAN_2XF32(A1, h1, x1);
        BBE_MULAN_2XF32(A2, h2, x2);
        BBE_MULAN_2XF32(A3, h3, x3);
        BBE_MULAN_2XF32(A0, h4, x4);
        BBE_MULAN_2XF32(A1, h5, x5);
        BBE_MULAN_2XF32(A2, h6, x6);
        BBE_MULAN_2XF32(A3, h7, x7);
        A0 = BBE_ADDN_2XF32(A0, A1);
        A2 = BBE_ADDN_2XF32(A2, A3);
        A0 = BBE_ADDN_2XF32(A0, A2);

        BBE_SVN_2XF32_IP(A0, pY, 2 * BBE_SIMD_WIDTH);
    }
    BBE_SVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(d2, pD, -2 * BBE_SIMD_WIDTH);
}
#endif

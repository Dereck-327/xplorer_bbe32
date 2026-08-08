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
    Interpolating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

#include "firinterpf_common.h"
#if HAVE_VFPU
/*-------------------------------------------------------------------------
Interpolating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with interpolation using real
IR stored in vector h. The complex data input is stored in vector x. The
filter output result is stored in vector y. The filter calculates N*D complex
output samples using M*D coefficients and requires last N+M-1 samples in the
delay line.

Representation:
firinterp   16-bit signed fixed-point format
            Filter coefficients are Q15
            Number of fractional bits for input/output samples is user-difined
firinterpf  IEEE-754 Std. single precision floating-point format for filter 
            coefficients and input/output samples

Parameters:
Input:
D           Interpolation ratio 
N           Length of input sample block
M           Length of subfilter. Total length of filter is M*D
h[M*D]      Filter coefficients; h[0] is to be multiplied by the newest 
            sample,Q15
x[N]        Input complex samples
Output:
y[N*D]      Output complex samples

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
N           Multiple of 8 (firinterp) or 4 (firinterpf)
M           2,4,8 or a positive multiple of 16 for D=2,3,4,6,12; or 
            a positive multiple of 8 for other D
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
subfilter lengths M=2,4,8,16 and 32 and interpolation factors D=2,3 and 4,
in any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firinterp[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

void firinterpf_6d_intlv(complex_float * y, int N)
{
    int n;
    const xb_vecN_2xf32 * restrict pYr = (const xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pYw = (      xb_vecN_2xf32 *)y;
    xb_vecN_2xf32 y0, y1, y2, y3, y4, y5;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5;
    xb_vecN_2xf32 z0, z1, z2, z3, z4, z5;

    __Pragma("loop_count min=1");
    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        BBE_LVN_2XF32_IP(y0, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y1, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y2, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y3, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y4, pYr, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y5, pYr, 2 * BBE_SIMD_WIDTH);

        BBE_DSELN_2XF32I(x1, x0, y1, y0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XF32I(x3, x2, y3, y2, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XF32I(x5, x4, y5, y4, BBE_DSELI_INTERLEAVE_4);
        z0 = BBE_SELN_2XF32I(x2, x0, BBE_SELI_EXTRACT_LO_HALVES);
        z2 = BBE_SELN_2XF32I(x4, x2, BBE_SELI_EXTRACT_HI_HALVES);
        z3 = BBE_SELN_2XF32I(x3, x1, BBE_SELI_EXTRACT_LO_HALVES);
        z5 = BBE_SELN_2XF32I(x5, x3, BBE_SELI_EXTRACT_HI_HALVES);

#if 0
        z1 = BBE_SELN_2XF32I(x4, x0, BBE_SELI_ROTATE_RIGHT_8);
        z4 = BBE_SELN_2XF32I(x5, x1, BBE_SELI_ROTATE_RIGHT_8);
        z1 = BBE_SHFLN_2XF32I(z1, BBE_SHFLI_SWAP_8);
        z4 = BBE_SHFLN_2XF32I(z4, BBE_SHFLI_SWAP_8);
#else
        z1 = BBE_SELN_2XF32I(y0, y4, BBE_SELI_INTERLEAVE_4_EVENODD);
        z4 = BBE_SELN_2XF32I(y1, y5, BBE_SELI_INTERLEAVE_4_EVENODD);
        BBE_DSELN_2XF32I(z4, z1, z4, z1, BBE_DSELI_INTERLEAVE_4);
#endif

        BBE_SVN_2XF32_IP(z0, pYw, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z1, pYw, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z2, pYw, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z3, pYw, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z4, pYw, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z5, pYw, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

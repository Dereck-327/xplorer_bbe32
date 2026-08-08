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
#if !HAVE_VFPU
DISCARD_FUN (void,firinterpf_proc_D6_M8,( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D ))
#else
void firinterpf_proc_D6_M8( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D )
{
    int n, d;
    const xb_vecN_2xf32 * restrict pX;
          xb_vecN_2xf32 * restrict pY;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)delayLine;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)coef;

    xb_vecN_2xf32 c0, h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 xx, x0, x1, x2, x3, x4, x5, x6, x7;
    xb_vecN_2xf32 d0, d1;
    xb_vecN_2xf32 y0, y1, y2, y3;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 6 && M == 8);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    for (d = 0; d < D; d++)
    {
        pX = (const xb_vecN_2xf32 *)x;
        pY = (      xb_vecN_2xf32 *)y + d;
        pH = (const xb_vecN_2xf32 *)coef + d;
        BBE_LVN_2XF32_XP(c0, pH, 0);
        BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(d1, pD, -2 * BBE_SIMD_WIDTH);
        x0 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_4);
        x1 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
        x2 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_12);

        __Pragma("loop_count min=1");
        for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
        {
            BBE_LVN_2XF32_IP(xx, pX, 2 * BBE_SIMD_WIDTH);
            x3 = d1;
            x4 = BBE_SELN_2XF32I(xx, d1, BBE_SELI_ROTATE_RIGHT_4);
            x5 = BBE_SELN_2XF32I(xx, d1, BBE_SELI_ROTATE_RIGHT_8);
            x6 = BBE_SELN_2XF32I(xx, d1, BBE_SELI_ROTATE_RIGHT_12);
            x7 = xx;

            //BBE_LVN_2XF32_XP(c0, pH, 0);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            y0 = BBE_MULN_2XF32(x0, h0);
            y1 = BBE_MULN_2XF32(x1, h1);
            y2 = BBE_MULN_2XF32(x2, h2);
            y3 = BBE_MULN_2XF32(x3, h3);
            BBE_MULAN_2XF32(y0, x4, h4);
            BBE_MULAN_2XF32(y1, x5, h5);
            BBE_MULAN_2XF32(y2, x6, h6);
            BBE_MULAN_2XF32(y3, x7, h7);
            y0 = BBE_ADDN_2XF32(y0, y1);
            y2 = BBE_ADDN_2XF32(y2, y3);
            y0 = BBE_ADDN_2XF32(y0, y2);

            BBE_SVN_2XF32_XP(y0, pY, D * 2 * BBE_SIMD_WIDTH);
            x0 = x4;
            x1 = x5;
            x2 = x6;
            d1 = xx;
        }
    }
    BBE_SVN_2XF32_I(x3, pD, 0);
    BBE_SVN_2XF32_I(d1, pD, 2 * BBE_SIMD_WIDTH);
    firinterpf_6d_intlv(y, N);
}// firinterpf_proc_D6_M8
#endif

const tFirFxdxns interpf_6d_8m = { &firinterpf_d_6_12, firinterpf_proc_D6_M8 };

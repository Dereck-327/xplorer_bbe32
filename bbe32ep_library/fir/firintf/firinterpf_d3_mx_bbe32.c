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
DISCARD_FUN (void,firinterpf_proc_D3_MX,( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D ))
#else
void firinterpf_proc_D3_MX( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D )
{
    int n, m;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pH;

    xb_vecN_2xf32 c0, h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 x0, x1, x2, x3;
    xb_vecN_2xf32 d0, d1, d2;
    xb_vecN_2xf32 y0, y1, y2;
    xb_vecN_2xf32 z1, z2;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 3 && !(M & 7));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        y0 = BBE_ZERON_2XF32();
        y1 = BBE_ZERON_2XF32();
        y2 = BBE_ZERON_2XF32();
        pD = (      xb_vecN_2xf32 *)delayLine;
        pH = (const xb_vecN_2xf32 *)coef;

        BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X(x0, pD, 4 * 2 * M);

        BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);

        for (m = 0; m < (M / (BBE_SIMD_WIDTH / 2)); m++)
        {
            BBE_LVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_12);
            x3 = d1;

            BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            BBE_MULAN_2XF32(y0, x0, h0);
            BBE_MULAN_2XF32(y1, x0, h1);
            BBE_MULAN_2XF32(y2, x0, h2);
            BBE_MULAN_2XF32(y0, x1, h3);
            BBE_MULAN_2XF32(y1, x1, h4);
            BBE_MULAN_2XF32(y2, x1, h5);
            BBE_MULAN_2XF32(y0, x2, h6);
            BBE_MULAN_2XF32(y1, x2, h7);

            BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            BBE_MULAN_2XF32(y2, x2, h0);
            BBE_MULAN_2XF32(y0, x3, h1);
            BBE_MULAN_2XF32(y1, x3, h2);
            BBE_MULAN_2XF32(y2, x3, h3);

            BBE_LVN_2XF32_XP(d2, pD, -2 * 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_12);
            x3 = d2;

            BBE_MULAN_2XF32(y0, x0, h4);
            BBE_MULAN_2XF32(y1, x0, h5);
            BBE_MULAN_2XF32(y2, x0, h6);
            BBE_MULAN_2XF32(y0, x1, h7);

            BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            BBE_MULAN_2XF32(y1, x1, h0);
            BBE_MULAN_2XF32(y2, x1, h1);
            BBE_MULAN_2XF32(y0, x2, h2);
            BBE_MULAN_2XF32(y1, x2, h3);
            BBE_MULAN_2XF32(y2, x2, h4);
            BBE_MULAN_2XF32(y0, x3, h5);
            BBE_MULAN_2XF32(y1, x3, h6);
            BBE_MULAN_2XF32(y2, x3, h7);

            BBE_SVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(d2, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
        }
        z1 = BBE_SELN_2XF32I(y1, y2, BBE_SELI_INTERLEAVE_2_EVENODD);
        z1 = BBE_SHFLN_2XF32I(z1, BBE_SHFLI_SWAP_2);
        z2 = BBE_SELN_2XF32I(y2, y1, BBE_SELI_INTERLEAVE_2_EVENODD);
        BBE_DSELN_2XF32I(x1, x0, z1, y0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XF32I(y1, y0, z2, x0, BBE_DSELI_INTERLEAVE_C3_STEP_0);
        BBE_DSELN_2XF32I_H(y1, y2, z2, x1, BBE_DSELI_INTERLEAVE_C3_STEP_1);

        BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y1, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y2, pY, 2 * BBE_SIMD_WIDTH);
    }
}// firinterpf_proc_D3_MX
#endif

const tFirFxdxns interpf_3d_xm = { &firinterpf_d_2_3_4_mx, firinterpf_proc_D3_MX };

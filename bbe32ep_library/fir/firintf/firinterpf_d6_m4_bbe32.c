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
DISCARD_FUN (void,firinterpf_proc_D6_M4,( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D ))
#else
void firinterpf_proc_D6_M4( complex_float * restrict y,
                             const complex_float * restrict x,
                             const float32_t     * restrict coef,
                                   float32_t     * restrict delayLine,
                             int M, int N, int D )
{
#if 0
    int n;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
    xb_vecN_2xf32 * restrict pY = (xb_vecN_2xf32 *)y;
    xb_vecN_2xf32 * restrict pD = (xb_vecN_2xf32 *)delayLine;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)coef;

    xb_vecN_2xf32 c0, c1, c2, h0, h1, h2, h3, h4, h5;
    xb_vecN_2xf32 xx, x0, x1, x2, x3;
    xb_vecN_2xf32 d0;
    xb_vecN_2xf32 y0, y1, y2, y3, y4, y5;
    xb_vecN_2xf32 z0, z1, z2, z3, z4, z5, x4, x5;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 6 && M == 4);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    d0 = BBE_LVN_2XF32_I(pD, 0);

    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        y0 = y1 = y2 = y3 = y4 = y5 = 0;
        BBE_LVN_2XF32_IP(xx, pX, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_4);
        x1 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_8);
        x2 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_12);
        x3 = xx;

        BBE_LVN_2XF32_XP(c0, pH, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(c1, pH, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(c2, pH, -2 * 2 * BBE_SIMD_WIDTH);

        h0 = BBE_REPN_2XF32(c0, 0);
        h1 = BBE_REPN_2XF32(c0, 1);
        h2 = BBE_REPN_2XF32(c0, 2);
        h3 = BBE_REPN_2XF32(c0, 3);
        h4 = BBE_REPN_2XF32(c0, 4);
        h5 = BBE_REPN_2XF32(c0, 5);
        BBE_MULAN_2XF32(y0, x0, h0);
        BBE_MULAN_2XF32(y1, x0, h1);
        BBE_MULAN_2XF32(y2, x0, h2);
        BBE_MULAN_2XF32(y3, x0, h3);
        BBE_MULAN_2XF32(y4, x0, h4);
        BBE_MULAN_2XF32(y5, x0, h5);

        h0 = BBE_REPN_2XF32(c0, 6);
        h1 = BBE_REPN_2XF32(c0, 7);
        h2 = BBE_REPN_2XF32(c1, 0);
        h3 = BBE_REPN_2XF32(c1, 1);
        h4 = BBE_REPN_2XF32(c1, 2);
        h5 = BBE_REPN_2XF32(c1, 3);
        BBE_MULAN_2XF32(y0, x1, h0);
        BBE_MULAN_2XF32(y1, x1, h1);
        BBE_MULAN_2XF32(y2, x1, h2);
        BBE_MULAN_2XF32(y3, x1, h3);
        BBE_MULAN_2XF32(y4, x1, h4);
        BBE_MULAN_2XF32(y5, x1, h5);

        h0 = BBE_REPN_2XF32(c1, 4);
        h1 = BBE_REPN_2XF32(c1, 5);
        h2 = BBE_REPN_2XF32(c1, 6);
        h3 = BBE_REPN_2XF32(c1, 7);
        h4 = BBE_REPN_2XF32(c2, 0);
        h5 = BBE_REPN_2XF32(c2, 1);
        BBE_MULAN_2XF32(y0, x2, h0);
        BBE_MULAN_2XF32(y1, x2, h1);
        BBE_MULAN_2XF32(y2, x2, h2);
        BBE_MULAN_2XF32(y3, x2, h3);
        BBE_MULAN_2XF32(y4, x2, h4);
        BBE_MULAN_2XF32(y5, x2, h5);

        h0 = BBE_REPN_2XF32(c2, 2);
        h1 = BBE_REPN_2XF32(c2, 3);
        h2 = BBE_REPN_2XF32(c2, 4);
        h3 = BBE_REPN_2XF32(c2, 5);
        h4 = BBE_REPN_2XF32(c2, 6);
        h5 = BBE_REPN_2XF32(c2, 7);
        BBE_MULAN_2XF32(y0, x3, h0);
        BBE_MULAN_2XF32(y1, x3, h1);
        BBE_MULAN_2XF32(y2, x3, h2);
        BBE_MULAN_2XF32(y3, x3, h3);
        BBE_MULAN_2XF32(y4, x3, h4);
        BBE_MULAN_2XF32(y5, x3, h5);



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

        BBE_SVN_2XF32_IP(z0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z1, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z2, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z3, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z4, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(z5, pY, 2 * BBE_SIMD_WIDTH);
        d0 = xx;
    }
    BBE_SVN_2XF32_I(d0, pD, 0);  
#endif

//34cycles
#if 0
    int n;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)delayLine;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)coef;

    xb_vecN_2xf32 c0, c1, c2, h0, h1, h2, h3, h4, h5;
    xb_vecN_2xf32 xx, x0, x1, x2, x3, x01, x12, x23;
    xb_vecN_2xf32 d0, d1, d2, d01, d12, d23;
    xb_vecN_2xf32 y0, y1, y2, y3, y4, y5;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 6 && M == 4);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    d0 = BBE_LVN_2XF32_I(pD, 0);

    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        BBE_LVN_2XF32_IP(xx, pX, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_SHFLN_2XF32I(xx, BBE_SHFLI_REP_0X4);
        x1 = BBE_SHFLN_2XF32I(xx, BBE_SHFLI_REP_1X4);
        x2 = BBE_SHFLN_2XF32I(xx, BBE_SHFLI_REP_2X4);
        x3 = BBE_SHFLN_2XF32I(xx, BBE_SHFLI_REP_3X4);
        x01 = BBE_SELN_2XF32I(x1, x0, BBE_SELI_ROTATE_RIGHT_8);
        x12 = BBE_SELN_2XF32I(x2, x1, BBE_SELI_ROTATE_RIGHT_8);
        x23 = BBE_SELN_2XF32I(x3, x2, BBE_SELI_ROTATE_RIGHT_8);
        d2 = BBE_SHFLN_2XF32I(d0, BBE_SHFLI_REP_3X4);
        d1 = BBE_SHFLN_2XF32I(d0, BBE_SHFLI_REP_2X4);
        d0 = BBE_SHFLN_2XF32I(d0, BBE_SHFLI_REP_1X4);
        d01 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
        d12 = BBE_SELN_2XF32I(d2, d1, BBE_SELI_ROTATE_RIGHT_8);
        d23 = BBE_SELN_2XF32I(x0, d2, BBE_SELI_ROTATE_RIGHT_8);


        BBE_LVN_2XF32_XP(c0, pH, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(c1, pH, 0);
        h0 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_DOUBLE_2_LO);
        h1 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_DOUBLE_2_HI);
        h5 = BBE_SHFLN_2XF32I(c1, BBE_SHFLI_DOUBLE_2_LO);
        h2 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_ROTATE_RIGHT_8);
        h3 = BBE_SELN_2XF32I(h5, h1, BBE_SELI_ROTATE_RIGHT_8);
        h4 = BBE_SELN_2XF32I(h1, h5, BBE_SELI_EXTRACT_HI_HALVES);
        h1 = BBE_SELN_2XF32I(h0, h1, BBE_SELI_EXTRACT_LO_HALVES);
        y0 = BBE_MULN_2XF32(d0 , h0);
        y1 = BBE_MULN_2XF32(d01, h1);
        y2 = BBE_MULN_2XF32(d1 , h2);
        y3 = BBE_MULN_2XF32(d2 , h0);
        y4 = BBE_MULN_2XF32(d23, h1);
        y5 = BBE_MULN_2XF32(x0 , h2);
        BBE_MULAN_2XF32(y0, d1 , h3);
        BBE_MULAN_2XF32(y1, d12, h4);
        BBE_MULAN_2XF32(y2, d2 , h5);
        BBE_MULAN_2XF32(y3, x0 , h3);
        BBE_MULAN_2XF32(y4, x01, h4);
        BBE_MULAN_2XF32(y5, x1 , h5);

        BBE_LVN_2XF32_XP(c1, pH, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(c2, pH, -2 * 2 * BBE_SIMD_WIDTH);
        h0 = BBE_SHFLN_2XF32I(c1, BBE_SHFLI_DOUBLE_2_HI);
        h1 = BBE_SHFLN_2XF32I(c2, BBE_SHFLI_DOUBLE_2_LO);
        h5 = BBE_SHFLN_2XF32I(c2, BBE_SHFLI_DOUBLE_2_HI);
        h2 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_ROTATE_RIGHT_8);
        h3 = BBE_SELN_2XF32I(h5, h1, BBE_SELI_ROTATE_RIGHT_8);
        h4 = BBE_SELN_2XF32I(h1, h5, BBE_SELI_EXTRACT_HI_HALVES);
        h1 = BBE_SELN_2XF32I(h0, h1, BBE_SELI_EXTRACT_LO_HALVES);
        BBE_MULAN_2XF32(y0, d2 , h0);
        BBE_MULAN_2XF32(y1, d23, h1);
        BBE_MULAN_2XF32(y2, x0 , h2);
        BBE_MULAN_2XF32(y3, x1 , h0);
        BBE_MULAN_2XF32(y4, x12, h1);
        BBE_MULAN_2XF32(y5, x2 , h2);
        BBE_MULAN_2XF32(y0, x0 , h3);
        BBE_MULAN_2XF32(y1, x01, h4);
        BBE_MULAN_2XF32(y2, x1 , h5);
        BBE_MULAN_2XF32(y3, x2 , h3);
        BBE_MULAN_2XF32(y4, x23, h4);
        BBE_MULAN_2XF32(y5, x3 , h5);

        BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y1, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y2, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y3, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y4, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_IP(y5, pY, 2 * BBE_SIMD_WIDTH);
        d0 = xx;
    }
    BBE_SVN_2XF32_I(d0, pD, 0);
#endif

#if 1
    int n, d;
    const xb_vecN_2xf32 * restrict pX;
          xb_vecN_2xf32 * restrict pY;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)delayLine;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)coef;

    xb_vecN_2xf32 c0, h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 xx, x0, x1, x2, x3;
    xb_vecN_2xf32 d0;
    xb_vecN_2xf32 y0, y1;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 6 && M == 4);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    for (d = 0; d < D / 2; d++)
    {
        pX = (const xb_vecN_2xf32 *)x;
        pY = (      xb_vecN_2xf32 *)y + 2 * d;
        BBE_LVN_2XF32_XP(c0, pH, 2 * BBE_SIMD_WIDTH);
        d0 = BBE_LVN_2XF32_I(pD, 0);

        __Pragma("loop_count min=1");
        for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
        {
            BBE_LVN_2XF32_IP(xx, pX, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(xx, d0, BBE_SELI_ROTATE_RIGHT_12);
            x3 = xx;

            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            y0 = BBE_MULN_2XF32(x0, h0);
            y1 = BBE_MULN_2XF32(x0, h4);
            BBE_MULAN_2XF32(y0, x1, h1);
            BBE_MULAN_2XF32(y1, x1, h5);
            BBE_MULAN_2XF32(y0, x2, h2);
            BBE_MULAN_2XF32(y1, x2, h6);
            BBE_MULAN_2XF32(y0, x3, h3);
            BBE_MULAN_2XF32(y1, x3, h7);

            BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_XP(y1, pY, (D - 1) * 2 * BBE_SIMD_WIDTH);
            d0 = xx;
        }
    }
    BBE_SVN_2XF32_I(d0, pD, 0);
    firinterpf_6d_intlv(y, N);
#endif
}// firinterpf_proc_D6_M4
#endif

const tFirFxdxns interpf_6d_4m = { &firinterpf_d_6_12, firinterpf_proc_D6_M4 };

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
    Decimating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "firdecf_common.h"

/*-------------------------------------------------------------------------
Decimating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with decimation using real IR 
stored in vector h. The complex data input is stored in vector x. The filter
output result is stored in vector y. The filter calculates N output samples
using M coefficients and requires last D*N+M-1 samples in the delay line.

NOTE:
To avoid aliasing, the IR should be synthesized in such a way that filter pass
band is limited by input sample frequency divided by 2*D.

Representation:
firdec   16-bit signed fixed-point format
         Filter coefficients are Q15
         Number of fractional bits for input/output samples is user-difined
firdecf  IEEE-754 Std. single precision floating-point format for filter 
         coefficients and input/output samples

Parameters:
Input:
D        Decimation factor
N        Length of output sample block
M        Length of filter
h[M]     Filter coefficients; h[0] is to be multiplied by the newest 
         sample
x[N*D]   Input complex samples
Output:
y[N]     Output complex samples

Restrictions:
x,y      Must not overlap
x,y      Aligned on 32-byte boundary
N        Multiple of 8 (firdec) or 4 (firdecf)
M        2,4,8 or a positive multiple of 16 for D=2,3,4; or 
         a positive multiple of 16 for D>4
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
filter lengths M=2,4,8,16 and 32 and decimation factors D=2,3 and 4, in
any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firdec[f]_init returns NULL handle.
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,firdecf_proc_D3_MX,( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D ))
#else
void firdecf_proc_D3_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D )
{
    int n, m;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pH_01;
    const xb_vecN_2xf32 * restrict pH_2;

    xb_vecN_2xf32 c0, h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5, x6, x7;
    xb_vecN_2xf32 s0, s1, s2, s3;
    xb_vecN_2xf32 d0, d1, d2;
    xb_vecN_2xf32 y0;
    valign vH;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 3 && !(M % 12));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        y0 = BBE_ZERON_2XF32();
        pD = (xb_vecN_2xf32 *)delayLine;
        pH_01 = (const xb_vecN_2xf32 *)coef;
        pH_2  = (const xb_vecN_2xf32 *)(coef + 2 * M / 3);
        vH = BBE_LAN_2XF32_PP(pH_2);

        BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);//r0,r1,r2,r3
        BBE_LVN_2XF32_IP(x1, pX, 2 * BBE_SIMD_WIDTH);//s0,s1,s2,s3
        BBE_LVN_2XF32_IP(x2, pX, 2 * BBE_SIMD_WIDTH);//t0,t1,t2,t3
        s0 = BBE_SELN_2XF32I(x0, x0, BBE_SELI_ROTATE_RIGHT_4);//r1,r2,r3,r0
        s1 = BBE_SELN_2XF32I(x1, s0, BBE_SELI_EXTRACT_LO_HALVES);//r1,r2,s0,s1
        s2 = BBE_SELN_2XF32I(x2, x1, BBE_SELI_ROTATE_RIGHT_4);//s1,s2,s3,t0
        s3 = BBE_SELN_2XF32I(s2, x2, BBE_SELI_INTERLEAVE_4_LO);//t0,s1,t1,s2
        x0 = BBE_SELN_2XF32I(s3, s0, BBE_SELI_EXTRACT_HI_HALVES);//r3,r0,t1,s2
        s0 = BBE_SHFLN_2XF32I(x0, BBE_SHFLI_SWAP_4);//r0,r3,s2,t1
        s3 = BBE_SELN_2XF32I(x2, s2, BBE_SELI_EXTRACT_HI_HALVES);//s3,t0,t2,t3
        s2 = BBE_SELN_2XF32I(s3, s1, BBE_SELI_EXTRACT_4_OF_8_OFF_4);//r2,s1,t0,t3
        s1 = BBE_SELN_2XF32I(s3, s1, BBE_SELI_EXTRACT_4_OF_8_OFF_0);//r1,s0,s3,t2
        BBE_SVN_2XF32_X(s0, pD, 4 * 2 * M);
        BBE_SVN_2XF32_X(s1, pD, 4 * 2 * M + 1 * 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X(s2, pD, 4 * 2 * M + 2 * 2 * BBE_SIMD_WIDTH);

        BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(d2, pD, 2 * BBE_SIMD_WIDTH);

        for (m = 0; m < (M / (3 * BBE_SIMD_WIDTH / 4)); m++)
        {
            BBE_LVN_2XF32_IP(s0, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(s1, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_XP(s2, pD, -5 * 2 * BBE_SIMD_WIDTH);


            BBE_LAVN_2XF32_XP(c0, vH, pH_2, BBE_SIMD_WIDTH);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);

            x0 = BBE_SELN_2XF32I(s0, d0, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(s0, d0, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(s0, d0, BBE_SELI_ROTATE_RIGHT_12);
            x3 = s0;

            BBE_MULAN_2XF32(y0, x0, h0);
            BBE_MULAN_2XF32(y0, x1, h1);
            BBE_MULAN_2XF32(y0, x2, h2);
            BBE_MULAN_2XF32(y0, x3, h3);


            BBE_LVN_2XF32_XP(c0, pH_01, 2 * BBE_SIMD_WIDTH);
            h0 = BBE_REPN_2XF32(c0, 0);
            h1 = BBE_REPN_2XF32(c0, 1);
            h2 = BBE_REPN_2XF32(c0, 2);
            h3 = BBE_REPN_2XF32(c0, 3);
            h4 = BBE_REPN_2XF32(c0, 4);
            h5 = BBE_REPN_2XF32(c0, 5);
            h6 = BBE_REPN_2XF32(c0, 6);
            h7 = BBE_REPN_2XF32(c0, 7);

            x0 = d1;
            x1 = BBE_SELN_2XF32I(s1, d1, BBE_SELI_ROTATE_RIGHT_4);
            x2 = BBE_SELN_2XF32I(s1, d1, BBE_SELI_ROTATE_RIGHT_8);
            x3 = BBE_SELN_2XF32I(s1, d1, BBE_SELI_ROTATE_RIGHT_12);
            x4 = d2;
            x5 = BBE_SELN_2XF32I(s2, d2, BBE_SELI_ROTATE_RIGHT_4);
            x6 = BBE_SELN_2XF32I(s2, d2, BBE_SELI_ROTATE_RIGHT_8);
            x7 = BBE_SELN_2XF32I(s2, d2, BBE_SELI_ROTATE_RIGHT_12);

            BBE_MULAN_2XF32(y0, x0, h0);
            BBE_MULAN_2XF32(y0, x1, h1);
            BBE_MULAN_2XF32(y0, x2, h2);
            BBE_MULAN_2XF32(y0, x3, h3);
            BBE_MULAN_2XF32(y0, x4, h4);
            BBE_MULAN_2XF32(y0, x5, h5);
            BBE_MULAN_2XF32(y0, x6, h6);
            BBE_MULAN_2XF32(y0, x7, h7);


            BBE_SVN_2XF32_IP(s0, pD, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(s1, pD, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(s2, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(d2, pD, 2 * BBE_SIMD_WIDTH);
        }
        BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    }
}// firdecf_proc_D3_MX()
#endif
const tFirFxdxns firdecf_3d_xm = { &firdecf_alloc_d3_mx, firdecf_proc_D3_MX };

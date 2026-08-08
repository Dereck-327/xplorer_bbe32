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
DISCARD_FUN (void,firdecf_proc_D3_M8,( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D ))
#else
void firdecf_proc_D3_M8( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D )
{
    int n;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD = (      xb_vecN_2xf32 *)delayLine;
    const xb_vecN_2xf32 * restrict pH = (const xb_vecN_2xf32 *)coef;

    xb_vecN_2xf32 c0, h0, h1, h2, h3, h4, h5, h6, h7;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5, x6, x7, xx;
    xb_vecN_2xf32 d0, d1, d2;
    xb_vecN_2xf32 s0, s1, s2, s3;
    xb_vecN_2xf32 y0, y1, y2, y3;
    valign vx;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 3 && M == 8);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    vx = BBE_LAN_2XF32_PP(pX);
    BBE_LAVN_2XF32_XP(x0, vx, pX, 8 * 1);
    BBE_LVN_2XF32_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(d0, pD, -2 * BBE_SIMD_WIDTH);
    x0 = BBE_SELN_2XF32I(x0, d0, BBE_SELI_ROTATE_RIGHT_4);//r0,r1,r2,r3
    d2 = BBE_SELN_2XF32I(d1, d1, BBE_SELI_ROTATE_RIGHT_8);
    d0 = BBE_SELN_2XF32I(d1, d1, BBE_SELI_ROTATE_RIGHT_12);
    
    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        BBE_LAN_2XF32_IP(x1, vx, pX);//s0,s1,s2,s3
        BBE_LAN_2XF32_IP(x2, vx, pX);//t0,t1,t2,t3
        BBE_LAN_2XF32_IP(xx, vx, pX);//u0,u1,u2,u3

        s0 = BBE_SELN_2XF32I(x0, x0, BBE_SELI_ROTATE_RIGHT_4);//r1,r2,r3,r0
        s1 = BBE_SELN_2XF32I(x1, s0, BBE_SELI_EXTRACT_LO_HALVES);//r1,r2,s0,s1
        s2 = BBE_SELN_2XF32I(x2, x1, BBE_SELI_ROTATE_RIGHT_4);//s1,s2,s3,t0
        s3 = BBE_SELN_2XF32I(s2, x2, BBE_SELI_INTERLEAVE_4_LO);//t0,s1,t1,s2
        x0 = BBE_SELN_2XF32I(s3, s0, BBE_SELI_EXTRACT_HI_HALVES);//r3,r0,t1,s2
        s0 = BBE_SHFLN_2XF32I(x0, BBE_SHFLI_SWAP_4);//r0,r3,s2,t1
        s3 = BBE_SELN_2XF32I(x2, s2, BBE_SELI_EXTRACT_HI_HALVES);//s3,t0,t2,t3
        s2 = BBE_SELN_2XF32I(s3, s1, BBE_SELI_EXTRACT_4_OF_8_OFF_4);//r2,s1,t0,t3
        s1 = BBE_SELN_2XF32I(s3, s1, BBE_SELI_EXTRACT_4_OF_8_OFF_0);//r1,s0,s3,t2


        BBE_LVN_2XF32_XP(c0, pH, 0);
        h0 = BBE_REPN_2XF32(c0, 0);
        h1 = BBE_REPN_2XF32(c0, 1);
        h2 = BBE_REPN_2XF32(c0, 2);
        h3 = BBE_REPN_2XF32(c0, 3);
        h4 = BBE_REPN_2XF32(c0, 4);
        h5 = BBE_REPN_2XF32(c0, 5);
        h6 = BBE_REPN_2XF32(c0, 6);
        h7 = BBE_REPN_2XF32(c0, 7);

        x0 = BBE_SELN_2XF32I(s2, d2, BBE_SELI_ROTATE_RIGHT_8);
        x1 = BBE_SELN_2XF32I(s0, d0, BBE_SELI_ROTATE_RIGHT_12);
        x2 = BBE_SELN_2XF32I(s1, d1, BBE_SELI_ROTATE_RIGHT_12);
        x3 = BBE_SELN_2XF32I(s2, d2, BBE_SELI_ROTATE_RIGHT_12);
        x4 = s0;
        x5 = s1;
        x6 = s2;
        x7 = BBE_SELN_2XF32I(xx, s0, BBE_SELI_ROTATE_RIGHT_4);//s1,t0,t3,u0

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

        BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
        x0 = xx;
        d0 = s0;
        d1 = s1;
        d2 = s2;
    }
    d1 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_INTERLEAVE_4_HI);
    d0 = BBE_SELN_2XF32I(d1, d2, BBE_SELI_EXTRACT_HI_HALVES);
    d1 = BBE_SELN_2XF32I(x0, x0, BBE_SELI_ROTATE_RIGHT_12);
    BBE_SVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(d1, pD, -2 * BBE_SIMD_WIDTH);
}// firdecf_proc_D3_M8()
#endif
const tFirFxdxns firdecf_3d_8m = { &firdecf_alloc_dx, firdecf_proc_D3_M8 };

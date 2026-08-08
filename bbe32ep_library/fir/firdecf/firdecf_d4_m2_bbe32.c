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
DISCARD_FUN (void,firdecf_proc_D4_M2,( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D ))
#else
void firdecf_proc_D4_M2( complex_float * restrict y,
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

    xb_vecN_2xf32 c0, h0, h1;
    xb_vecN_2xf32 x0, x1, x2, x3;
    xb_vecN_2xf32 d0, d1;
    xb_vecN_2xf32 y0;
    valign vX;

    NASSERT(N > 0 && !(N & 3));
    NASSERT(D == 4 && M == 2);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    vX = BBE_LAN_2XF32_PP(pX);
    BBE_LAVN_2XF32_XP(x0, vX, pX, 8 * 3);
    d0 = BBE_LVN_2XF32_I(pD, 0);
    x0 = BBE_SELN_2XF32I(x0, d0, BBE_SELI_PACK_4);
    BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
    h0 = BBE_REPN_2XF32(c0, 0);
    h1 = BBE_REPN_2XF32(c0, 1);

    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        BBE_LAN_2XF32_IP(x1, vX, pX);
        BBE_LAN_2XF32_IP(x2, vX, pX);
        BBE_LAN_2XF32_IP(x3, vX, pX);

        d0 = BBE_SELN_2XF32I(x2, x0, BBE_SELI_INTERLEAVE_4_LO);
        d1 = BBE_SELN_2XF32I(x3, x1, BBE_SELI_INTERLEAVE_4_LO);
        BBE_DSELN_2XF32I(x1, x0, d1, d0, BBE_DSELI_INTERLEAVE_4);

        y0 = BBE_MULN_2XF32(x0, h0);
        BBE_MULAN_2XF32(y0, x1, h1);

        BBE_SVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_LAN_2XF32_IP(x0, vX, pX);
    }
    BBE_SVN_2XF32_I(x0, pD, 0);
}// firdecf_proc_D4_M2()
#endif
const tFirFxdxns firdecf_4d_2m = { &firdecf_alloc_dx, firdecf_proc_D4_M2 };

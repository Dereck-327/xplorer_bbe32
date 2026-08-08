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
#include "firdec_common.h"

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

/* processing function D==2 && !(M&15) && (M>16) */
static void firdec_proc_D2_MX(  int16_t * restrict y,
                         const int16_t * restrict x,
                         const int16_t * restrict coef,
                               int16_t * restrict delayLine,
                         int M, int N, int D )
{
    xb_vecNx16 * restrict Y;
    xb_vecNx16 * restrict pD;
    xb_vecNx16 * restrict pH;
    const xb_vecNx16 *          X;

    xb_vecNx40 w0, w1;

    xb_vecNx16 cf;
    xb_vecNx16 x0, x1, x2, x3;
    xb_vecNx16 d0, d1, d2, d3;
    xb_vecNx16 y0;

    xb_vecNx16 p00, p01, p02, p03, p04, p05, p06, p07;
    xb_vecNx16 p10, p11, p12, p13, p14, p15, p16, p17;
    uint32_t   q00, q01, q02, q03;
    uint32_t   q10, q11, q12, q13;

    int n, m;

    NASSERT(N>0 && !(N & 7));

    NASSERT(D == 2 && !(M & 15));

    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    Y = (xb_vecNx16*)y;
    X = (const xb_vecNx16*)x;

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        w0 = 0;
        w1 = 0;
        pH = (xb_vecNx16*)coef;
        cf = BBE_LVNX16_I(pH, 0);
        pD = (xb_vecNx16 *)(delayLine);

        // Load 8x2 input samples, CQ15
        BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, X, 2 * BBE_SIMD_WIDTH);
        BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_SVNX16_X(x0, pD, 2 * 2 * M);
        BBE_SVNX16_X(x1, pD, 2 * 2 * M + 2 * BBE_SIMD_WIDTH);

        x0 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 0);
        x1 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 2 * BBE_SIMD_WIDTH);

        d0 = x0;
        d1 = x1;

        __Pragma("ymemory(pD)");
        __Pragma("ymemory(pH)");
        for (m = 0; m<(M / (BBE_SIMD_WIDTH)); m++)
        {
            // Load 8x2 input samples, CQ15
            d2 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH);
            d3 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH + 2 * BBE_SIMD_WIDTH);

            x0 = d0;
            x1 = d2;
            x2 = d1;
            x3 = d3;

            //
            // Coefficients bank 1
            // 

            BBE_SELPCNX16I(p01, p00, x1, x0, 1);
            BBE_SELPCNX16I(p03, p02, x1, x0, 3);
            BBE_SELPCNX16I(p05, p04, x1, x0, 5);
            BBE_SELPCNX16I(p07, p06, x1, x0, 7);

            cf = BBE_LVNX16_I(pH, 0);

            q00 = BBE_EXTRNX16C(cf, 4);
            q01 = BBE_EXTRNX16C(cf, 5);
            q02 = BBE_EXTRNX16C(cf, 6);
            q03 = BBE_EXTRNX16C(cf, 7);


            BBE_MULANX16PR(w0, p01, p00, q00);
            BBE_MULANX16PR(w1, p03, p02, q01);
            BBE_MULANX16PR(w0, p05, p04, q02);
            BBE_MULANX16PR(w1, p07, p06, q03);

            BBE_SELPCNX16I(p11, p10, x3, x2, 0);
            BBE_SELPCNX16I(p13, p12, x3, x2, 2);
            BBE_SELPCNX16I(p15, p14, x3, x2, 4);
            BBE_SELPCNX16I(p17, p16, x3, x2, 6);

            BBE_LVNX16_IP(cf, pH, 2 * BBE_SIMD_WIDTH);
            // Load 8 real filter coefficients, Q15.
            q10 = BBE_EXTRNX16C(cf, 0);
            q11 = BBE_EXTRNX16C(cf, 1);
            q12 = BBE_EXTRNX16C(cf, 2);
            q13 = BBE_EXTRNX16C(cf, 3);

            BBE_MULANX16PR(w1, p17, p16, q13);
            BBE_MULANX16PR(w0, p15, p14, q12);
            BBE_MULANX16PR(w1, p13, p12, q11);
            BBE_MULANX16PR(w0, p11, p10, q10);

            d0 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH);
            d1 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH + 2 * BBE_SIMD_WIDTH);

            BBE_SVNX16_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(d1, pD, 2 * BBE_SIMD_WIDTH);
        }
        //
        // Save 8 output samples.
        //

        w0 = BBE_ADDNX40(w0, w1);
        // CQ15 <- CQ30 - 15 w/ rounding and saturation.
        y0 = BBE_PACKQNX40(w0);

        BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    }
}
const tFirFxdxns firdec_2d_x_8n ={&firdec_alloc_d2_mx,firdec_proc_D2_MX};

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

/* processing function D==3, M==16 */
static void firdec_proc_D3_M16 ( int16_t * restrict y,
                          const int16_t * restrict x,
                          const int16_t * restrict coef,
                                int16_t * restrict delayLine,
                            int M, int N, int D )
{
    xb_vecNx16 * restrict Y;
    const xb_vecNx16 *          X;
    xb_vecNx16 * restrict S;

    xb_vecNx40 w0;

    xb_vecNx16 cf, cf1;
    xb_vecNx16 x0, x1, x2;
    xb_vecNx16 s0, s1, s2;
    xb_vecNx16 d0, d1, d2;
    xb_vecNx16 y0;

    xb_vecNx16 p00, p01, p02, p03, p04, p05;
    xb_vecNx16 p10, p11, p12, p13, p14, p15;
    xb_vecNx16 p20, p21, p22, p23, p24, p25;

    uint32_t   q00, q01, q02;
    uint32_t   q10, q11, q12;
    uint32_t   q20, q21, q22;

    int n;

    NASSERT(N>0 && N % 8 == 0);

    NASSERT(M == 18 && D == 3);

    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    Y = (xb_vecNx16*)y;
    X = (const xb_vecNx16*)x;
    S = (xb_vecNx16*)delayLine;

    //
    // Load the delay line state.
    //

    d0 = BBE_LVNX16_I(S, 0 * 2 * BBE_SIMD_WIDTH);
    d1 = BBE_LVNX16_I(S, 1 * 2 * BBE_SIMD_WIDTH);
    d2 = BBE_LVNX16_I(S, 2 * 2 * BBE_SIMD_WIDTH);

    //
    // Process data.
    //

    cf = BBE_LVNX16_I((const xb_vecNx16*)coef, 0);
    cf1 = BBE_LPNX16_I(coef, 16 * 2);

        __Pragma("ymemory(S)");
        __Pragma("ymemory(X)");
    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        // Load 8x3 input samples, CQ15
        BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, X, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, X, 2 * BBE_SIMD_WIDTH);

        // Transposition 8x3 -> 3x8 (3 banks, each of 8 samples).
        BBE_DSELNX16I(s2, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_C3_STEP_0);
        BBE_DSELNX16I_H(s2, x1, x1, x2, BBE_DSELI_DEINTERLEAVE_C3_STEP_1);
        BBE_DSELNX16I(s1, s0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);

        d2 = BBE_LVNX16_I(S, 2 * 2 * BBE_SIMD_WIDTH);

        BBE_SVNX16_I(s2, S, 2 * 2 * BBE_SIMD_WIDTH);

        //
        // Coefficients bank 2
        //

        BBE_SELPCNX16I(p01, p00, s0, d0, 3);
        BBE_SELPCNX16I(p03, p02, s0, d0, 5);
        BBE_SELPCNX16I(p05, p04, s0, d0, 7);

        d0 = s0;

        // Load 6 real filter coefficients, Q15.
        q00 = BBE_EXTRNX16C(cf, 6);
        q01 = BBE_EXTRNX16C(cf, 7);
        q02 = BBE_EXTRNX16C(cf1, 0);

        w0 = BBE_MULNX16PR(p01, p00, q00);
        BBE_MULANX16PR(w0, p03, p02, q01);
        BBE_MULANX16PR(w0, p05, p04, q02);

        //                            
        // Coefficients bank 0
        //                            

        BBE_SELPCNX16I(p11, p10, s1, d1, 2);
        BBE_SELPCNX16I(p13, p12, s1, d1, 4);
        BBE_SELPCNX16I(p15, p14, s1, d1, 6);

        d1 = s1;

        q10 = BBE_EXTRNX16C(cf, 0);
        q11 = BBE_EXTRNX16C(cf, 1);
        q12 = BBE_EXTRNX16C(cf, 2);

        BBE_MULANX16PR(w0, p11, p10, q10);
        BBE_MULANX16PR(w0, p13, p12, q11);
        BBE_MULANX16PR(w0, p15, p14, q12);

        //
        // Coefficients bank 1
        //

        BBE_SELPCNX16I(p21, p20, s2, d2, 2);
        BBE_SELPCNX16I(p23, p22, s2, d2, 4);
        BBE_SELPCNX16I(p25, p24, s2, d2, 6);

        q20 = BBE_EXTRNX16C(cf, 3);
        q21 = BBE_EXTRNX16C(cf, 4);
        q22 = BBE_EXTRNX16C(cf, 5);

        BBE_MULANX16PR(w0, p21, p20, q20);
        BBE_MULANX16PR(w0, p23, p22, q21);
        BBE_MULANX16PR(w0, p25, p24, q22);

        //
        // Save 8 output samples.
        //
        // CQ15 <- CQ30 - 15 w/ rounding and saturation.
        y0 = BBE_PACKQNX40(w0);

        BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    }

    //
    // Save the delay line state.
    //

    BBE_SVNX16_I(d0, (xb_vecNx16*)delayLine, 0 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(d1, (xb_vecNx16*)delayLine, 1 * 2 * BBE_SIMD_WIDTH);
} // firdec_proc_D3_M16()
const tFirFxdxns firdec_3d_16_8n = { &firdec_alloc_d3_m16, firdec_proc_D3_M16 };

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

/* processing function D==3, M==32 */
static void firdec_proc_D3_M32 ( int16_t * restrict y,
                          const int16_t * restrict x,
                          const int16_t * restrict coef,
                                int16_t * restrict delayLine,
                            int M, int N, int D )
{
          xb_vecNx16 * restrict Y;
    const xb_vecNx16 *          X;
          xb_vecNx16 * restrict S;

    xb_vecNx40 w0;

    xb_vecNx16 cf0, cf1, cf2;
    xb_vecNx16 x0, x1, x2;
    xb_vecNx16 s0, s1, s2;
    xb_vecNx16 d00, d10, d01, d11, d02, d12;
    xb_vecNx16 y0;

    xb_vecNx16 p00, p01, p02, p03, p04, p05, p06, p07, p08, p09, p0a, p0b;
    xb_vecNx16 /*p10, p11,*/ p12, p13, p14, p15, p16, p17, p18, p19, p1a, p1b;
    xb_vecNx16 p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p2a, p2b;

    uint32_t   q00, q01, q02, q03, q04, q05;
    uint32_t   /*q10,*/ q11, q12, q13, q14, q15;
    uint32_t   q20, q21, q22, q23, q24, q25;

    int n;

    NASSERT(N>0 && N % 8 == 0);

    NASSERT(M == 36 && D == 3);

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

    d00 = BBE_LVNX16_I(S, 0 * 2 * BBE_SIMD_WIDTH);
    d01 = BBE_LVNX16_I(S, 1 * 2 * BBE_SIMD_WIDTH);
    d02 = BBE_LVNX16_I(S, 2 * 2 * BBE_SIMD_WIDTH);

    d10 = BBE_LVNX16_I(S, 3 * 2 * BBE_SIMD_WIDTH);
    d11 = BBE_LVNX16_I(S, 4 * 2 * BBE_SIMD_WIDTH);
    d12 = BBE_LVNX16_I(S, 5 * 2 * BBE_SIMD_WIDTH);

    //
    // Process data.
    //

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

        cf0 = BBE_LVNX16_I((const xb_vecNx16*)coef, 0 * 2 * BBE_SIMD_WIDTH);
        cf1 = BBE_LVNX16_I((const xb_vecNx16*)coef, 1 * 2 * BBE_SIMD_WIDTH);
        cf2 = BBE_LVNX16_I((const xb_vecNx16*)coef, 32 * 2);//cf2 = BBE_LV4X16_I(coef, 32 * 2);

        d00 = BBE_LVNX16_I(S, 0 * 2 * BBE_SIMD_WIDTH);
        d01 = BBE_LVNX16_I(S, 1 * 2 * BBE_SIMD_WIDTH);
        d02 = BBE_LVNX16_I(S, 2 * 2 * BBE_SIMD_WIDTH);

        d10 = BBE_LVNX16_I(S, 3 * 2 * BBE_SIMD_WIDTH);
        d11 = BBE_LVNX16_I(S, 4 * 2 * BBE_SIMD_WIDTH);
        d12 = BBE_LVNX16_I(S, 5 * 2 * BBE_SIMD_WIDTH);

        BBE_SVNX16_I(d10, S, 0 * 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_I(d11, S, 1 * 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_I(d12, S, 2 * 2 * BBE_SIMD_WIDTH);

        BBE_SVNX16_I(s0, S, 3 * 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_I(s1, S, 4 * 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_I(s2, S, 5 * 2 * BBE_SIMD_WIDTH);

        //
        // Coefficients bank 2
        //
#if 0
        BBE_SELPCNX16I(p01, p00, d10, d00, 5);
        BBE_SELPCNX16I(p03, p02, d10, d00, 7);
        BBE_SELPCNX16I(p05, p04, s0, d10, 1);
        BBE_SELPCNX16I(p07, p06, s0, d10, 3);
        BBE_SELPCNX16I(p09, p08, s0, d10, 5);
        BBE_SELPCNX16I(p0b, p0a, s0, d10, 7);
#else
        BBE_SELPCNX16I(p01, p00, d10, d00, 5);
        p02 = BBE_SELNX16I(d10, d00, BBE_SELI_ROTATE_RIGHT_14);
        p03 = d10;
        //BBE_SELPCNX16I( p03, p02, d10, d00, 7 );
        BBE_SELPCNX16I(p05, p04, s0, d10, 1);
        BBE_SELPCNX16I(p07, p06, s0, d10, 3);
        BBE_SELPCNX16I(p09, p08, s0, d10, 5);
        BBE_SELPCNX16I(p0b, p0a, s0, d10, 7);
#endif

        // Load 12 real filter coefficients, Q15.
        q00 = BBE_EXTRNX16C(cf1, 4);
        q01 = BBE_EXTRNX16C(cf1, 5);
        q02 = BBE_EXTRNX16C(cf1, 6);
        q03 = BBE_EXTRNX16C(cf1, 7);
        q04 = BBE_EXTRNX16C(cf2, 0);
        q05 = BBE_EXTRNX16C(cf2, 1);

        w0 = BBE_MULNX16PR(p01, p00, q00);
        BBE_MULANX16PR(w0, p03, p02, q01);
        BBE_MULANX16PR(w0, p05, p04, q02);
        BBE_MULANX16PR(w0, p07, p06, q03);
        BBE_MULANX16PR(w0, p09, p08, q04);
        BBE_MULANX16PR(w0, p0b, p0a, q05);

        //                            
        // Coefficients bank 0
        //                            

#if 1
        //BBE_SELPCNX16I( p11, p10, d11, d01, 4 ); 
        BBE_SELPCNX16I(p13, p12, d11, d01, 6);
        BBE_SELPCNX16I(p15, p14, s1, d11, 0);
        BBE_SELPCNX16I(p17, p16, s1, d11, 2);
        BBE_SELPCNX16I(p19, p18, s1, d11, 4);
        BBE_SELPCNX16I(p1b, p1a, s1, d11, 6);
#else
        //BBE_SELPCNX16I( p11, p10, d11, d01, 4 ); 
        BBE_SELPCNX16I(p13, p12, d11, d01, 6);
        //BBE_SELPCNX16I( p15, p14,  s1, d11, 0 );
        p14 = d11;
        p15 = BBE_SELNX16I(s1, d11, BBE_SELI_ROTATE_RIGHT_2);
        BBE_SELPCNX16I(p17, p16, s1, d11, 2);
        BBE_SELPCNX16I(p19, p18, s1, d11, 4);
        BBE_SELPCNX16I(p1b, p1a, s1, d11, 6);
#endif

        //q10 = BBE_EXTRNX16C( cf0, 0 );
        q11 = BBE_EXTRNX16C(cf0, 1);
        q12 = BBE_EXTRNX16C(cf0, 2);
        q13 = BBE_EXTRNX16C(cf0, 3);
        q14 = BBE_EXTRNX16C(cf0, 4);
        q15 = BBE_EXTRNX16C(cf0, 5);

        //BBE_MULANX16PR( w0, p11, p10, q10 );
        BBE_MULANX16PR(w0, p13, p12, q11);
        BBE_MULANX16PR(w0, p15, p14, q12);
        BBE_MULANX16PR(w0, p17, p16, q13);
        BBE_MULANX16PR(w0, p19, p18, q14);
        BBE_MULANX16PR(w0, p1b, p1a, q15);

        //
        // Coefficients bank 1
        //
#if 0
        BBE_SELPCNX16I(p21, p20, d12, d02, 4);
        BBE_SELPCNX16I(p23, p22, d12, d02, 6);
        BBE_SELPCNX16I(p25, p24, s2, d12, 0);
        BBE_SELPCNX16I(p27, p26, s2, d12, 2);
        BBE_SELPCNX16I(p29, p28, s2, d12, 4);
        BBE_SELPCNX16I(p2b, p2a, s2, d12, 6);
#else
        BBE_SELPCNX16I(p21, p20, d12, d02, 4);
        BBE_SELPCNX16I(p23, p22, d12, d02, 6);
        //BBE_SELPCNX16I( p25, p24,  s2, d12, 0 );
        p24 = d12;
        p25 = BBE_SELNX16I(s2, d12, BBE_SELI_ROTATE_RIGHT_2);
        BBE_SELPCNX16I(p27, p26, s2, d12, 2);
        BBE_SELPCNX16I(p29, p28, s2, d12, 4);
        BBE_SELPCNX16I(p2b, p2a, s2, d12, 6);
#endif

        q20 = BBE_EXTRNX16C(cf0, 6);
        q21 = BBE_EXTRNX16C(cf0, 7);
        q22 = BBE_EXTRNX16C(cf1, 0);
        q23 = BBE_EXTRNX16C(cf1, 1);
        q24 = BBE_EXTRNX16C(cf1, 2);
        q25 = BBE_EXTRNX16C(cf1, 3);

        BBE_MULANX16PR(w0, p21, p20, q20);
        BBE_MULANX16PR(w0, p23, p22, q21);
        BBE_MULANX16PR(w0, p25, p24, q22);
        BBE_MULANX16PR(w0, p27, p26, q23);
        BBE_MULANX16PR(w0, p29, p28, q24);
        BBE_MULANX16PR(w0, p2b, p2a, q25);

        //
        // Save 8 output samples.
        //
        // CQ15 <- CQ30 - 15 w/ rounding and saturation.
        y0 = BBE_PACKQNX40(w0);

        BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
    }
} // firdec_proc_D3_M32()
const tFirFxdxns firdec_3d_32_8n = {&firdec_alloc_d3_m32, firdec_proc_D3_M32};

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
    Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "cxfirf_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR Filter

Computes a complex FIR filter (direct-form) using complex IR stored in 
vector h. The complex data input is stored in vector x. The filter output
result is stored in vector y. The filter calculates N output samples using
M coefficients and requires last M+N-1 samples in the delay line. 

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by cxfir[f]_algDelay(M).

Representation:
cxfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
cxfirf  IEEE-754 Std. single precision floating-point format for filter
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[2*M]  Filter coefficients; h[0]+j*h[1] is to be multiplied by the newest
        sample,
N       Length of sample block
M       Length of filter
x[N]    Complex input samples
Output:
y[N]    Complex output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cxfir) or 4 (cxfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, cxfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,cxfirf_process_4m,(complex_float * restrict y, const complex_float * restrict x, const float32_t * restrict coef, float32_t * restrict delayLine, int M, int N))
#else
void cxfirf_process_4m(complex_float * restrict y, const complex_float * restrict x, const float32_t * restrict coef, float32_t * restrict delayLine, int M, int N)
{
    int n, m;
    const xb_vecN_2xf32 * restrict pX = (const xb_vecN_2xf32 *)x;
          xb_vecN_2xf32 * restrict pY = (      xb_vecN_2xf32 *)y;
          xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pH;

    xb_vecN_2xf32 d0, d1, d2, c0;
    xb_vecN_2xf32 x0, x1, x2, x3;
    xb_vecN_2xf32 h0, h1, h2, h3;
    xb_vecN_2xf32 A0, A1, A2, A3;

    NASSERT(N > 0 && N % 4 == 0);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delayLine);
    NASSERT_ALIGN32(coef);

    //__Pragma("loop_count min=1");
    for (n = 0; n < N / (BBE_SIMD_WIDTH / 4); n++)
    {
        pD = (xb_vecN_2xf32 *)(delayLine);
        pH = (const xb_vecN_2xf32 *)coef;

        BBE_LVN_2XF32_IP(d2, pX, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X(d2, pD, 2 * 4 * M);

        d0 = BBE_LVN_2XF32_I(pD, 0 * 2 * BBE_SIMD_WIDTH);

        A0 = BBE_ZERON_2XF32();
        A1 = BBE_ZERON_2XF32();
        A2 = BBE_ZERON_2XF32();
        A3 = BBE_ZERON_2XF32();

        //__Pragma("loop_count min=1");
        for (m = 0; m < M / (BBE_SIMD_WIDTH / 4); m++)
        {
            d1 = BBE_LVN_2XF32_I(pD, 1 * 2 * BBE_SIMD_WIDTH);

            BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
            h0 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_REP_0X4);
            h1 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_REP_1X4);
            h2 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_REP_2X4);
            h3 = BBE_SHFLN_2XF32I(c0, BBE_SHFLI_REP_3X4);

            x0 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(d1, d0, BBE_SELI_ROTATE_RIGHT_12);
            x3 = d1;

            BBE_MULMASN_2XF32(A0, x0, h0, 0, 4);
            BBE_MULMASN_2XF32(A1, x0, h0, 1, 11);
            BBE_MULMASN_2XF32(A2, x1, h1, 0, 4);
            BBE_MULMASN_2XF32(A3, x1, h1, 1, 11);
            BBE_MULMASN_2XF32(A0, x2, h2, 0, 4);
            BBE_MULMASN_2XF32(A1, x2, h2, 1, 11);
            BBE_MULMASN_2XF32(A2, x3, h3, 0, 4);
            BBE_MULMASN_2XF32(A3, x3, h3, 1, 11);

            BBE_SVN_2XF32_I(d1, pD, 0 * 2 * BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_IP(d0, pD, 2 * BBE_SIMD_WIDTH);
        }
        A0 = BBE_ADDN_2XF32(A0, A1);
        A2 = BBE_ADDN_2XF32(A2, A3);
        A0 = BBE_ADDN_2XF32(A0, A2);
        BBE_SVN_2XF32_IP(A0, pY, 2 * BBE_SIMD_WIDTH);
    }
} //cxfir_process_8m()
#endif

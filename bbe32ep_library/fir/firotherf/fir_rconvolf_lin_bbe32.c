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
    Circular/Linear Convolution for Real Data
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fir.h"

/*-------------------------------------------------------------------------
Circular/Linear Convolution for Real Data

Compute circular or linear convolution between real-valued vectors x (of
length N) and y (of length M) resulting in vector r of N (circular) or
N+M-1 (linear) real values.

MATLAB code for circular convolution:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-(1:M),N)).*y(1:M));
  end

MATLAB code for linear convolution:
  for n=1:N+M-1
    ix = max(n-N+1,1):min(n,M);
    r(n) = 2*sum(x(n+1-ix).*y(ix));
  end

Representation: IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]    Left-hand data sequence
y[M]    Right-hand data sequence
Output:
r[]     Output data. Size is N for circular convolution, or N+M-1 
        for linear convolution
Restrictions:
x,y,r   Must not overlap
x,y,r   Must be aligned on 32-byte boundary
N,M     Multiples of 8
N>=M    N must be greater than or equal to M
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,fir_rconvolf_lin,(float32_t * restrict r,
                const float32_t * restrict x,
                const float32_t * restrict y,
                int N, int M))
#else
void fir_rconvolf_lin(float32_t * restrict r,
                const float32_t * restrict x,
                const float32_t * restrict y,
                int N, int M)
{
    int n, m;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pY;
          xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pS;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5, x6, x7, x_pre, x_cur;
    xb_vecN_2xf32 y0, y1, y2, y3, y4, y5, y6, y7, y_cur;
    // Masks for predicated vector stores.
    static const uint8_t cst[] = { 0x7f, 0xff };
    vboolN_2 vb0;
    uint32_t ix;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(r, (2 * BBE_SIMD_WIDTH));
    NASSERT(N > 0 && N % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(M > 0 && M % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(N >= M);

    pR = (xb_vecN_2xf32 *)r;
    pX = (const xb_vecN_2xf32 *)x;

    for (n = 0; n < M / (BBE_SIMD_WIDTH / 2); n++)
    {
        xb_vecN_2xf32 A = 0.f;
        x_pre = 0.f;
        pS = pX;
        pY = (const xb_vecN_2xf32 *)((uintptr_t)y + n * 2 * BBE_SIMD_WIDTH);

        for (m = 0; m < n + 1; m++)
        {
            BBE_LVN_2XF32_IP(y_cur, pY, -2 * BBE_SIMD_WIDTH);
            y_cur = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REVERSE_2);
            y0 = BBE_REPN_2XF32(y_cur, 0);
            y1 = BBE_REPN_2XF32(y_cur, 1);
            y2 = BBE_REPN_2XF32(y_cur, 2);
            y3 = BBE_REPN_2XF32(y_cur, 3);
            y4 = BBE_REPN_2XF32(y_cur, 4);
            y5 = BBE_REPN_2XF32(y_cur, 5);
            y6 = BBE_REPN_2XF32(y_cur, 6);
            y7 = BBE_REPN_2XF32(y_cur, 7);

            BBE_LVN_2XF32_IP(x_cur, pS, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_2);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_6);
            x3 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x4 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_10);
            x5 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x6 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_14);
            x7 = x_cur;
            x_pre = x_cur;

            BBE_MULAN_2XF32(A, y0, x0);
            BBE_MULAN_2XF32(A, y1, x1);
            BBE_MULAN_2XF32(A, y2, x2);
            BBE_MULAN_2XF32(A, y3, x3);
            BBE_MULAN_2XF32(A, y4, x4);
            BBE_MULAN_2XF32(A, y5, x5);
            BBE_MULAN_2XF32(A, y6, x6);
            BBE_MULAN_2XF32(A, y7, x7);
        }
        BBE_SVN_2XF32_IP(A, pR, 2 * BBE_SIMD_WIDTH);
    }

    pY = (const xb_vecN_2xf32 *)((uintptr_t)y + sizeof(float32_t) * M - 2 * BBE_SIMD_WIDTH);

    for (n = 0; n < (N - M) / (BBE_SIMD_WIDTH / 2); n++)
    {
        xb_vecN_2xf32 A = 0.f;
        BBE_LVN_2XF32_IP(x_pre, pX, 2 * BBE_SIMD_WIDTH);
        pS = pX;

        for (m = 0; m < M / (BBE_SIMD_WIDTH / 2); m++)
        {
            BBE_LVN_2XF32_IP(y_cur, pY, -2 * BBE_SIMD_WIDTH);
            y_cur = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REVERSE_2);
            y0 = BBE_REPN_2XF32(y_cur, 0);
            y1 = BBE_REPN_2XF32(y_cur, 1);
            y2 = BBE_REPN_2XF32(y_cur, 2);
            y3 = BBE_REPN_2XF32(y_cur, 3);
            y4 = BBE_REPN_2XF32(y_cur, 4);
            y5 = BBE_REPN_2XF32(y_cur, 5);
            y6 = BBE_REPN_2XF32(y_cur, 6);
            y7 = BBE_REPN_2XF32(y_cur, 7);

            BBE_LVN_2XF32_IP(x_cur, pS, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_2);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_6);
            x3 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x4 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_10);
            x5 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x6 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_14);
            x7 = x_cur;
            x_pre = x_cur;

            BBE_MULAN_2XF32(A, y0, x0);
            BBE_MULAN_2XF32(A, y1, x1);
            BBE_MULAN_2XF32(A, y2, x2);
            BBE_MULAN_2XF32(A, y3, x3);
            BBE_MULAN_2XF32(A, y4, x4);
            BBE_MULAN_2XF32(A, y5, x5);
            BBE_MULAN_2XF32(A, y6, x6);
            BBE_MULAN_2XF32(A, y7, x7);
        }
        pY = (const xb_vecN_2xf32*)((uintptr_t)pY + sizeof(float32_t) * M);
        //pY = (const xb_vecN_2xf32*)XT_ADDX4(M, (uintptr_t)pY);

        BBE_SVN_2XF32_IP(A, pR, 2 * BBE_SIMD_WIDTH);
    }

    pY = (const xb_vecN_2xf32 *)y;

    for (n = M / (BBE_SIMD_WIDTH / 2) - 1; n >= 0; n--)
    {
        xb_vecN_2xf32 A = 0.f;
        pX = (const xb_vecN_2xf32*)((uintptr_t)x + sizeof(float32_t) * N - 2 * BBE_SIMD_WIDTH);
        x_cur = 0.f;
        BBE_LVN_2XF32_IP(y_cur, pY, 2 * BBE_SIMD_WIDTH);
        pS = pY;

        for (m = 0; m < n + 1; m++)
        {
            y_cur = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REVERSE_2);
            y0 = BBE_REPN_2XF32(y_cur, 0);
            y1 = BBE_REPN_2XF32(y_cur, 1);
            y2 = BBE_REPN_2XF32(y_cur, 2);
            y3 = BBE_REPN_2XF32(y_cur, 3);
            y4 = BBE_REPN_2XF32(y_cur, 4);
            y5 = BBE_REPN_2XF32(y_cur, 5);
            y6 = BBE_REPN_2XF32(y_cur, 6);
            y7 = BBE_REPN_2XF32(y_cur, 7);

            BBE_LVN_2XF32_IP(x_pre, pX, -2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_2);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_6);
            x3 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x4 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_10);
            x5 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x6 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_14);
            x7 = x_cur;
            x_cur = x_pre;

            BBE_MULAN_2XF32(A, y0, x0);
            BBE_MULAN_2XF32(A, y1, x1);
            BBE_MULAN_2XF32(A, y2, x2);
            BBE_MULAN_2XF32(A, y3, x3);
            BBE_MULAN_2XF32(A, y4, x4);
            BBE_MULAN_2XF32(A, y5, x5);
            BBE_MULAN_2XF32(A, y6, x6);
            BBE_MULAN_2XF32(A, y7, x7);

            BBE_LVN_2XF32_IP(y_cur, pS, 2 * BBE_SIMD_WIDTH);
        }
        ix = XT_MIN(1, n);
        vb0 = BBE_LBN_2_I((const vboolN_2*)(cst + ix), 0);
        BBE_SVN_2XF32T_IP(A, pR, 2 * BBE_SIMD_WIDTH, vb0);
    }
} //fir_rconvolf_lin()
#endif

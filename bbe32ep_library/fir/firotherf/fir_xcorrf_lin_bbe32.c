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
    Circular/Linear Cross-Correlation for Complex Data
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
Circular/Linear Cross-Correlation for Complex Data

Estimates the circular/linear cross-correlation between complex-valued
vectors x (of length N) and y (of length M) resulting in vector r. It
is similar to convolution, but y is read in opposite direction.

MATLAB code for circular cross-correlation:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-1+(0:M-1),N)).*conj(y(1:M)));
  end

MATLAB code for linear cross-correlation:
  for n=1:M+N-1
    ix = max(n-M+1,1):min(n,N);
    r(n) = 2*sum(x(ix).*conj(y(M-(n-ix))));
  end;

Representation:
fir_xcorr   Signed fixed-point format
            Input data are 16-bit Q15, output data are 32-bit Q31
fir_xcorrf  IEEE-754 Std. single precision floating-point format for
            input/output data

Parameters:
Input:
x[N]        Left-hand complex data sequence
y[M]        Right-hand complex data sequence
N           Length of x
M           Length of y
Output:
r[]         Complex output data. Size is N for circular cross-correlation,
            or (N+M-1) for linear cross-correlation
Restrictions:
x,y,r       Must not overlap
x,y,r       Must be aligned on 32-byte boundary
N,M         Multiples of 8 (fir_xcorr_circ, fir_xcorr_lin) or 4 
            (fir_xcorrf_circ, fir_xcorrf_lin)
N>=M        N must be greater than or equal to M

NOTES:
1.fir_xcorr[f]_lin returns the same output as MATLAB xcorr but omits first 
N-M zeros
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,fir_xcorrf_lin,(complex_float * restrict r,
              const complex_float * restrict x,
              const complex_float * restrict y,
              int N, int M))
#else
void fir_xcorrf_lin(complex_float * restrict r,
              const complex_float * restrict x,
              const complex_float * restrict y,
              int N, int M)
{  
    int n, m;
    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pY;
    xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pS;
    xb_vecN_2xf32 x0, x1, x2, x3, x_pre, x_cur;
    xb_vecN_2xf32 y0, y1, y2, y3, y_cur;
    // Masks for predicated vector stores.
    static const uint8_t cst[] = { 0x3f, 0xff };
    vboolN_2 vb0;
    uint32_t ix;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(r, (2 * BBE_SIMD_WIDTH));
    NASSERT(N > 0 && N % (BBE_SIMD_WIDTH / 4) == 0);
    NASSERT(M > 0 && M % (BBE_SIMD_WIDTH / 4) == 0);
    NASSERT(N >= M);

    pR = (xb_vecN_2xf32 *)r;
    pX = (const xb_vecN_2xf32 *)x;
    pY = (const xb_vecN_2xf32 *)((uintptr_t)y + 2 * sizeof(float32_t) * M - 2 * BBE_SIMD_WIDTH);

    for (n = 0; n < M / (BBE_SIMD_WIDTH / 4); n++)
    {
        xb_vecN_2xf32 A0 = 0.f, A1 = 0.f, A2 = 0.f, A3 = 0.f;
        x_pre = 0.f;
        pS = pX;

        for (m = 0; m < n + 1; m++)
        {
            BBE_LVN_2XF32_IP(y_cur, pY, 2 * BBE_SIMD_WIDTH);
            y0 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_0X4);
            y1 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_1X4);
            y2 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_2X4);
            y3 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_3X4);

            BBE_LVN_2XF32_IP(x_cur, pS, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x3 = x_cur;
            x_pre = x_cur;

            BBE_MULMASN_2XF32(A0, y0, x0, 0, 4);
            BBE_MULMASN_2XF32(A1, y0, x0, 2, 11);
            BBE_MULMASN_2XF32(A2, y1, x1, 0, 4);
            BBE_MULMASN_2XF32(A3, y1, x1, 2, 11);
            BBE_MULMASN_2XF32(A0, y2, x2, 0, 4);
            BBE_MULMASN_2XF32(A1, y2, x2, 2, 11);
            BBE_MULMASN_2XF32(A2, y3, x3, 0, 4);
            BBE_MULMASN_2XF32(A3, y3, x3, 2, 11);
        }
        pY = (const xb_vecN_2xf32*)((uintptr_t)pY - (n + 2) * 2 * BBE_SIMD_WIDTH);

        A0 = BBE_ADDN_2XF32(A0, A1);
        A2 = BBE_ADDN_2XF32(A2, A3);
        A0 = BBE_ADDN_2XF32(A0, A2);

        BBE_SVN_2XF32_IP(A0, pR, 2 * BBE_SIMD_WIDTH);
    }

    pY = (const xb_vecN_2xf32 *)y;

    for (n = 0; n < (N - M) / (BBE_SIMD_WIDTH / 4); n++)
    {
        xb_vecN_2xf32 A0 = 0.f, A1 = 0.f, A2 = 0.f, A3 = 0.f;
        BBE_LVN_2XF32_IP(x_pre, pX, 2 * BBE_SIMD_WIDTH);
        pS = pX;

        for (m = 0; m < M / (BBE_SIMD_WIDTH / 4); m++)
        {
            BBE_LVN_2XF32_IP(y_cur, pY, 2 * BBE_SIMD_WIDTH);
            y0 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_0X4);
            y1 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_1X4);
            y2 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_2X4);
            y3 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_3X4);

            BBE_LVN_2XF32_IP(x_cur, pS, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x3 = x_cur;
            x_pre = x_cur;

            BBE_MULMASN_2XF32(A0, y0, x0, 0, 4);
            BBE_MULMASN_2XF32(A1, y0, x0, 2, 11);
            BBE_MULMASN_2XF32(A2, y1, x1, 0, 4);
            BBE_MULMASN_2XF32(A3, y1, x1, 2, 11);
            BBE_MULMASN_2XF32(A0, y2, x2, 0, 4);
            BBE_MULMASN_2XF32(A1, y2, x2, 2, 11);
            BBE_MULMASN_2XF32(A2, y3, x3, 0, 4);
            BBE_MULMASN_2XF32(A3, y3, x3, 2, 11);
        }
        pY = (const xb_vecN_2xf32*)((uintptr_t)pY - 2 * sizeof(float32_t) * M);

        A0 = BBE_ADDN_2XF32(A0, A1);
        A2 = BBE_ADDN_2XF32(A2, A3);
        A0 = BBE_ADDN_2XF32(A0, A2);

        BBE_SVN_2XF32_IP(A0, pR, 2 * BBE_SIMD_WIDTH);
    }

    pY = (const xb_vecN_2xf32 *)((uintptr_t)y + 2 * sizeof(float32_t) * M - 2 * BBE_SIMD_WIDTH);

    for (n = M / (BBE_SIMD_WIDTH / 4) - 1; n >= 0; n--)
    {
        xb_vecN_2xf32 A0 = 0.f, A1 = 0.f, A2 = 0.f, A3 = 0.f;
        pX = (const xb_vecN_2xf32*)((uintptr_t)x + 2 * sizeof(float32_t) * N - 2 * BBE_SIMD_WIDTH);
        x_cur = 0.f;
        pS = pY;
        pY = (const xb_vecN_2xf32*)((uintptr_t)pY - 2 * BBE_SIMD_WIDTH);

        for (m = 0; m < n + 1; m++)
        {
            BBE_LVN_2XF32_IP(y_cur, pS, -2 * BBE_SIMD_WIDTH);
            y0 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_0X4);
            y1 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_1X4);
            y2 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_2X4);
            y3 = BBE_SHFLN_2XF32I(y_cur, BBE_SHFLI_REP_3X4);

            BBE_LVN_2XF32_IP(x_pre, pX, -2 * BBE_SIMD_WIDTH);
            x0 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_4);
            x1 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_8);
            x2 = BBE_SELN_2XF32I(x_cur, x_pre, BBE_SELI_ROTATE_RIGHT_12);
            x3 = x_cur;
            x_cur = x_pre;

            BBE_MULMASN_2XF32(A0, y0, x0, 0, 4);
            BBE_MULMASN_2XF32(A1, y0, x0, 2, 11);
            BBE_MULMASN_2XF32(A2, y1, x1, 0, 4);
            BBE_MULMASN_2XF32(A3, y1, x1, 2, 11);
            BBE_MULMASN_2XF32(A0, y2, x2, 0, 4);
            BBE_MULMASN_2XF32(A1, y2, x2, 2, 11);
            BBE_MULMASN_2XF32(A2, y3, x3, 0, 4);
            BBE_MULMASN_2XF32(A3, y3, x3, 2, 11);
        }
        A0 = BBE_ADDN_2XF32(A0, A1);
        A2 = BBE_ADDN_2XF32(A2, A3);
        A0 = BBE_ADDN_2XF32(A0, A2);

        ix = XT_MIN(1, n);
        vb0 = BBE_LBN_2_I((const vboolN_2*)(cst + ix), 0);
        BBE_SVN_2XF32T_IP(A0, pR, 2 * BBE_SIMD_WIDTH, vb0);
    }
} // fir_xcorrf_lin()
#endif

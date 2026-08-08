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
    Circular/Linear Convolution for Complex Data
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
Circular/Linear Convolution for Complex Data

Compute circular or linear convolution between complex-valued vectors x (of
length N) and y (of length M) resulting in vector r of N (circular) or
N+M-1 (linear) complex values.

MATLAB code for circular convolution:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-(1:M),N)).*y(1:M));
  end

MATLAB code for linear convolution:
  for n=1:N+M-1
    ix = max(n-N+1,1):min(n,M);
    r(n) = 2*sum(x(n+1-ix).*y(ix));
  end

Representation:
fir_convol   Signed fixed-point format
             Input data are 16-bit Q15, output data are 32-bit Q31
fir_convolf  IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]         Left-hand complex data sequence
y[M]         Right-hand complex data sequence
Output:
r[]          Complex output data. Size is N for circular convolution, 
             or (N+M-1) for linear convolution
Restrictions:
x,y,r        Must not overlap
x,y,r        Must be aligned on 32-byte boundary
N,M          Multiples of 8 (fir_convol_circ, fir_convol_lin) or 4 
             (fir_convolf_circ, fir_convolf_lin)
N>=M         N must be greater than or equal to M
-------------------------------------------------------------------------*/

#if !HAVE_MULPC
DISCARD_FUN(void, fir_convol_lin, (complex_fract32 * restrict r,
                             const complex_fract16 * restrict x,
                             const complex_fract16 * restrict y,
                             int N, int M))
#else

void fir_convol_lin ( complex_fract32 * restrict r,
                const complex_fract16 * restrict x,
                const complex_fract16 * restrict y,
                int N, int M )
{
    xb_vecNx16 * restrict R;
    const xb_vecNx16 *          X;
    const xb_vecNx16 *          Y;
    const xb_vecNx16 *          S;

    xb_vecNx16 x0, x1, y0, r0, r1;

    xb_vecNx40 w0;

    vboolN vb0;

    uint32_t ix;

    int n, m;

    // Masks for predicated vector stores.
    static const uint16_t cst[] = { 0x0fff, 0xffff };

    NASSERT_ALIGN32(r);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(y);

    NASSERT(N>0 && !(N & 7));
    NASSERT(M>0 && !(M & 7));
    NASSERT(N >= M);

    //----------------------------------------------------------------------------
    // Part I. Compute r(1:M).

    R = (xb_vecNx16*)r;
    X = (const xb_vecNx16*)x;

    //
    // Outer loop is unrolled 8 times, i.e. we compute 8 convolution samples at
    // each iteration.
    //

    for (n = 0; n<M / (BBE_SIMD_WIDTH / 2); n++)
    {
        // Right-hand data pointer looks at y(1+8*n).
        Y = (const xb_vecNx16*)((uintptr_t)y + n * 2 * BBE_SIMD_WIDTH);
        // Left-hand pointer data is reset to x(1) at each iteration.
        S = X;

        x0 = 0;

        w0 = 0;

        //
        // Inner loop is unrolled 8 times, that is, we compute 8 products for
        // each of 16 convolutions.
        //

        __Pragma("ymemory(S)");
        __Pragma("ymemory(Y)");
        for (m = 0; m<n + 1; m++)
        {
            // Load x(8*m+1:8)
            BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);
            // Load y(8*n-8*m+1:8)
            BBE_LVNX16_IP(y0, Y, -2 * BBE_SIMD_WIDTH);

            //
            // Compute convolution for loaded vectors and update 16 accumulators.
            // x0,x1 contain x(1:16) with x(1) at 32-bit LSW of x0 and x(16) at MSW
            // of x1; y0 contains y(1:8) with y(1) at 32-bit LSW and y(8) at MSW.
            // 8 80-bit complex accumulators w(1:8) are updated as follows:
            //   for n=1:8
            //     w(n) = w(n)+x(n+(8:-1:1)).*y(1:8);
            //   end
            //

            {
                xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31;

                xb_vecNx16 t0;
                xb_vecNx16 q0, q1, q2, q3;

                BBE_SELPCNX16I(p01, p00, x1, x0, 1);
                BBE_SELPCNX16I(p11, p10, x1, x0, 3);
                BBE_SELPCNX16I(p21, p20, x1, x0, 5);
                BBE_SELPCNX16I(p31, p30, x1, x0, 7);

                t0 = BBE_SHFLNX16I(y0, BBE_SHFLI_REVERSE_2);

                q0 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_0X4);
                q1 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_1X4);
                q2 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4);
                q3 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_3X4);

                BBE_MULANX16PC_0(w0, p00, q0);
                BBE_MULANX16PC_1(w0, p01, q0);
                BBE_MULANX16PC_0(w0, p10, q1);
                BBE_MULANX16PC_1(w0, p11, q1);
                BBE_MULANX16PC_0(w0, p20, q2);
                BBE_MULANX16PC_1(w0, p21, q2);
                BBE_MULANX16PC_0(w0, p30, q3);
                BBE_MULANX16PC_1(w0, p31, q3);
            }

            x0 = x1;
        }

        // CQ31 <- CQ30 + 1
        w0 = BBE_ADDNX40(w0, w0);

        // 32-bit saturation!
        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        // Save 8 convolution results, CQ31
        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);
    }

    //----------------------------------------------------------------------------
    // Part II. Compute r(M+1:N).

    // Right-hand data pointer will be reset to y(M-7) at each iteration of the
    // outer loop.
    Y = (const xb_vecNx16*)((uintptr_t)y + 4 * M - 2 * BBE_SIMD_WIDTH);

    for (n = 0; n<(N - M) / (BBE_SIMD_WIDTH / 2); n++)
    {
        // Load x(8*n+1:8)
        BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);

        // Latch the left-hand pointer and use the shuttle pointer to access the
        // left-hand sequence in the inner loop.
        S = X;

        w0 = 0;

        __Pragma("ymemory(S)");
        __Pragma("ymemory(Y)");
        __Pragma("loop_count min=1");
        for (m = 0; m<M / (BBE_SIMD_WIDTH / 2); m++)
        {
            // Load x(8*n+8*(m+1)+1:8)
            BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);
            // Load y(M-8*(m+1)+1:8)
            BBE_LVNX16_IP(y0, Y, -2 * BBE_SIMD_WIDTH);

            //
            // Compute convolution for loaded vectors and update 8 accumulators.
            //

            {
                xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31;

                xb_vecNx16 t0;
                xb_vecNx16 q0, q1, q2, q3;

                BBE_SELPCNX16I(p01, p00, x1, x0, 1);
                BBE_SELPCNX16I(p11, p10, x1, x0, 3);
                BBE_SELPCNX16I(p21, p20, x1, x0, 5);
                BBE_SELPCNX16I(p31, p30, x1, x0, 7);

                t0 = BBE_SHFLNX16I(y0, BBE_SHFLI_REVERSE_2);

                q0 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_0X4);
                q1 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_1X4);
                q2 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4);
                q3 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_3X4);

                BBE_MULANX16PC_0(w0, p00, q0);
                BBE_MULANX16PC_1(w0, p01, q0);
                BBE_MULANX16PC_0(w0, p10, q1);
                BBE_MULANX16PC_1(w0, p11, q1);
                BBE_MULANX16PC_0(w0, p20, q2);
                BBE_MULANX16PC_1(w0, p21, q2);
                BBE_MULANX16PC_0(w0, p30, q3);
                BBE_MULANX16PC_1(w0, p31, q3);
            }

            x0 = x1;
        }

        // Right-hand sequence is re-read (backwards) on each iteration of the
        // outer loop.
        Y = (const xb_vecNx16*)XT_ADDX4(M, (uintptr_t)Y);

        // CQ31 <- CQ30 + 1
        w0 = BBE_ADDNX40(w0, w0);

        // 32-bit saturation!
        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        // Save 8 convolution results, CQ31
        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);
    }

    //----------------------------------------------------------------------------
    // Part III. Compute r(N+1:N+M-1). Unlike first two parts, here we read the
    // left-hand data sequence x() in backward direction, and the right-hand
    // sequence y() - in forward direction.

    Y = (const xb_vecNx16*)y;

    for (n = M / (BBE_SIMD_WIDTH / 2) - 1; n >= 0; n--)
    {
        // Left-hand data pointer is reset to x(N-7) on each iteration.
        X = (const xb_vecNx16*)((uintptr_t)x + 4 * N - 2 * BBE_SIMD_WIDTH);

        x1 = 0;

        // Load y(M-8*(n+1)+1:8).
        BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
        // Latch the right-hand pointer and use the shuttle pointer to access the
        // right-hand sequence in the inner loop.
        S = Y;

        w0 = 0;

        __Pragma("ymemory(S)");
        __Pragma("ymemory(X)");
        for (m = 0; m<n + 1; m++)
        {
            // Load x(N-8*(m+1)+1:8)
            BBE_LVNX16_IP(x0, X, -2 * BBE_SIMD_WIDTH);

            //
            // Compute convolution for loaded vectors and update 8 accumulators.
            //

            {
                xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31;

                xb_vecNx16 t0;
                xb_vecNx16 q0, q1, q2, q3;

                BBE_SELPCNX16I(p01, p00, x1, x0, 1);
                BBE_SELPCNX16I(p11, p10, x1, x0, 3);
                BBE_SELPCNX16I(p21, p20, x1, x0, 5);
                BBE_SELPCNX16I(p31, p30, x1, x0, 7);

                t0 = BBE_SHFLNX16I(y0, BBE_SHFLI_REVERSE_2);

                q0 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_0X4);
                q1 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_1X4);
                q2 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_2X4);
                q3 = BBE_SHFLNX16I(t0, BBE_SHFLI_REP_3X4);

                BBE_MULANX16PC_0(w0, p00, q0);
                BBE_MULANX16PC_1(w0, p01, q0);
                BBE_MULANX16PC_0(w0, p10, q1);
                BBE_MULANX16PC_1(w0, p11, q1);
                BBE_MULANX16PC_0(w0, p20, q2);
                BBE_MULANX16PC_1(w0, p21, q2);
                BBE_MULANX16PC_0(w0, p30, q3);
                BBE_MULANX16PC_1(w0, p31, q3);
            }

            x1 = x0;

            // Load y(M-8*n+8*m+1:8)
            BBE_LVNX16_IP(y0, S, 2 * BBE_SIMD_WIDTH);
        }

        // CQ31 <- CQ30 + 1
        w0 = BBE_ADDNX40(w0, w0);

        // 32-bit saturation!
        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        ix = XT_MIN(1, n);

        // Mask index is 1 for all iteration but the last.
        vb0 = BBE_LBN_I((const vboolN*)cst + ix, 0);

        // Save 8 or 7 convolution results, CQ31. 7 results are saved on the last
        // iteration.
        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16T_IP(r1, R, 2 * BBE_SIMD_WIDTH, vb0);
    }
} /* fir_convol_lin() */

#endif

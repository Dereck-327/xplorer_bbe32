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
    Autocorrelation for a Ñomplex Data Vector
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
Autocorrelation for a Ñomplex Data Vector

Estimates the auto-correlation of complex-valued vector x, positive side
only. Returns autocorrelation of length N. For an input vector of N complex
samples x[0..N-1] the computation follows the MATLAB code given below:

  for n = 1:N
    r(n) = sum(x(n:N).*conj(x(1:(N+1-n))));
  end

Representation:
fir_acorr   Signed fixed-point format
            Input data are 16-bit Q15, output data are 32-bit Q31
fir_acorrf  IEEE-754 Std. single precision floating-point format for
            input/output data

Parameters:
Input:
x[N]        Complex input data
N           Length of x
Output:
r[N]        Complex output data

Restrictions:
x,r         Must not overlap
x,r         Must be aligned on 32-byte boundary
N           Multiple of 8 (fir_acorr) or 4 (fir_acorrf)
-------------------------------------------------------------------------*/

#if !HAVE_MULPC
DISCARD_FUN(void, fir_acorr, (complex_fract32 * restrict r, const complex_fract16 * restrict x, int N))
#else

void fir_acorr(complex_fract32 * restrict r, const complex_fract16 * restrict x, int N)
{
          xb_vecNx16 * restrict R;
    const xb_vecNx16 *          X0; // Left-hand sequence pointer
    const xb_vecNx16 *          X1; // Shuttle pointer for the left-hand sequence
    const xb_vecNx16 *          Y;  // Right-hand sequence pointer
    int n, m, M;

    NASSERT(!(N & 7));
    NASSERT(r && x);
    NASSERT_ALIGN32(r);
    NASSERT_ALIGN32(x);

    // Inner loop iterations number decreases on each iteration of the outer loop,
    // because the overlapping of left-hand and right-hand data sequences shortens
    // when the index of the autocorrelation sequence increases. Initial value is
    // also decreased by 1 because the last SIMD_WIDTH/2 accumulator updates are
    // pinched off the inner loop.
    M = N / (BBE_SIMD_WIDTH / 2) - 1;

    R = (xb_vecNx16*)r;
    X0 = (const xb_vecNx16*)x;

    //
    // We compute SIMD_WIDTH/2 autocorrelation values at a time.
    //

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++, M--)
    {
        xb_vecNx16 x0, x1, y0, r0, r1;
        xb_vecNx40 w0;

        w0 = 0;

        // Load x[n*SIMD_WIDTH/2+0..SIMD_WIDTH/2-1], CQ15.
        BBE_LVNX16_IP(x0, X0, 2 * BBE_SIMD_WIDTH);
        // Latch the left-hand sequence pointer to reuse it on the next iteration of
        // the outer loop. The inner loop reads the rest of sequence through the
        // shuttle pointer.
        X1 = X0;

        // Reset the right-hand sequence pointer to x[0].
        Y = (const xb_vecNx16*)x;

        __Pragma("ymemory(X1)");
        __Pragma("ymemory(Y)");
        for (m = 0; m<M; m++)
        {
            // Load x[(n+m+1)*SIMD_WIDTH/2+0..SIMD_WIDTH/2-1], CQ15
            BBE_LVNX16_IP(x1, X1, 2 * BBE_SIMD_WIDTH);
            // Load x[m*SIMD_WIDTH/2+0..SIMD_WIDTH/2-1], CQ15
            BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);

            // Update w0, w1 with correlation between x0 and and y0.
            // w0 contains accumulators for even-numbered autocorrelation values
            // (CQ30) at even complex positions. w1 - odd-numbered values at odd
            // positions. We could use a single wvec, but splitting accumulators
            // into two wvecs avoids pipeline stalls. 
            {
                xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31;

                xb_vecNx16 q0, q1, q2, q3;

                BBE_SELPCNX16I(p01, p00, x1, x0, 0);
                BBE_SELPCNX16I(p11, p10, x1, x0, 2);
                BBE_SELPCNX16I(p21, p20, x1, x0, 4);
                BBE_SELPCNX16I(p31, p30, x1, x0, 6);

                y0 = BBE_CONJSNX16C(y0);

                q0 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_0X4);
                q1 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_1X4);
                q2 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_2X4);
                q3 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_3X4);

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

        {
            // To handle the last SIMD_WIDTH/2 values of the left-hand sequence, we
            // have to pad it with SIMD_WIDTH/2 zeros.
            x1 = 0;

            // Load x[N-(n+1)*SIMD_WIDTH/2+0..SIMD_WIDTH/2-1], CQ15
            BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);

            {
                xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31;

                xb_vecNx16 q0, q1, q2, q3;

                BBE_SELPCNX16I(p01, p00, x1, x0, 0);
                BBE_SELPCNX16I(p11, p10, x1, x0, 2);
                BBE_SELPCNX16I(p21, p20, x1, x0, 4);
                BBE_SELPCNX16I(p31, p30, x1, x0, 6);

                y0 = BBE_CONJSNX16C(y0);

                q0 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_0X4);
                q1 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_1X4);
                q2 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_2X4);
                q3 = BBE_SHFLNX16I(y0, BBE_SHFLI_REP_3X4);

                BBE_MULANX16PC_0(w0, p00, q0);
                BBE_MULANX16PC_1(w0, p01, q0);
                BBE_MULANX16PC_0(w0, p10, q1);
                BBE_MULANX16PC_1(w0, p11, q1);
                BBE_MULANX16PC_0(w0, p20, q2);
                BBE_MULANX16PC_1(w0, p21, q2);
                BBE_MULANX16PC_0(w0, p30, q3);
                BBE_MULANX16PC_1(w0, p31, q3);
            }
        }

        // CQ31 <- CQ30 + CQ30
        w0 = BBE_ADDNX40(w0, w0);

        // 32-bit signed saturation.
        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        // Store SIMD_WIDTH/2 32-bit autocorrelation values to the output array.
        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);
    }
} /* fir_acorr() */

#endif

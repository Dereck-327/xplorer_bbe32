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

#if !HAVE_MULPC
DISCARD_FUN(void, fir_xcorr_lin, (complex_fract32 * restrict r,
                            const complex_fract16 * restrict x,
                            const complex_fract16 * restrict y,
                            int N, int M))
#else

#define FIR_XCORR_NXN( w0, x1, x0, y0, ofs )           \
  {                                                    \
    xb_vecNx16 p00, p01, p10, p11, p20, p21, p30, p31; \
                                                       \
    xb_vecNx16 q0, q1, q2, q3;                         \
                                                       \
    BBE_SELPCNX16I( p01, p00, x1, x0,  0 + ofs );      \
    BBE_SELPCNX16I( p11, p10, x1, x0,  2 + ofs );      \
    BBE_SELPCNX16I( p21, p20, x1, x0,  4 + ofs );      \
    BBE_SELPCNX16I( p31, p30, x1, x0,  6 + ofs );      \
                                                       \
    y0 = BBE_CONJSNX16C( y0 );                         \
                                                       \
    q0 = BBE_SHFLNX16I( y0, BBE_SHFLI_REP_0X4 );       \
    q1 = BBE_SHFLNX16I( y0, BBE_SHFLI_REP_1X4 );       \
    q2 = BBE_SHFLNX16I( y0, BBE_SHFLI_REP_2X4 );       \
    q3 = BBE_SHFLNX16I( y0, BBE_SHFLI_REP_3X4 );       \
                                                       \
    BBE_MULANX16PC_0( w0, p00, q0 );                   \
    BBE_MULANX16PC_1( w0, p01, q0 );                   \
    BBE_MULANX16PC_0( w0, p10, q1 );                   \
    BBE_MULANX16PC_1( w0, p11, q1 );                   \
    BBE_MULANX16PC_0( w0, p20, q2 );                   \
    BBE_MULANX16PC_1( w0, p21, q2 );                   \
    BBE_MULANX16PC_0( w0, p30, q3 );                   \
    BBE_MULANX16PC_1( w0, p31, q3 );                   \
  }

void fir_xcorr_lin ( complex_fract32 * restrict r,
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

    Y = (const xb_vecNx16*)((uintptr_t)y + 4 * M - 2 * BBE_SIMD_WIDTH);

    for (n = 0; n<M / (BBE_SIMD_WIDTH / 2); n++)
    {
        S = X;

        x0 = 0;

        w0 = 0;

        __Pragma("ymemory(S)");
        __Pragma("ymemory(Y)");
        for (m = 0; m<n + 1; m++)
        {
            BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);

            FIR_XCORR_NXN(w0, x1, x0, y0, 1);

            x0 = x1;
        }

        Y = (const xb_vecNx16*)((uintptr_t)Y - (n + 2) * 2 * BBE_SIMD_WIDTH);

        w0 = BBE_ADDNX40(w0, w0);

        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);
    }

    //----------------------------------------------------------------------------
    // Part II. Compute r(M+1:N).

    Y = (const xb_vecNx16*)y;

    for (n = 0; n<(N - M) / (BBE_SIMD_WIDTH / 2); n++)
    {
        BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);

        S = X;

        w0 = 0;

        __Pragma("ymemory(S)");
        __Pragma("ymemory(Y)");
        __Pragma("loop_count min=1");
        for (m = 0; m<M / (BBE_SIMD_WIDTH / 2); m++)
        {
            BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);

            FIR_XCORR_NXN(w0, x1, x0, y0, 1);

            x0 = x1;
        }

        Y = (const xb_vecNx16*)((uintptr_t)Y - 4 * M);

        w0 = BBE_ADDNX40(w0, w0);

        r1 = BBE_MOVSVWH(w0);
        r0 = BBE_MOVSVWL(w0);

        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);
    }

    //----------------------------------------------------------------------------
    // Part III. Compute r(N+1:N+M-1). Unlike first two parts, here we read the
    // data sequences in backward direction.

    Y = (const xb_vecNx16*)((uintptr_t)y + 4 * M - 2 * BBE_SIMD_WIDTH);

    for (n = M / (BBE_SIMD_WIDTH / 2) - 1; n >= 0; n--)
    {
        X = (const xb_vecNx16*)((uintptr_t)x + 4 * N - 2 * BBE_SIMD_WIDTH);

        x1 = 0;

        S = Y;

        Y = (const xb_vecNx16*)((uintptr_t)Y - 2 * BBE_SIMD_WIDTH);

        w0 = 0;

        __Pragma("ymemory(S)");
        __Pragma("ymemory(X)");
        for (m = 0; m<n + 1; m++)
        {
            BBE_LVNX16_IP(x0, X, -2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(y0, S, -2 * BBE_SIMD_WIDTH);

            FIR_XCORR_NXN(w0, x1, x0, y0, 1);

            x1 = x0;
        }

        w0 = BBE_ADDNX40(w0, w0);

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
} /* fir_xcorr_lin() */

#endif

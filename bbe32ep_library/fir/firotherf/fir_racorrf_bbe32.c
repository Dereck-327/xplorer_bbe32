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
    Autocorrelation for a Real Data Vector
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
Autocorrelation for a Real Data Vector

Estimates the auto-correlation of real-valued vector x, positive side 
only. Returns autocorrelation of length N. For an input vector of N 
samples x[0..N-1] the computation follows the MATLAB code given below:

  for n = 1:N
    r(n) = sum(x(n:N).*x(1:(N+1-n)));
  end

Representation: IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]         Input data
N            Length of x
Output:
r[N]         Output data

Restrictions:
x,r          Must not overlap
x,r          Must be aligned on 32-byte boundary
N            Multiple of 8 (fir_racorrf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,fir_racorrf,(float32_t * restrict r, const float32_t * restrict x, int N))
#else
void fir_racorrf(float32_t * restrict r, const float32_t * restrict x, int N)
{
    int n, m;
    int size;
    const xb_vecN_2xf32 * restrict pXnm;
    const xb_vecN_2xf32 * restrict pXm;
          xb_vecN_2xf32 * restrict pR;
    xb_vecN_2xf32 x0, x1, x2, x3, x4, x5, x6, x7, xm, xnm_pre, xnm_cur;
    xb_vecN_2xf32 xn0, xn1, xn2, xn3, xn4, xn5, xn6, xn7;
    valign vx;

    ASSERT(x && r && N >= 0);
    NASSERT_ALIGN(r, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N >= 0 && N % (BBE_SIMD_WIDTH / 2) == 0);

    pR = (xb_vecN_2xf32 *)r;

    for (n = 0; n < N; n += BBE_SIMD_WIDTH / 2)
    {
        xb_vecN_2xf32 A = 0.f;
        size = (N - n) * sizeof(float32_t);
        pXm = (const xb_vecN_2xf32 *)x;
        pXnm = (const xb_vecN_2xf32 *)(x + n);
        vx = BBE_LAN_2XF32_PP(pXnm);
        BBE_LAVN_2XF32_XP(xnm_pre, vx, pXnm, size);

        for (m = 0; m < N - n; m += BBE_SIMD_WIDTH / 2)
        {
            BBE_LVN_2XF32_IP(xm, pXm, 2 * BBE_SIMD_WIDTH);
            x0 = BBE_REPN_2XF32(xm, 0);
            x1 = BBE_REPN_2XF32(xm, 1);
            x2 = BBE_REPN_2XF32(xm, 2);
            x3 = BBE_REPN_2XF32(xm, 3);
            x4 = BBE_REPN_2XF32(xm, 4);
            x5 = BBE_REPN_2XF32(xm, 5);
            x6 = BBE_REPN_2XF32(xm, 6);
            x7 = BBE_REPN_2XF32(xm, 7);

            size -= 2 * BBE_SIMD_WIDTH;
            BBE_LAVN_2XF32_XP(xnm_cur, vx, pXnm, size);
            xn0 = xnm_pre;
            xn1 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_2);
            xn2 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_4);
            xn3 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_6);
            xn4 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_8);
            xn5 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_10);
            xn6 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_12);
            xn7 = BBE_SELN_2XF32I(xnm_cur, xnm_pre, BBE_SELI_ROTATE_RIGHT_14);
            xnm_pre = xnm_cur;

            BBE_MULAN_2XF32(A, x0, xn0);
            BBE_MULAN_2XF32(A, x1, xn1);
            BBE_MULAN_2XF32(A, x2, xn2);
            BBE_MULAN_2XF32(A, x3, xn3);
            BBE_MULAN_2XF32(A, x4, xn4);
            BBE_MULAN_2XF32(A, x5, xn5);
            BBE_MULAN_2XF32(A, x6, xn6);
            BBE_MULAN_2XF32(A, x7, xn7);
        }
        BBE_SVN_2XF32_IP(A, pR, 2 * BBE_SIMD_WIDTH);
    }
} //fir_racorrf()
#endif

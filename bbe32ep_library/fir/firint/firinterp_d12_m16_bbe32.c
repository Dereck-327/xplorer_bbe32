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
    Interpolating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

#include "firinterp_common.h"

/*-------------------------------------------------------------------------
Interpolating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with interpolation using real
IR stored in vector h. The complex data input is stored in vector x. The
filter output result is stored in vector y. The filter calculates N*D complex
output samples using M*D coefficients and requires last N+M-1 samples in the
delay line.

Representation:
firinterp   16-bit signed fixed-point format
            Filter coefficients are Q15
            Number of fractional bits for input/output samples is user-difined
firinterpf  IEEE-754 Std. single precision floating-point format for filter 
            coefficients and input/output samples

Parameters:
Input:
D           Interpolation ratio 
N           Length of input sample block
M           Length of subfilter. Total length of filter is M*D
h[M*D]      Filter coefficients; h[0] is to be multiplied by the newest 
            sample,Q15
x[N]        Input complex samples
Output:
y[N*D]      Output complex samples

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
N           Multiple of 8 (firinterp) or 4 (firinterpf)
M           2,4,8 or a positive multiple of 16 for D=2,3,4,6,12; or 
            a positive multiple of 8 for other D
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
subfilter lengths M=2,4,8,16 and 32 and interpolation factors D=2,3 and 4,
in any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firinterp[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Filter processing function for D=12 M=8 and N=8*n. */
static void filter_proc_12d_16_8n (   
                              void* handle,
                                 int16_t *  restrict  y,
                           const int16_t *  restrict  x,
                           const int16_t *  restrict  coef,
                                 int16_t *  restrict  delayLine,
                                      int   M,
                                      int   N,
                                      int   D
                          )
{
    int n;
    int16_t k;
    const xb_vecNx16 *  restrict pX = (const xb_vecNx16 *)x;
    const        int *  restrict pH = (const int        *)coef;
    xb_vecNx16 *  restrict pD = (xb_vecNx16 *)delayLine;
    xb_vecNx16 *  restrict pY = (xb_vecNx16 *)y;

    xb_vecNx16 x0, x1, y0, y1;
    xb_vecNx16 d0, d1, d2;
    uint32_t   c00, c01, c02, c03, c10, c11, c12, c13;
    uint32_t   c04, c05, c06, c07, c14, c15, c16, c17;

    xb_vecNx40 A0, A1;
    vsaN       shft;
    NASSERT(N>0 && N % 8 == 0);
    NASSERT(M == 16);
    NASSERT(D == 12);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);
    shft = BBE_MOVVSA32(11);

    for (k = 0; k<6; k++)
    {
        d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
        d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
        pX = (const xb_vecNx16 *)x;
        pH = (const int        *)(coef + 4 * k);
        pY = (xb_vecNx16 *)(y + 2 * BBE_SIMD_WIDTH*k);
        for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
        {
            c00 = XT_L32I(pH, 0);
            c10 = XT_L32I(pH, 4);
            c01 = XT_L32I(pH, 48);
            c11 = XT_L32I(pH, 52);
            c02 = XT_L32I(pH, 96);
            c12 = XT_L32I(pH, 100);
            c03 = XT_L32I(pH, 144);
            c13 = XT_L32I(pH, 148);
            c04 = XT_L32I(pH, 192);
            c14 = XT_L32I(pH, 196);
            c05 = XT_L32I(pH, 240);
            c15 = XT_L32I(pH, 244);
            c06 = XT_L32I(pH, 288);
            c16 = XT_L32I(pH, 292);
            c07 = XT_L32I(pH, 336);
            c17 = XT_L32I(pH, 340);

            BBE_LVNX16_IP(d2, pX, 2 * BBE_SIMD_WIDTH);
            BBE_SELPCNX16I(x1, x0, d2, d1, 7);
            A0 = BBE_MULNX16PR(x0, x1, c00);
            A1 = BBE_MULNX16PR(x0, x1, c10);
            BBE_SELPCNX16I(x1, x0, d2, d1, 5);
            BBE_MULANX16PR(A0, x0, x1, c01);
            BBE_MULANX16PR(A1, x0, x1, c11);
            BBE_SELPCNX16I(x1, x0, d2, d1, 3);
            BBE_MULANX16PR(A0, x0, x1, c02);
            BBE_MULANX16PR(A1, x0, x1, c12);
            BBE_SELPCNX16I(x1, x0, d2, d1, 1);
            BBE_MULANX16PR(A0, x0, x1, c03);
            BBE_MULANX16PR(A1, x0, x1, c13);
            BBE_SELPCNX16I(x1, x0, d1, d0, 7);
            BBE_MULANX16PR(A0, x0, x1, c04);
            BBE_MULANX16PR(A1, x0, x1, c14);
            BBE_SELPCNX16I(x1, x0, d1, d0, 5);
            BBE_MULANX16PR(A0, x0, x1, c05);
            BBE_MULANX16PR(A1, x0, x1, c15);
            BBE_SELPCNX16I(x1, x0, d1, d0, 3);
            BBE_MULANX16PR(A0, x0, x1, c06);
            BBE_MULANX16PR(A1, x0, x1, c16);
            BBE_SELPCNX16I(x1, x0, d1, d0, 1);
            BBE_MULANX16PR(A0, x0, x1, c07);
            BBE_MULANX16PR(A1, x0, x1, c17);

            y0 = BBE_PACKVNX40(A0, shft);
            y1 = BBE_PACKVNX40(A1, shft);

            BBE_SVNX16_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(y1, pY, 11 * 2 * BBE_SIMD_WIDTH);

            d0 = d1;
            d1 = d2;
        }
    }

  firinterp_12d_intlv(y, N);

  BBE_SVNX16_IP(d0,pD,2*BBE_SIMD_WIDTH);
  BBE_SVNX16_IP(d1,pD,2*BBE_SIMD_WIDTH);
};
const tFirFxdxns interp_12d_16_8n   ={&firinterp_dx,filter_proc_12d_16_8n };


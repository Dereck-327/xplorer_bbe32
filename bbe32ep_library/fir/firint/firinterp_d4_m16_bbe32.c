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

/* Filter processing function for D=4 M==16 and N%8==0 */
static void filter_proc_4d_16_8n (    
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
    const xb_vecNx16 *  restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 *  restrict pH = (const xb_vecNx16 *)coef;
    xb_vecNx16 *  restrict pD = (xb_vecNx16 *)delayLine;
    xb_vecNx16 *  restrict pY = (xb_vecNx16 *)y;

    xb_vecNx16 y0, y1, y2, y3;
    xb_vecNx16 d0, d1, d2;
    xb_vecNx16 t0, t1, t2, t3;
    uint32_t   c00, c01, c02, c03, c10, c11, c12, c13;
    uint32_t   c04, c05, c06, c07, c14, c15, c16, c17;
    uint32_t   c20, c21, c22, c23, c24, c25, c26, c27;
    uint32_t   c30, c31, c32, c33, c34, c35, c36, c37;
    xb_vecNx40 A0, A1, A2, A3;
    vsaN       shft;
    xb_vecNx16 CoefVec;
    NASSERT(N % 8 == 0);
    NASSERT(M == 16);
    NASSERT(D == 4);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    shft = BBE_MOVVSA32(13);
    if (N<=0) return;
    __Pragma("loop_count min=1");
    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        xb_vecNx16 p00, p01, p02, p03, p04, p05, p06, p07;
        xb_vecNx16 p10, p11, p12, p13, p14, p15, p16, p17;

        CoefVec = BBE_LVNX16_I(pH, 0);
        c00 = BBE_EXTRNX16C(CoefVec, 0);
        c10 = BBE_EXTRNX16C(CoefVec, 1);
        c20 = BBE_EXTRNX16C(CoefVec, 2);
        c30 = BBE_EXTRNX16C(CoefVec, 3);
        c01 = BBE_EXTRNX16C(CoefVec, 4);
        c11 = BBE_EXTRNX16C(CoefVec, 5);
        c21 = BBE_EXTRNX16C(CoefVec, 6);
        c31 = BBE_EXTRNX16C(CoefVec, 7);

        CoefVec = BBE_LVNX16_I(pH, 2 * BBE_SIMD_WIDTH);
        c02 = BBE_EXTRNX16C(CoefVec, 0);
        c12 = BBE_EXTRNX16C(CoefVec, 1);
        c22 = BBE_EXTRNX16C(CoefVec, 2);
        c32 = BBE_EXTRNX16C(CoefVec, 3);
        c03 = BBE_EXTRNX16C(CoefVec, 4);
        c13 = BBE_EXTRNX16C(CoefVec, 5);
        c23 = BBE_EXTRNX16C(CoefVec, 6);
        c33 = BBE_EXTRNX16C(CoefVec, 7);

        CoefVec = BBE_LVNX16_I(pH, 2 * 2 * BBE_SIMD_WIDTH);
        c04 = BBE_EXTRNX16C(CoefVec, 0);
        c14 = BBE_EXTRNX16C(CoefVec, 1);
        c24 = BBE_EXTRNX16C(CoefVec, 2);
        c34 = BBE_EXTRNX16C(CoefVec, 3);
        c05 = BBE_EXTRNX16C(CoefVec, 4);
        c15 = BBE_EXTRNX16C(CoefVec, 5);
        c25 = BBE_EXTRNX16C(CoefVec, 6);
        c35 = BBE_EXTRNX16C(CoefVec, 7);

        CoefVec = BBE_LVNX16_I(pH, 3 * 2 * BBE_SIMD_WIDTH);
        c06 = BBE_EXTRNX16C(CoefVec, 0);
        c16 = BBE_EXTRNX16C(CoefVec, 1);
        c26 = BBE_EXTRNX16C(CoefVec, 2);
        c36 = BBE_EXTRNX16C(CoefVec, 3);
        c07 = BBE_EXTRNX16C(CoefVec, 4);
        c17 = BBE_EXTRNX16C(CoefVec, 5);
        c27 = BBE_EXTRNX16C(CoefVec, 6);
        c37 = BBE_EXTRNX16C(CoefVec, 7);

        d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
        d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(d2, pX, 2 * BBE_SIMD_WIDTH);

        BBE_SELPCNX16I(p01, p00, d1, d0, 1);
        BBE_SELPCNX16I(p03, p02, d1, d0, 3);
        BBE_SELPCNX16I(p05, p04, d1, d0, 5);
        BBE_SELPCNX16I(p07, p06, d1, d0, 7);

        BBE_SELPCNX16I(p11, p10, d2, d1, 1);
        BBE_SELPCNX16I(p13, p12, d2, d1, 3);
        BBE_SELPCNX16I(p15, p14, d2, d1, 5);
        BBE_SELPCNX16I(p17, p16, d2, d1, 7);

        BBE_SVNX16_I(d1, pD, 0 * 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_I(d2, pD, 1 * 2 * BBE_SIMD_WIDTH);

        A0 = BBE_MULNX16PR(p06, p07, c04);
        BBE_MULANX16PR(A0, p04, p05, c05);
        BBE_MULANX16PR(A0, p02, p03, c06);
        BBE_MULANX16PR(A0, p00, p01, c07);
        BBE_MULANX16PR(A0, p16, p17, c00);
        BBE_MULANX16PR(A0, p14, p15, c01);
        BBE_MULANX16PR(A0, p12, p13, c02);
        BBE_MULANX16PR(A0, p10, p11, c03);

        A1 = BBE_MULNX16PR(p06, p07, c14);
        BBE_MULANX16PR(A1, p04, p05, c15);
        BBE_MULANX16PR(A1, p02, p03, c16);
        BBE_MULANX16PR(A1, p00, p01, c17);
        BBE_MULANX16PR(A1, p16, p17, c10);
        BBE_MULANX16PR(A1, p14, p15, c11);
        BBE_MULANX16PR(A1, p12, p13, c12);
        BBE_MULANX16PR(A1, p10, p11, c13);

        A2 = BBE_MULNX16PR(p06, p07, c24);
        BBE_MULANX16PR(A2, p04, p05, c25);
        BBE_MULANX16PR(A2, p02, p03, c26);
        BBE_MULANX16PR(A2, p00, p01, c27);
        BBE_MULANX16PR(A2, p16, p17, c20);
        BBE_MULANX16PR(A2, p14, p15, c21);
        BBE_MULANX16PR(A2, p12, p13, c22);
        BBE_MULANX16PR(A2, p10, p11, c23);

        A3 = BBE_MULNX16PR(p06, p07, c34);
        BBE_MULANX16PR(A3, p04, p05, c35);
        BBE_MULANX16PR(A3, p02, p03, c36);
        BBE_MULANX16PR(A3, p00, p01, c37);
        BBE_MULANX16PR(A3, p16, p17, c30);
        BBE_MULANX16PR(A3, p14, p15, c31);
        BBE_MULANX16PR(A3, p12, p13, c32);
        BBE_MULANX16PR(A3, p10, p11, c33);

        y0 = BBE_PACKVNX40(A0, shft);
        y1 = BBE_PACKVNX40(A1, shft);
        y2 = BBE_PACKVNX40(A2, shft);
        y3 = BBE_PACKVNX40(A3, shft);
        BBE_DSELNX16I(t2, t0, y1, y0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(t3, t1, y3, y2, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(y1, y0, t1, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(y3, y2, t3, t2, BBE_DSELI_INTERLEAVE_4);

        BBE_SVNX16_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(y1, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(y2, pY, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(y3, pY, 2 * BBE_SIMD_WIDTH);
    }
}
const tFirFxdxns interp_4d_16_8n   ={&firinterp_dx,filter_proc_4d_16_8n} ;

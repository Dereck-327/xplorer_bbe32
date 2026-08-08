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

/* Filter processing function for D=4 M==32 and N%8==0 */
static void filter_proc_4d_32_8n (    
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
    int n, k;
    const xb_vecNx16 *  restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 *  restrict pH = (const xb_vecNx16 *)coef;
    xb_vecNx16 *  restrict pD = (xb_vecNx16 *)delayLine;
    xb_vecNx16 *  restrict pY = (xb_vecNx16 *)y;
    xb_vecNx16 *  restrict pY_ = (xb_vecNx16 *)y;

    xb_vecNx16 y0, y1, y2, y3;
    xb_vecNx16 t0, t1, t2, t3;
    xb_vecNx16 d0, d1, d2, d3, d4;
    uint32_t   c00, c01, c02, c03;
    uint32_t   c04, c05, c06, c07;
    uint32_t   c08, c09, c0a, c0b;
    uint32_t   c0c, c0d, c0e, c0f;
    xb_vecNx40 A0;
    vsaN       shft;
    NASSERT(N>0 && N % 8 == 0);
    NASSERT(M == 32);
    NASSERT(D == 4);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    shft = BBE_MOVVSA32(13);
    for (k = 0; k < D; k++)
    {
        pX = (const xb_vecNx16 *)x;
        pY = (xb_vecNx16 *)(y + k*BBE_SIMD_WIDTH);
        d0 = BBE_LVNX16_I(pD, 0 * 2 * BBE_SIMD_WIDTH);
        d1 = BBE_LVNX16_I(pD, 1 * 2 * BBE_SIMD_WIDTH);
        d2 = BBE_LVNX16_I(pD, 2 * 2 * BBE_SIMD_WIDTH);
        d3 = BBE_LVNX16_I(pD, 3 * 2 * BBE_SIMD_WIDTH);

        for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
        {
            xb_vecNx16 p00,p01,p02,p03,p04,p05,p06,p07;
            xb_vecNx16 p08,p09,p0a,p0b,p0c,p0d,p0e,p0f;
            xb_vecNx16 p10,p11,p12,p13,p14,p15,p16,p17;
            xb_vecNx16 p18,p19,p1a,p1b,p1c,p1d,p1e,p1f;

            {
                xb_vecNx16 CoefVec;
                CoefVec = BBE_LVNX16_I(pH,0*2*BBE_SIMD_WIDTH);
                c00 = BBE_EXTRNX16C(CoefVec,0);
                c01 = BBE_EXTRNX16C(CoefVec,1); 
                c02 = BBE_EXTRNX16C(CoefVec,2); 
                c03 = BBE_EXTRNX16C(CoefVec,3); 
                c04 = BBE_EXTRNX16C(CoefVec,4); 
                c05 = BBE_EXTRNX16C(CoefVec,5); 
                c06 = BBE_EXTRNX16C(CoefVec,6); 
                c07 = BBE_EXTRNX16C(CoefVec,7);

                CoefVec = BBE_LVNX16_I(pH,2*BBE_SIMD_WIDTH);
                c08 = BBE_EXTRNX16C(CoefVec,0);
                c09 = BBE_EXTRNX16C(CoefVec,1); 
                c0a = BBE_EXTRNX16C(CoefVec,2); 
                c0b = BBE_EXTRNX16C(CoefVec,3); 
                c0c = BBE_EXTRNX16C(CoefVec,4); 
                c0d = BBE_EXTRNX16C(CoefVec,5); 
                c0e = BBE_EXTRNX16C(CoefVec,6); 
                c0f = BBE_EXTRNX16C(CoefVec,7);

            }
            BBE_LVNX16_IP(d4,pX,2*BBE_SIMD_WIDTH);  

            BBE_SELPCNX16I(p01,p00,d1,d0,1);
            BBE_SELPCNX16I(p03,p02,d1,d0,3);
            BBE_SELPCNX16I(p05,p04,d1,d0,5);
            BBE_SELPCNX16I(p07,p06,d1,d0,7);

            BBE_SELPCNX16I(p09,p08,d2,d1,1);
            BBE_SELPCNX16I(p0b,p0a,d2,d1,3);
            BBE_SELPCNX16I(p0d,p0c,d2,d1,5);
            BBE_SELPCNX16I(p0f,p0e,d2,d1,7);

            BBE_SELPCNX16I(p11,p10,d3,d2,1);
            BBE_SELPCNX16I(p13,p12,d3,d2,3);
            BBE_SELPCNX16I(p15,p14,d3,d2,5);
            BBE_SELPCNX16I(p17,p16,d3,d2,7);

            BBE_SELPCNX16I(p19,p18,d4,d3,1);
            BBE_SELPCNX16I(p1b,p1a,d4,d3,3);
            BBE_SELPCNX16I(p1d,p1c,d4,d3,5);
            BBE_SELPCNX16I(p1f,p1e,d4,d3,7);

            A0 = BBE_MULNX16PR(p0e,p0f,c08);
            BBE_MULANX16PR(A0,p0c,p0d,c09);
            BBE_MULANX16PR(A0,p0a,p0b,c0a);
            BBE_MULANX16PR(A0,p08,p09,c0b);
            BBE_MULANX16PR(A0,p06,p07,c0c);
            BBE_MULANX16PR(A0,p04,p05,c0d);
            BBE_MULANX16PR(A0,p02,p03,c0e);
            BBE_MULANX16PR(A0,p00,p01,c0f);

            BBE_MULANX16PR(A0,p1e,p1f,c00); 
            BBE_MULANX16PR(A0,p1c,p1d,c01);
            BBE_MULANX16PR(A0,p1a,p1b,c02);
            BBE_MULANX16PR(A0,p18,p19,c03);
            BBE_MULANX16PR(A0,p16,p17,c04);
            BBE_MULANX16PR(A0,p14,p15,c05);
            BBE_MULANX16PR(A0,p12,p13,c06);
            BBE_MULANX16PR(A0,p10,p11,c07);

            d0 = d1;
            d1 = d2;
            d2 = d3;
            d3 = BBE_LVNX16_I(pX,-2*BBE_SIMD_WIDTH);
            y0 = BBE_PACKVNX40(A0,shft);
            BBE_SVNX16_XP(y0,pY,2*D*BBE_SIMD_WIDTH);
        }
        pH = pH + 2;
    }
    BBE_SVNX16_I(d0,pD,0*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(d1,pD,1*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(d2,pD,2*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(d3,pD,3*2*BBE_SIMD_WIDTH);
    pY = (      xb_vecNx16 *) y;
    pY_ = (      xb_vecNx16 *) y;
    __Pragma("no_reorder")
    for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
    {
        BBE_LVNX16_IP(y0,pY_,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(y1,pY_,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(y2,pY_,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(y3,pY_,2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(t2,t0,y1,y0,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(t3,t1,y3,y2,BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(y1,y0,t1,t0,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(y3,y2,t3,t2,BBE_DSELI_INTERLEAVE_4);

        BBE_SVNX16_IP(y0,pY,2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP(y1,pY,2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP(y2,pY,2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP(y3,pY,2*BBE_SIMD_WIDTH );
    }
}
const tFirFxdxns interp_4d_32_8n = { &firinterp_dx_m_32, filter_proc_4d_32_8n };

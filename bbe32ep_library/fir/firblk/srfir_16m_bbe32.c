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
    Block Real FIR filter w/ Symmetric Impulse Response
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "srfir_common.h"

/*-------------------------------------------------------------------------
Block Real FIR filter w/ Symmetric Impulse Response

Passes real input data through a direct-form FIR filter with real coefficients
and symmetric impulse response. The filter calculates N output samples using
(M+1)/2 coefficients and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by srfir_algDelay(M)

Precision: 16-bit data, 16-bit coefficients, 16-bit outputs, everything in Q15

Parameters:
Input:
objmem      Allocated memory block
h[(M+1)/2]  Filter coefficients (half the number of filter taps); h[0] is to
            be multiplied by the newest sample,  Q15
N           Length of sample block
M           Length of filter
x[N]        Input samples, Q15
Output:
y[N]        Output samples, Q15

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
M           2,4,8 or a multiple of 16
N           Multiple of 16  
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8,16 and 32.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srfir_init returns NULL handle.
---------------------------------------------------------------------------*/


#if !(HAVE_ABCD && HAVE_ADDPN && 1)
DISCARD_FUN(void, srfir_process_16m,  (void* _srfir, int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N))
#else

#define MIN(a,b)  ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)  ( (a)>(b) ? (a) : (b) )

/*    M%16==0     */
void srfir_process_16m (void* _srfir, int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N)
{
    static const int16_t ALIGN(32) ss[2][BBE_SIMD_WIDTH] = {
        { 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7 },
        { 9, 10, 11, 12, 13, 14, 15, 16, 8, 9, 10, 11, 12, 13, 14, 15 }
    };
    srfir_t *srfir = (srfir_t *)_srfir;
    xb_vecNx16* restrict W0;
    xb_vecNx16* restrict W1;
    xb_vecNx16* restrict R0;
    xb_vecNx16* restrict R1;
    xb_vecNx16* restrict Y;
    const xb_vecNx16* restrict pX;
    const xb_vecNx16* restrict H;
    const xb_vecNx16* restrict SEL = (const xb_vecNx16*)ss;

    xb_vecNx40 A;
    xb_vecNx16 t0;
    xb_vecNx16 x0, x1, d0, d1, d2, d3, d4, h0;
    xb_vecNx16 p0, p1;
    unsigned int h01, h23, h45, h67;
    int n, m;
    valign h_align;
    vsaN   shft;
    vselN  sel0, sel1;

    NASSERT(N>0 && N % 16 == 0);
    NASSERT(M>0 && M % 16 == 0);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delay);
    NASSERT_ALIGN32(h);
    NASSERT_ALIGN32(srfir->p0);
    NASSERT_ALIGN32(srfir->p1);

    shft = BBE_MOVVSA32(14);
    pX = (const xb_vecNx16*)x;
    Y = (xb_vecNx16*)y;
    W0 = (xb_vecNx16*)srfir->p0;
    W1 = (xb_vecNx16*)srfir->p1;
    // Load the delay line state.
    // move delay line and put new samples
    WUR_CBEGIN((uintptr_t)delay);
    WUR_CEND((uintptr_t)(delay + 2 * 2 * (M + BBE_SIMD_WIDTH)));

    t0 = BBE_LVNX16_I(SEL, 0);
    sel0 = BBE_MOVVSELNX16(t0, 0);
    t0 = BBE_LVNX16_I(SEL, 2 * BBE_SIMD_WIDTH);
    sel1 = BBE_MOVVSELNX16(t0, 0);

    {
        int ix, sh;
        int16_t * p0;
        int16_t * p1;
        xb_vecNx16 t0;

        sh = 2 * BBE_SIMD_WIDTH * 3 / 4 - 1; //sh=23

        t0 = BBE_LPNX16_I(pX, 0);

        p0 = (int16_t*)MIN((uint32_t)W0, (uint32_t)W1);
        p1 = (int16_t*)MAX((uint32_t)W0, (uint32_t)W1);

        p0[2 * M + sh] = x[0];

        ix = p1 - delay;
        if (ix >= (2 * M + 2 * BBE_SIMD_WIDTH * 2)) p1 -= (2 * M + 2 * BBE_SIMD_WIDTH * 2);

        p1[sh] = x[0];
    }

    H = (const xb_vecNx16*)h;

    __Pragma("ymemory(pX)")
    __Pragma("ymemory(H)")
    for (n = 0; n<N / BBE_SIMD_WIDTH; n++)
    {
        h_align = BBE_LAVNX16_PP(H);

        BBE_LVNX16_IP(p0, pX, 2 * BBE_SIMD_WIDTH);
        p1 = BBE_LVNX16_I(pX, 0);

        d1 = BBE_SELNX16(p1, p0, sel1);
        d0 = BBE_SELNX16(p1, p0, sel0);

        BBE_SVNX16_IP(d0, W0, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(d0, W1, 2 * BBE_SIMD_WIDTH);

        BBE_SVNX16_IC(d1, W0);
        BBE_SVNX16_IC(d1, W1);

        R1 = (xb_vecNx16 *)XT_MAXU((uintptr_t)W0, (uintptr_t)W1);
        R0 = (xb_vecNx16 *)XT_MINU((uintptr_t)W0, (uintptr_t)W1);

        BBE_LVNX16_IP(d0, R0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(d1, R0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(d2, R0, 2 * BBE_SIMD_WIDTH);

        BBE_LVNX16_XP(d3, R1, -6 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(d4, R1, -2 * BBE_SIMD_WIDTH);

        d0 = BBE_SELNX16I(d1, d0, BBE_SELI_EXTRACT_LO_HALVES);

        BBE_MOVSAV(d0);
        BBE_MOVSDV(p0);
        BBE_ADDPNX16RRUMBCIAD(x0, x1, d2, d4);
        A = BBE_MOVWA32(1L << 13);

        __Pragma("ymemory(R0)")
        __Pragma("ymemory(R1)")
        __Pragma("ymemory(H)")
        __Pragma("loop_count min=1")
        for (m = 0; m<M / (BBE_SIMD_WIDTH); m++)
        {
            // take IR
            BBE_LAVNX16_XP(h0, h_align, H, BBE_SIMD_WIDTH);
            h01 = BBE_EXTRNX16C(h0, 0);
            h23 = BBE_EXTRNX16C(h0, 1);
            h45 = BBE_EXTRNX16C(h0, 2);
            h67 = BBE_EXTRNX16C(h0, 3);

            BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h01);
            BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h23);
            BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h45);

            BBE_LVNX16_IP(d0, R0, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(d1, R1, -2 * BBE_SIMD_WIDTH);

            BBE_ADDPNX16RRUMBC(x1, x0, d0, d1); BBE_MULANX16PR(A, x0, x1, h67);
        }

        H = (const xb_vecNx16*)((uintptr_t)H - M);

        x0 = BBE_PACKVNX40(A, shft);
        BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);
    }
    srfir->p0 = (int16_t*)W0;
    srfir->p1 = (int16_t*)W1;
}
#endif

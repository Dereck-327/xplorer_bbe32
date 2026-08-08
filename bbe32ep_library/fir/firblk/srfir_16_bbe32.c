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
DISCARD_FUN(void, srfir_process_16,  (void* _srfir, int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N))
#else

/*    M==16     */
void srfir_process_16 (void* _srfir, int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N)
{
    xb_vecNx16* restrict D;
    xb_vecNx16* restrict Y;
    const xb_vecNx16* restrict X;
    const xb_vecNx16* restrict X_;
    xb_vecNx40 A;
    xb_vecNx16 t0, x0, x1;
    uint32_t h01, h23, h45, h67;
    xb_vecNx16 a, b, c, d;
    int n;
    valign x_align;
    vsaN   shft;

    NASSERT(N>0 && !(N & 15));
    NASSERT(M>0 && M == 16);
    NASSERT_ALIGN16(y);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(delay);
    NASSERT_ALIGN16(h);

    X = (const xb_vecNx16*)x;
    X_ = (const xb_vecNx16*)(x + 1);
    // take IR
    t0 = BBE_LVNX16_I((const xb_vecNx16*)h, 0);
    h01 = BBE_EXTRNX16C(t0, 0);
    h23 = BBE_EXTRNX16C(t0, 1);
    h45 = BBE_EXTRNX16C(t0, 2);
    h67 = BBE_EXTRNX16C(t0, 3);
    D = (xb_vecNx16*)delay;
    Y = (xb_vecNx16*)y;
    // load delay line
    c = BBE_LVNX16_I(D, 0 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(d, X, 2 * BBE_SIMD_WIDTH);
    // begin filtering
    a = BBE_SELNX16I(d, c, BBE_SELI_ROTATE_RIGHT_1);
    x_align = BBE_LANX16_PP(X_);
    BBE_LANX16_IP(b, x_align, X_);
    BBE_MOVSBV(b);
    BBE_MOVSCV(c);
    shft = BBE_MOVVSA32(14);
    for (n = 0; n<N - BBE_SIMD_WIDTH; n += BBE_SIMD_WIDTH)
    {
        BBE_MOVSAV(a);
        BBE_MOVSDV(d);

        a = b;
        c = d;
        BBE_LVNX16_IP(d, X, 2 * BBE_SIMD_WIDTH);
        BBE_LANX16_IP(b, x_align, X_);

        BBE_ADDPNX16RRU(x1, x0); A = BBE_MULNX16PR(x0, x1, h01);
        BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h23);
        BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h45);
        BBE_ADDPNX16RRUMBC(x1, x0, b, c); BBE_MULANX16PR(A, x0, x1, h67);
        x0 = BBE_PACKVNX40(A, shft);
        BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);
    }
    BBE_MOVSAV(a);
    BBE_MOVSDV(d);
    BBE_ADDPNX16RRU(x1, x0); A = BBE_MULNX16PR(x0, x1, h01);
    BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h23);
    BBE_ADDPNX16RRU(x1, x0); BBE_MULANX16PR(A, x0, x1, h45);
    BBE_ADDPNX16RRUMBC(x1, x0, b, d); BBE_MULANX16PR(A, x0, x1, h67);
    x0 = BBE_PACKVNX40(A, shft);
    BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);

    // save last sample in the delay line
    BBE_SVNX16_I(d, D, 0);
}
#endif

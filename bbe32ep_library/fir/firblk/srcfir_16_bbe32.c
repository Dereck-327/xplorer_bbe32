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
    Block Complex FIR filter w/ Real Symmetric Impulse Response
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "srcfir_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR filter w/ Real Symmetric Impulse Response

Passes complex input data through a direct-form FIR filter with real 
coefficients and symmetric impulse response. The filter calculates N output 
samples using (M+1)/2 coefficients and requires last M+N-1 samples in the delay 
line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in complex samples) depends on FIR order M and is defined by srcfir_algDelay(M)

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
M           17,33 or a multiple of 16
N           Multiple of 8
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=16,17,32 and 33.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srcfir_init returns NULL handle.
---------------------------------------------------------------------------*/

#if !(HAVE_ABCD && HAVE_ADDPN && 1)
DISCARD_FUN(void, srcfir_process_16,(void* _srcfir,
                       int16_t * restrict delay, complex_fract16 * restrict y, 
                       const complex_fract16 * restrict x, const int16_t * restrict h,
                       int M, int N))
#else
//--------------------------------------------------------
// processing function for M==16
//--------------------------------------------------------
void srcfir_process_16(void* _srcfir,
                       int16_t * restrict delay, complex_fract16 * restrict y, 
                       const complex_fract16 * restrict x, const int16_t * restrict h,
                       int M, int N)
{
    xb_vecNx16* restrict D;
    xb_vecNx16* restrict Y;
    const xb_vecNx16* restrict X0;
    const xb_vecNx16* restrict X1;
    xb_vecNx40 A;
    xb_vecNx16 x0, x1, d0, d1, d2, t0, t1;
    vsaN vsa0;
    uint32_t h01, h23, h45, h67;
    int n;
    valign X1_va;

    NASSERT(N>0 && !(N & 7));
    NASSERT(M>0 && M == 16);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delay);
    NASSERT_ALIGN32(h);

    vsa0 = BBE_MOVVSA32(14);

    X0 = (const xb_vecNx16*)(x + 0);
    X1 = (const xb_vecNx16*)(x + 1);

    X1_va = BBE_LANX16_PP(X1);

    // take IR
    d0 = BBE_LVNX16_I((const xb_vecNx16*)h, 0);
    h01 = BBE_EXTRNX16C(d0, 0);
    h23 = BBE_EXTRNX16C(d0, 1);
    h45 = BBE_EXTRNX16C(d0, 2);
    h67 = BBE_EXTRNX16C(d0, 3);

    D = (xb_vecNx16*)delay;
    Y = (xb_vecNx16*)y;

    // load delay line
    d0 = BBE_LVNX16_I(D, 0 * 2 * BBE_SIMD_WIDTH);
    d1 = BBE_LVNX16_I(D, 1 * 2 * BBE_SIMD_WIDTH);

    // Pre-load the first data vector
    BBE_LVNX16_IP(d2, X0, 2 * BBE_SIMD_WIDTH);

    t0 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_2);
    t1 = BBE_SELNX16I(d2, d1, BBE_SELI_ROTATE_RIGHT_2);

    // begin filtering
    BBE_MOVSAV(t0);
    BBE_MOVSBV(t1);
    BBE_MOVSCV(d1);
    BBE_MOVSDV(d2);

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        // shift delay line and load the next data vector
        d0 = d1;
        d1 = d2;

        BBE_LVNX16_IP(d2, X0, 2 * BBE_SIMD_WIDTH);

        t0 = t1;

        BBE_LANX16_IP(t1, X1_va, X1);

        BBE_ADDPNX16RCU(x1, x0); A = BBE_MULNX16PR(x0, x1, h01);
        BBE_ADDPNX16RCU(x1, x0); BBE_MULANX16PR(A, x0, x1, h23);
        BBE_ADDPNX16RCU(x1, x0); BBE_MULANX16PR(A, x0, x1, h45);
        BBE_ADDPNX16RCUMBCIAD(x1, x0, t1, d1); BBE_MULANX16PR(A, x0, x1, h67);

        x0 = BBE_PACKVNX40(A, vsa0);
        BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);

        BBE_MOVSAV(t0);
        BBE_MOVSDV(d2);
    }

    // save last 2 samples in the delay line
    BBE_SVNX16_I(d0, D, 0);
    BBE_SVNX16_I(d1, D, 2 * BBE_SIMD_WIDTH);
}
#endif

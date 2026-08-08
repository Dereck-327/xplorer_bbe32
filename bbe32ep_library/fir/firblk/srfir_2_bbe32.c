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

/*    M==2     */
void srfir_process_2 (void* _srfir, int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N)
{
    xb_vecNx16* restrict D;
    xb_vecNx16* restrict Y;
    const xb_vecNx16* restrict X;
    xb_vecNx40 A;
    xb_vecNx16 x0;
    uint32_t h01;
    uint16_t h0;
    xb_vecNx16 a, c, d;
    int n;

    NASSERT(N>0 && !(N & 15));
    NASSERT(M>0 && M == 2);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delay);
    NASSERT_ALIGN32(h);

    X = (const xb_vecNx16*)x;
    // take IR
    h0 = h[0];
    h01 = (h0 << 16) | (h0);

    D = (xb_vecNx16*)delay;
    Y = (xb_vecNx16*)y;

    // load delay line
    c = BBE_LVNX16_I(D, 0);
    __Pragma("ymemory(X)")
    __Pragma("loop_count min=1")
    for (n = 0; n<N / (BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(d, X, 0*2 * BBE_SIMD_WIDTH);
        a = BBE_SELNX16I(d, c, BBE_SELI_ROTATE_LEFT_1);

        A = BBE_MULNX16PR(a, d, h01); 

        x0 = BBE_PACKQNX40(A);
        BBE_SVNX16_IP(x0, Y, 2 * BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(c, X, 2 * BBE_SIMD_WIDTH);
    }
    //update the delay line
    BBE_SVNX16_I(c, D, 0);
}

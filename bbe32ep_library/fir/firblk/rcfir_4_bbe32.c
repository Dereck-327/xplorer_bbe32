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
    Block Complex FIR Filter with Real Coefficients
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "rcfir_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR Filter with Real Coefficients

Computes a complex FIR filter (direct-form) using real IR stored in vector h.
The complex data input is stored in vector x. The filter output result is
stored in vector y. The filter calculates N output samples using M coefficients
and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in complex samples) depends on FIR order M and defined by rcfir[f]_algDelay(M)

Representation:
rcfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
rcfirf  IEEE-754 Std. single precision floating-point format for filter 
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[M]    Filter coefficients; h[0] is to be multiplied by the newest sample
N       Length of sample block
M       Length of filter
x[N]    Complex input samples
Output:
y[N]    Complex output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (rcfir) or 4 (rcfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, rcfir[f]_init returnd NULL handle.
-------------------------------------------------------------------------*/

/*    processing function for M==4 */
void rcfir_process_4(int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N)
{
    uint32_t h0, h1;
    int n;
    const xb_vecNx16* restrict X;
          xb_vecNx16* restrict Y;
    xb_vecNx16 d0,d1;
    xb_vecNx16 p0,p1,p2,p3;
    xb_vecNx16 y0;
    xb_vecNx16 h_;
    xb_vecNx40 w0;

    NASSERT(N>0 && !(N&7));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(delay);
    NASSERT_ALIGN32(h);

    X=(const xb_vecNx16*)x;
    Y=(      xb_vecNx16*)y;

    h_ = BBE_LVNX16_I((const xb_vecNx16*)h, 0);
    h0 = BBE_EXTRNX16C(h_, 0);
    h1 = BBE_EXTRNX16C(h_, 1);

    d0 = BBE_LVNX16_I((xb_vecNx16*)delay, 0);

    __Pragma("ymemory(X)")
    for (n=0; n<N/(BBE_SIMD_WIDTH/2); n++)
    {
        BBE_LVNX16_IP(d1, X, 2*BBE_SIMD_WIDTH);

        BBE_SELPCNX16I(p1, p0, d1, d0, 7);
        BBE_SELPCNX16I(p3, p2, d1, d0, 5);
        
        w0 = BBE_MULNX16PR(  p0, p1, h0);
        BBE_MULANX16PR(w0, p2, p3, h1);

        y0 = BBE_PACKQNX40(w0);
        BBE_SVNX16_IP(y0, Y, 2*BBE_SIMD_WIDTH);

        d0 = d1;
    }

    BBE_SVNX16_I(d0, (xb_vecNx16*)delay, 0);

}

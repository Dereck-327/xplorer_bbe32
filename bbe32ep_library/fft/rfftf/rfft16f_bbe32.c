/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
NatureDSP_Baseband library. FFT
Radix-2 forward FFT on complex_float data
C code optimized for BBE32
IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Tables of a twiddles */
#include "fft_fp_tw.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"


#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU
/* first stage rfft16f, ordinary DFT4 stage with real input  */
inline_ int first_stage_rfft16f_fp(const complex_float *tw, complex_float *x, /*complex_float *y, */ xb_vecN_4xcf32 *y2)
{
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 t0, t1, t2, t3, tw1, tw2, tw3;
    xb_vecN_4xcf32 d, s, a0, a1, a2, a3;

    BBE_LVN_4XCF32_XP(t0, p_src, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_XP(t1, p_src, 2 * BBE_SIMD_WIDTH);

    BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

    /* DFT4 real input - complex output */
    BBE_ADDSUBN_4XCF32(d, s, t0, t1);
    /* Make complex values with im = 0 */
    BBE_DSELN_4XCF32I(a1, a0, 0.0f, s, BBE_DSELI_INTERLEAVE_2);
    /* a3 = d(0:3) + 1j*d(4:end) */
    a3 = BBE_SHFLN_4XCF32I(d, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);
    BBE_ADDSUBN_4XCF32(a2, a0, a0, a1);
    /* a1 = conj(a3), since DFT4 input is real */
    a1 = BBE_MULJN_4XCF32(tw1, a3);
    a2 = BBE_MULN_4XCF32(a2, tw2);
    a3 = BBE_MULN_4XCF32(a3, tw3);

    BBE_DSELN_4XCF32I(t1, t0, a2, a0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(t3, t2, a3, a1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(a1, a0, t2, t0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(a3, a2, t3, t1, BBE_DSELI_INTERLEAVE_4);

    y2[0] = a0;
    y2[1] = a1;
    y2[2] = a2;
    y2[3] = a3;

    return 0;
}


/* last stage rfft16f, ordinary DFT4 stage, provide to output N/2+1 complex_float samples   */
inline_ int last_stage_rfft16f_fp(xb_vecN_4xcf32 *x2, complex_float *y)
{

    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    valign u = BBE_ZALIGN();

    xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
    xb_vecN_4xcf32 d0, d1, s0, s1;

    _x0 = x2[0];
    _x1 = x2[1];
    _x2 = x2[2];
    _x3 = x2[3];
    // DFT4
    BBE_ADDSUBN_4XCF32(d0, s0, _x0, _x2);
    BBE_ADDSUBN_4XCF32(d1, s1, _x1, _x3);
    d1 = BBE_SHFLN_4XCF32I(d1, BBE_SHFLI_SWAP_2);
    d1 = BBE_CONJN_4XCF32(d1);
    BBE_ADDSUBN_4XCF32(_x2, _x0, s0, s1);
    BBE_ADDSUBN_4XCF32(_x3, _x1, d0, d1);
    // Store N/2+1 complex samples (N=16)
    BBE_SAVN_4XCF32_XP(_x0, u, p_dst, 2 * BBE_SIMD_WIDTH);
    BBE_SAVN_4XCF32_XP(_x1, u, p_dst, 2 * BBE_SIMD_WIDTH);
    BBE_SAVN_4XCF32_XP(_x2, u, p_dst, sizeof(complex_float));
    BBE_SAN_4XCF32POS_FP(u, p_dst);

    return 0;
}

#endif //#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU

/*-------------------------------------------------------------------------
Radix-2 forward floating point FFT on real data

Description: These functions make FFT on real data of length N=2^n, n=4..15.
The algorithm exploits the symmetry properties of the FFT: first, a complex
FFT of half the original size is applied to input data, then the resulting
spectrum undergoes a postprocessing procedure which results in complex spectrum
of real input data.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:      
    x[N]          Real input signal
  Output:      
    y[(N/2+1)]    Output spectrum samples. 
  Returned value:
                  zero
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(int, rfft16f, (complex_float * restrict y, float32_t * restrict x) )
#else
int rfft16f(complex_float * restrict y, float32_t * restrict x)
{
    xb_vecN_4xcf32 tmp[4];
    complex_float *twiddle_table = (complex_float*)cfftf16_twd1;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT(x != (float32_t*)y);

    first_stage_rfft16f_fp(twiddle_table, (complex_float*)x, tmp);
    last_stage_rfft16f_fp(tmp, y);

    return 0;
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU

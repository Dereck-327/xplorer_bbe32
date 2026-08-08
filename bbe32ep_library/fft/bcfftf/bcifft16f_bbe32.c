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
/* Portable data types. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Common utility declarations. */
#include "common.h"
/* Twiddles tables for float point FFT */
#include "fft_fp_tw.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"

#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU
/* Last stage for bcifft16f */
inline_ int blk_last_stage_iDFT4_FP(complex_float *x, complex_float *y, int N, int L)
{
    int i;
    int l; 
    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    for (l = 0; l < L; l++)
    {

        for (i = 0; i<count; i++)
        {
            xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
            BBE_LVN_4XCF32_XP(_x0, p_src, stride);
            BBE_LVN_4XCF32_XP(_x1, p_src, stride);
            BBE_LVN_4XCF32_XP(_x2, p_src, stride);
            BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + 2 * BBE_SIMD_WIDTH);

            IDFT4_FP(_x0, _x1, _x2, _x3); 


            _SCALE_N_4XCF32(_x0, 1.0f / N);
            _SCALE_N_4XCF32(_x1, 1.0f / N);
            _SCALE_N_4XCF32(_x2, 1.0f / N);
            _SCALE_N_4XCF32(_x3, 1.0f / N);

            BBE_SVN_4XCF32_XP(_x0, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x1, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x2, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x3, p_dst, -3 * stride + 2 * BBE_SIMD_WIDTH);

        } //for (i = 0; i<count; i++)

        p_src = (xb_vecN_4xcf32*)(N*sizeof(*x) - 2 * BBE_SIMD_WIDTH*count + (uintptr_t)p_src);
        p_dst = (xb_vecN_4xcf32*)(N*sizeof(*x) - 2 * BBE_SIMD_WIDTH*count + (uintptr_t)p_dst);
    } // for (l = 0; l < L; l++)   
    __Pragma("no_reorder");
    return 0;
}



/* First stage radix4 for bcifft16f */
inline_ int blk_first_stage_iDFT4_N16_FP(
    const complex_float *tw,
    complex_float *x,
    complex_float *y,
    int N,
    int L)
{
    int  l;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = 3 + (xb_vecN_4xcf32 *)(x);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);

    xb_vecN_4xcf32  tw1, tw2, tw3;

    NASSERT(count == 1);
    NASSERT(N == 16);

    BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

    valign uu0, uu1;
    xb_vecNx16 t0, t1, t2, t3;

    for (l = 0; l < L; l++)
    {

        /*
        22 cycles per pipeline stage in steady state with unroll=2
        */
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
        BBE_LVN_4XCF32_XP(_x3, p_src, -2 * 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(_x1, p_src, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(_x2, p_src, -2 * 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(_x0, p_src, 7 * 2 * BBE_SIMD_WIDTH);

        IDFT4_FP(_x0, _x1, _x2, _x3);

        _x1 = BBE_MULJN_4XCF32(_x1, tw1);
        _x2 = BBE_MULJN_4XCF32(_x2, tw2);
        _x3 = BBE_MULJN_4XCF32(_x3, tw3);


        t0 = BBE_MOVNX16_FROMN_4XCF32(_x0);
        t1 = BBE_MOVNX16_FROMN_4XCF32(_x1);
        t2 = BBE_MOVNX16_FROMN_4XCF32(_x2);
        t3 = BBE_MOVNX16_FROMN_4XCF32(_x3);

        BBE_DSELNX16I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_4);

        uu0 = BBE_MOVUVR(t0);
        uu1 = BBE_MOVUVR(t1);
        BBE_SVINTLARNX16_XP(t2, uu0, p_dst, 2 * 2 * BBE_SIMD_WIDTH, 1);
        BBE_SVINTLARNX16_XP(t3, uu1, p_dst, -2 * BBE_SIMD_WIDTH, 1);

        BBE_SALIGNVRNX16_XP(t3, uu0, p_dst, 2 * 2 * BBE_SIMD_WIDTH);
        BBE_SALIGNVRNX16_XP(t3, uu1, p_dst, 2 * BBE_SIMD_WIDTH);

    } //  for (l = 0; l < L; l++)
    __Pragma("no_reorder");
    return 0;
}
#endif //#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point inverse FFT on complex data

Description: These functions make inverse FFT on L blocks, each of N=2^n
complex samples, where n=4..7

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input spectrum
  Output:          
    y[L][N]   Output complex data
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(void, bcifft16f, (complex_float * restrict y, complex_float * restrict x, int L) )
#else
void bcifft16f(complex_float * restrict y, complex_float * restrict x, int L)
{
    const complex_float *tw_tab = (const complex_float*)cfftf16_twd1;

    NASSERT(x != y);
    NASSERT_ALIGN(x, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(y, BBE_SIMD_WIDTH * 2);

    blk_first_stage_iDFT4_N16_FP(tw_tab, x, y, 16, L);
    blk_last_stage_iDFT4_FP(y, y, 16, L); 
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU



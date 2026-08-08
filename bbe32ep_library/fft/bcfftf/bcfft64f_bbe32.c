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
/* Last stage radix-4 */
inline_ int blk_last_stage_DFT4_FP(complex_float *x, complex_float *y, int N, int L)
{
    int i;
    int l; 
    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);

    xb_vecN_4xcf32 * p_src2 = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst2 = (xb_vecN_4xcf32 *)(y);

    ASSERT( (count&1)==0 ); 
    for (i = 0; i<count/2; i++)
    {
        p_src = (xb_vecN_4xcf32*)(x + 4 * i);
        p_dst = (xb_vecN_4xcf32*)(y + 4 * i);
        p_src2 = (xb_vecN_4xcf32*)(x + 4 * (i + count / 2) );
        p_dst2 = (xb_vecN_4xcf32*)(y + 4 * (i + count / 2) );
        /* 9 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {
            xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
            /*
            p_src = (xb_vecN_4xcf32*)(x + N*l + 4 * i);
            p_dst = (xb_vecN_4xcf32*)(y + N*l + 4 * i);
            */
            BBE_LVN_4XCF32_XP(_x0, p_src, stride);
            BBE_LVN_4XCF32_XP(_x1, p_src, stride);
            BBE_LVN_4XCF32_XP(_x2, p_src, stride);
            BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + N*sizeof(*x));

            DFT4_MULM_FP(_x0, _x1, _x2, _x3);

            BBE_SVN_4XCF32_XP(_x0, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x1, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x2, p_dst, stride);
            BBE_SVN_4XCF32_XP(_x3, p_dst, -3 * stride + N*sizeof(*x));

            BBE_LVN_4XCF32_XP(_x0, p_src2, stride);
            BBE_LVN_4XCF32_XP(_x1, p_src2, stride);
            BBE_LVN_4XCF32_XP(_x2, p_src2, stride);
            BBE_LVN_4XCF32_XP(_x3, p_src2, -3 * stride + N*sizeof(*x));

            DFT4_MULM_FP(_x0, _x1, _x2, _x3);

            BBE_SVN_4XCF32_XP(_x0, p_dst2, stride);
            BBE_SVN_4XCF32_XP(_x1, p_dst2, stride);
            BBE_SVN_4XCF32_XP(_x2, p_dst2, stride);
            BBE_SVN_4XCF32_XP(_x3, p_dst2, -3 * stride + N*sizeof(*x));
        } // for (l = 0; l < L; l++)
    } // for (i = 0; i<count/2; i++)

    __Pragma("no_reorder");
    return 0;
} /* blk_last_stage_DFT4_FP */
#endif //#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point forward FFT on complex data

Description: These functions make forward FFT on L blocks, each of N=2^n
complex samples, where n=4..7. 

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input signal
  Output:          
    y[L][N]   Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(void, bcfft64f, (complex_float * restrict y, complex_float * restrict x, int L) )
#else
void bcfft64f(complex_float * restrict y, complex_float * restrict x, int L)
{
    const int N = 64;
    const complex_float *tw_tab  = (const complex_float*)cfftf64_twd1;
    const complex_float *tw_tab2 = (const complex_float*)cfftf64_twd2;

    NASSERT(x != y);
    NASSERT_ALIGN(x, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(y, BBE_SIMD_WIDTH * 2);

    blk_first_stage_DFT4_FP(tw_tab, x, y, N, L);
    blk_inner_stage_DFT4_FP(tw_tab2, y, x, N, L); 
    blk_last_stage_DFT4_FP(x, y, N, L);
}
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || ! HAVE_VFPU




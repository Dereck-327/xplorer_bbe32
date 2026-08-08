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
Internal components for the floating point FFT
*/
#ifndef __FFT_FP_COMMON_H__
#define __FFT_FP_COMMON_H__

#include "NatureDSP_types.h"
#include "common.h"

#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU

#ifdef IS_INV_FFTF

inline_ int stage_last_iDFT4xIN_4_FP(complex_float *x, /* Input */
    complex_float *y, /* Output, can coincide with x */
    int N             /* Transform size */
    ) ATTRIBUTE_ALWAYS_INLINE;

inline_ int first_stage_iDFT4_FP(const complex_float *tw, /* Twiddles table */
    complex_float *x,  /* Input */
    complex_float *y,  /* Output */
    int N)  ATTRIBUTE_ALWAYS_INLINE;

#else //#ifdef IS_INV_FFTF
inline_ int stage_last_DFT4xIN_4_FP(complex_float *x, /* Input */
    complex_float *y, /* Output, can coincide with x */
    int N             /* Transform size */
    ) ATTRIBUTE_ALWAYS_INLINE;

inline_ int first_stage_DFT4_FP (   const complex_float *tw, /* Twiddles table */
                                    complex_float *x,  /* Input */
                                    complex_float *y,  /* Output */
                                    int N)  ATTRIBUTE_ALWAYS_INLINE;

#endif //#ifdef IS_INV_FFTF

#ifdef IS_INV_FFTF
inline_ int stage_last_iDFT8xIN_8_FP
#else
inline_ int stage_last_DFT8xIN_8_FP
#endif
(complex_float *x,  /* Input */
complex_float *y,  /* Output, can coincide with x */
int N             /* Transform size */) ATTRIBUTE_ALWAYS_INLINE;

inline_ void rifft_spec_conv_fp(complex_float *x,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N) ATTRIBUTE_ALWAYS_INLINE;

#ifdef IS_INV_FFTF
inline_ int stage_inner_R3_iDFT4xIv_FP(
#else
inline_ int stage_inner_R3_DFT4xIv_FP(
#endif
    const complex_float *_tw, /* Twiddles table*/
    complex_float *x                                            /*input*/,
    complex_float *y                                            /*output*/,
    const int N,
    const int v                                                 /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    ) ATTRIBUTE_ALWAYS_INLINE;

#ifdef IS_INV_FFTF
inline_ int first_stage_iDFT4_FP_v2
#else
inline_ int first_stage_DFT4_FP_v2
#endif
(const complex_float *tw, complex_float *x, complex_float *y, int N) ATTRIBUTE_ALWAYS_INLINE;

#ifdef IS_INV_FFTF
inline_ int stage_inner_merged_iDFT4_FP
#else
inline_ int stage_inner_merged_DFT4_FP
#endif
(const complex_float *_tw,
complex_float *x   /*input*/,
complex_float *y   /*output*/,
const int N,
const int v /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
)  ATTRIBUTE_ALWAYS_INLINE;

#ifdef IS_INV_FFTF
inline_ int stage_inner_iDFT4_v16_FP
#else
inline_ int stage_inner_DFT4_v16_FP
#endif
(const complex_float *_tw,
complex_float *x   /*input*/,
complex_float *y   /*output*/,
const int N,
const int v /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
)  ATTRIBUTE_ALWAYS_INLINE;

inline_ void rfft_spec_conv_fp(complex_float *y,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N) ATTRIBUTE_ALWAYS_INLINE;



/*  Multiply by {-j, -j, -j, -j} */
inline_ void  mul_1j(xb_vecN_4xcf32 *x /* in/out */)
{
    xb_vecN_2xf32 z0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(*x));

   
    z0 = BBE_MULMN_2XF32(z0, BBE_CONSTN_2XF32(1), 2, 6);
    *x = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(z0));
}

/* IDFT4xI4 */
#define IDFT4_FP(_x0, _x1, _x2, _x3)\
{                                                         \
    xb_vecN_4xcf32 _d0, _d1, _s0, _s1;                    \
    BBE_ADDSUBN_4XCF32(_d0, _s0, _x0, _x2);               \
    BBE_ADDSUBN_4XCF32(_d1, _s1, _x1, _x3);               \
    _d1 = BBE_SHFLN_4XCF32I(_d1, BBE_SHFLI_SWAP_2);       \
    _d1 = BBE_CONJN_4XCF32(_d1);                          \
    BBE_ADDSUBN_4XCF32(_x2, _x0, _s0, _s1);               \
    BBE_ADDSUBN_4XCF32(_x1, _x3, _d0, _d1);               \
}

/* DFT4xI4 floating point with MULM */

#define DFT4_MULM_FP(_x0, _x1, _x2, _x3)\
{                                                         \
    xb_vecN_4xcf32 _d0, _d1, _s0, _s1;                    \
    BBE_ADDSUBN_4XCF32(_d0, _s0, _x0, _x2);               \
    BBE_ADDSUBN_4XCF32(_d1, _s1, _x1, _x3);               \
    mul_1j(&_d1);                                         \
                                                          \
    BBE_ADDSUBN_4XCF32(_x2, _x0, _s0, _s1);               \
    BBE_ADDSUBN_4XCF32(_x3, _x1, _d0, _d1);               \
}

/* DFT4xI4 floating point */
#define DFT4_FP(_x0, _x1, _x2, _x3)\
{                                                         \
    xb_vecN_4xcf32 _d0, _d1, _s0, _s1;                    \
    BBE_ADDSUBN_4XCF32(_d0, _s0, _x0, _x2);               \
    BBE_ADDSUBN_4XCF32(_d1, _s1, _x1, _x3);               \
    _d1 = BBE_SHFLN_4XCF32I(_d1, BBE_SHFLI_SWAP_2);       \
    _d1 = BBE_CONJN_4XCF32(_d1);                          \
    BBE_ADDSUBN_4XCF32(_x2, _x0, _s0, _s1);               \
    BBE_ADDSUBN_4XCF32(_x3, _x1, _d0, _d1);               \
}


/*
Calculate x*=factor:
x - complex float vector (xb_vecN_2xf32);
factor - float, single precision.
*/
#define _SCALE_N_4XCF32(__x/*in/out*/, factor)                                    \
{                                                                               \
    xb_vecN_2xf32 tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(__x)); \
    tmp = BBE_MULN_2XF32(tmp, (factor));                                        \
    __x = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));               \
}

/*
    Last stage DFT8-based. 
*/
#ifdef IS_INV_FFTF
inline_ int stage_last_iDFT8xIN_8_FP
#else
inline_ int stage_last_DFT8xIN_8_FP
#endif
                                    (complex_float *x,  /* Input */
                                    complex_float *y,  /* Output, can coincide with x */
                                     int N             /* Transform size */ )
{
    int i;
    const int stride = N / 8 * sizeof(*x);
    int count = N / (8 * BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x));

    xb_vecN_4xcf32 * p_src0 = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_src1 = (xb_vecN_4xcf32 *)(x + stride / sizeof(*x));
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    xb_vecN_4xcf32 a0, a1, b0, b1, c0, c1, d0, d1;

    /* Algorithm:
    function Y = dft8(x)
    %
    A = fft(x(1:2:end));


    % B = fft(x(2:2:end));

    B = zeros(4,1);
    x0 = x(2);
    x1 = x(4);
    x2 = x(6);
    x3 = x(8);


    s0 = x0+x2;
    s1 = x1+x3

    d0 = x0-x2;
    d1 = x1-x3

    B(1) = s0+s1;
    B(3) = s0-s1;
    B(2) = (-1j)*(d0 + d1) + (d0 - d1); 
    B(4) = (-1j)*(d0 + d1) - (d0 - d1); 


    B = B .* [ 1;
    1/sqrt(2);
    -1j;
    1/sqrt(2);]

    Y = zeros(8,1)
    Y(1:4) = A + B;
    Y(5:8) = A - B;
    */
    /* cifftf 20 cycles per pipeline stage in steady state with unroll=1 */
    /* cfftf: 16 cycles per pipeline stage in steady state with unroll=1 */
    for (i = 0; i<count; i++)
    {
        xb_vecN_4xcf32 sd, dd;

        BBE_LVN_4XCF32_XP(a0, p_src0, 2 * stride);
        BBE_LVN_4XCF32_XP(b0, p_src0, 2 * stride);
        BBE_LVN_4XCF32_XP(c0, p_src0, 2 * stride);
        BBE_LVN_4XCF32_XP(d0, p_src0, -6 * stride + 2 * BBE_SIMD_WIDTH);

        BBE_LVN_4XCF32_XP(a1, p_src1, 2 * stride);
        BBE_LVN_4XCF32_XP(b1, p_src1, 2 * stride);
        BBE_LVN_4XCF32_XP(c1, p_src1, 2 * stride);
        BBE_LVN_4XCF32_XP(d1, p_src1, -6 * stride + 2 * BBE_SIMD_WIDTH);

        DFT4_FP(a0, b0, c0, d0);
        {
            xb_vecN_4xcf32 _d0, _d1, _s0, _s1;
            BBE_ADDSUBN_4XCF32(_d0, _s0, a1, c1);
            BBE_ADDSUBN_4XCF32(_d1, _s1, b1, d1);

            BBE_ADDSUBN_4XCF32(dd, sd, _d0, _d1);
            // sd *= -j
            sd = BBE_SHFLN_4XCF32I(sd, BBE_SHFLI_SWAP_2);
            sd = BBE_CONJN_4XCF32(sd);

            BBE_ADDSUBN_4XCF32(d1, b1, sd, dd);
#ifdef IS_INV_FFTF
            _SCALE_N_4XCF32(b1, 0.707106781186547f/N);
            _SCALE_N_4XCF32(d1, 0.707106781186547f/N);
#else
            _SCALE_N_4XCF32(b1, 0.707106781186547f);
            _SCALE_N_4XCF32(d1, 0.707106781186547f);
#endif

            BBE_ADDSUBN_4XCF32(c1, a1, _s0, _s1);
        }

        // c1 *= 1j; 
        c1 = BBE_SHFLN_4XCF32I(c1, BBE_SHFLI_SWAP_2);
        c1 = BBE_CONJN_4XCF32(c1);

#ifdef IS_INV_FFTF
        _SCALE_N_4XCF32(a0, 1.0f / N);
        _SCALE_N_4XCF32(d0, 1.0f / N);
        _SCALE_N_4XCF32(c0, 1.0f / N);
        //_SCALE_N_4XCF32(b1, 1.0f / N);
        _SCALE_N_4XCF32(a1, 1.0f / N);
        //_SCALE_N_4XCF32(d1, 1.0f / N);
        _SCALE_N_4XCF32(c1, 1.0f / N);
        _SCALE_N_4XCF32(b0, 1.0f / N);
#endif

        BBE_ADDSUBN_4XCF32(a1, a0, a0, a1); 
        BBE_ADDSUBN_4XCF32(b1, b0, b0, b1); 
        BBE_ADDSUBN_4XCF32(c1, c0, c0, c1); 
        BBE_ADDSUBN_4XCF32(d1, d0, d0, d1); 



#ifdef IS_INV_FFTF
        BBE_SVN_4XCF32_XP(a0, p_dst, stride);
        BBE_SVN_4XCF32_XP(d1, p_dst, stride);
        BBE_SVN_4XCF32_XP(c1, p_dst, stride);
        BBE_SVN_4XCF32_XP(b1, p_dst, stride);

        BBE_SVN_4XCF32_XP(a1, p_dst, stride);
        BBE_SVN_4XCF32_XP(d0, p_dst, stride);
        BBE_SVN_4XCF32_XP(c0, p_dst, stride);
        BBE_SVN_4XCF32_XP(b0, p_dst, -7 * stride + 2 * BBE_SIMD_WIDTH);
#else
        BBE_SVN_4XCF32_XP(a0, p_dst, stride);
        BBE_SVN_4XCF32_XP(b0, p_dst, stride);
        BBE_SVN_4XCF32_XP(c0, p_dst, stride);
        BBE_SVN_4XCF32_XP(d0, p_dst, stride);

        BBE_SVN_4XCF32_XP(a1, p_dst, stride);
        BBE_SVN_4XCF32_XP(b1, p_dst, stride);
        BBE_SVN_4XCF32_XP(c1, p_dst, stride);
        BBE_SVN_4XCF32_XP(d1, p_dst, -7 * stride + 2 * BBE_SIMD_WIDTH);
#endif

    }
    return 0;
} //stage_last_DFT8xIN_8_FP

/*
    Last stage DFT4-based.
*/
inline_ int stage_last_DFT4xIN_4_FP( complex_float *x, /* Input */
                                     complex_float *y, /* Output, can coincide with x */
                                     int N             /* Transform size */ 
                                     )
{
    int i;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    /*  17 cycles per pipeline stage in steady state with unroll=4 */
    for (i = 0; i<count; i++)
    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
        BBE_LVN_4XCF32_XP(_x0, p_src, stride);
        BBE_LVN_4XCF32_XP(_x1, p_src, stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + 2 * BBE_SIMD_WIDTH);

        {
            xb_vecN_4xcf32 d0, d1, s0, s1;

            BBE_ADDSUBN_4XCF32(d0, s0, _x0, _x2);
            BBE_ADDSUBN_4XCF32(d1, s1, _x1, _x3);
            mul_1j(&d1); 
            //d1 = BBE_SHFLN_4XCF32I(d1, BBE_SHFLI_SWAP_2);
            //d1 = BBE_CONJN_4XCF32(d1);
            BBE_ADDSUBN_4XCF32(_x2, _x0, s0, s1);
            BBE_ADDSUBN_4XCF32(_x3, _x1, d0, d1);
        }

        BBE_SVN_4XCF32_XP(_x0, p_dst, stride);
        BBE_SVN_4XCF32_XP(_x1, p_dst, stride);
        BBE_SVN_4XCF32_XP(_x2, p_dst, stride);
        BBE_SVN_4XCF32_XP(_x3, p_dst, -3 * stride + 2 * BBE_SIMD_WIDTH);

    } //for (i = 0; i<count; i++)

    return 0;
} //stage_last_DFT4xIN_4_FP

inline_ int stage_last_iDFT4xIN_4_FP(complex_float *x, complex_float *y, int N)
{
    int i;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    /* 12 cycles per pipeline stage in steady state with unroll = 2 */
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

    return 0;
} // stage_last_iDFT4xIN_4_FP

/*
    Inner stage , DFT4 based,
    x,y  - Must not overlap and must be aligned on 32-byte boundary
*/
#ifdef IS_INV_FFTF
inline_ int stage_inner_R3_iDFT4xIv_FP(
#else
inline_ int stage_inner_R3_DFT4xIv_FP(
#endif
    const complex_float *_tw, /* Twiddles table*/
    complex_float *x                                            /*input*/,
    complex_float *y                                            /*output*/,
    const int N,
    const int v                                                 /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    )
{
    int i, j;
    xb_vecN_4xcf32 *tw = (xb_vecN_4xcf32 *)_tw;

    xb_vecN_4xcf32 *px = (xb_vecN_4xcf32 *)x;
    xb_vecN_4xcf32 *py = (xb_vecN_4xcf32 *)y;

    valign u = BBE_LAVNX16_PP((xb_vecNx16*)tw);

    const int stride_bytes = N / 4 * sizeof(*x);
    const int num_bfls = N / (4 * v);
    int M = v / (BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x)); // Number of iterations in the inner loop

    if (M == 1)
    {
        /* Second stage */
        px = (xb_vecN_4xcf32*)(stride_bytes + (uintptr_t)px);
        for (i = 0; i < num_bfls; i++)
        {
            xb_vecN_4xcf32  tw1, tw2, tw3;
            xb_vecN_4xcf32 tmp2;
            xb_vecN_4xcf32 x0, x1, x2, x3;

            BBE_LAVN_4XCF32_XP(tmp2, u, tw, 3 * sizeof(*x));
            /* 9 cycles per pipeline stage in steady state */
            tw1 = BBE_REPN_4XCF32(tmp2, 0);
            tw2 = BBE_REPN_4XCF32(tmp2, 1);
            tw3 = BBE_REPN_4XCF32(tmp2, 2);

            /* 8 cycles per pipeline stage in steady state with unroll = 1 */
            BBE_LVN_4XCF32_XP(x1, px, 2 * stride_bytes);
            BBE_LVN_4XCF32_XP(x3, px, -3 * stride_bytes);
            BBE_LVN_4XCF32_XP(x0, px, 2 * stride_bytes);
            BBE_LVN_4XCF32_XP(x2, px, -1 * stride_bytes + 2 * BBE_SIMD_WIDTH);
#ifdef IS_INV_FFTF
            IDFT4_FP(x0, x1, x2, x3);
            x1 = BBE_MULJN_4XCF32(x1, tw1);
            x2 = BBE_MULJN_4XCF32(x2, tw2);
            x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
            DFT4_FP(x0, x1, x2, x3);
            x1 = BBE_MULN_4XCF32(x1, tw1);
            x2 = BBE_MULN_4XCF32(x2, tw2);
            x3 = BBE_MULN_4XCF32(x3, tw3);
#endif



            BBE_SVN_4XCF32_IP(x0, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x1, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x2, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x3, py, 2 * BBE_SIMD_WIDTH);
        } //for (i = 0; i < num_bfls; i++)
    }
    else if (M == 4)
    {   
        /* Third stage */
        for (j = 0; j < M; j++)
        {
            tw = (xb_vecN_4xcf32*)_tw;
            u = BBE_LAVNX16_PP((xb_vecNx16*)tw);
            py = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*j + (uintptr_t)y);
            px = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*j + stride_bytes + (uintptr_t)x);
            for (i = 0; i < num_bfls; i++)
            {
                xb_vecN_4xcf32  tw1, tw2, tw3;
                xb_vecN_4xcf32 x0, x1, x2, x3;
                xb_vecN_4xcf32 tmp2;

                /* 9 cycles per pipeline stage in steady state with unroll=1 */
                BBE_LAVN_4XCF32_XP(tmp2, u, tw, 3 * sizeof(*x));

                tw1 = BBE_REPN_4XCF32(tmp2, 0);
                tw2 = BBE_REPN_4XCF32(tmp2, 1);
                tw3 = BBE_REPN_4XCF32(tmp2, 2);

            /*  px = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*(j + i*M) + (uintptr_t)x);
                py = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*(j + i*4*M) + (uintptr_t)y); */

                BBE_LVN_4XCF32_XP(x1, px, 2 * stride_bytes);
                BBE_LVN_4XCF32_XP(x3, px, -3 * stride_bytes);
                BBE_LVN_4XCF32_XP(x0, px, 2 * stride_bytes);
                BBE_LVN_4XCF32_XP(x2, px, -1 * stride_bytes + M * 2 * BBE_SIMD_WIDTH);

#ifdef IS_INV_FFTF
                IDFT4_FP(x0, x1, x2, x3);
                x1 = BBE_MULJN_4XCF32(x1, tw1);
                x2 = BBE_MULJN_4XCF32(x2, tw2);
                x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
                DFT4_FP(x0, x1, x2, x3);
                x1 = BBE_MULN_4XCF32(x1, tw1);
                x2 = BBE_MULN_4XCF32(x2, tw2);
                x3 = BBE_MULN_4XCF32(x3, tw3);
#endif

                BBE_SVN_4XCF32_IP(x0, py, 16 * sizeof(*x));
                BBE_SVN_4XCF32_IP(x1, py, 16 * sizeof(*x));
                BBE_SVN_4XCF32_IP(x2, py, 16 * sizeof(*x));
                BBE_SVN_4XCF32_IP(x3, py, 16 * sizeof(*x));
            } // for (i = 0; i < num_bfls; i++)
        } //for (j = 0; j<M; j++)
    } // else if (M == 4)
    else
    {
        px = (xb_vecN_4xcf32*)(stride_bytes + (uintptr_t)px);
        for (i = 0; i < num_bfls; i++)
        {
            xb_vecN_4xcf32  tw1, tw2, tw3;
            xb_vecN_4xcf32 tmp2;

            BBE_LAVN_4XCF32_XP(tmp2, u, tw, 3 * sizeof(*x));

            tw1 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_0X4);
            tw2 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_1X4);
            tw3 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_2X4);

            /* 17 cycles per pipeline stage in steady state with unroll=2 */
            for (j = 0; j<M; j++)
            {
                xb_vecN_4xcf32 x0, x1, x2, x3;
                BBE_LVN_4XCF32_XP(x1, px, 2 * stride_bytes);
                BBE_LVN_4XCF32_XP(x3, px, -3 * stride_bytes);
                BBE_LVN_4XCF32_XP(x0, px, 2 * stride_bytes);
                BBE_LVN_4XCF32_XP(x2, px, -1 * stride_bytes + 2 * BBE_SIMD_WIDTH);

#ifdef IS_INV_FFTF
                IDFT4_FP(x0, x1, x2, x3);
                x1 = BBE_MULJN_4XCF32(x1, tw1);
                x2 = BBE_MULJN_4XCF32(x2, tw2);
                x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
                DFT4_FP(x0, x1, x2, x3);
                x1 = BBE_MULN_4XCF32(x1, tw1);
                x2 = BBE_MULN_4XCF32(x2, tw2);
                x3 = BBE_MULN_4XCF32(x3, tw3);
#endif


                BBE_SVN_4XCF32_XP(x0, py, v * sizeof(*x));
                BBE_SVN_4XCF32_XP(x1, py, v * sizeof(*x));
                BBE_SVN_4XCF32_XP(x2, py, v * sizeof(*x));
                BBE_SVN_4XCF32_XP(x3, py, 2 * BBE_SIMD_WIDTH - 3 * v * sizeof(*x));

            } //for (j = 0; j<M; j++)
            py += (4 - 1) * v * sizeof(*x) / sizeof(*py);
        } //for (i = 0; i < num_bfls; i++)
    }
    __Pragma("no_reorder");
    return 0;
} /* stage_inner_R3_DFT4xIv_FP */

/* First stage DFT4-based, best for N<128 */
#ifdef IS_INV_FFTF
inline_ int first_stage_iDFT4_FP
#else
inline_ int first_stage_DFT4_FP
#endif
                                   (const complex_float *tw, /* Twiddles table */
                                      complex_float *x,  /* Input */
                                      complex_float *y,  /* Output */
                                      int N)
{
    int i;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);

    xb_vecN_4xcf32 t0, t1, t2, t3, tw1, tw2, tw3;

    /* 13 cycles per pipeline stage in steady state with unroll=1 */
    for (i = 0; i<count; i++)
    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
        BBE_LVN_4XCF32_XP(_x0, p_src, stride);
        BBE_LVN_4XCF32_XP(_x1, p_src, stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + 2 * BBE_SIMD_WIDTH);

        BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);
#ifdef IS_INV_FFTF
        IDFT4_FP(_x0, _x1, _x2, _x3); 
        _x1 = BBE_MULJN_4XCF32(_x1, tw1);
        _x2 = BBE_MULJN_4XCF32(_x2, tw2);
        _x3 = BBE_MULJN_4XCF32(_x3, tw3);
#else
        DFT4_FP(_x0, _x1, _x2, _x3);
        _x1 = BBE_MULN_4XCF32(_x1, tw1);
        _x2 = BBE_MULN_4XCF32(_x2, tw2);
        _x3 = BBE_MULN_4XCF32(_x3, tw3);
#endif



        BBE_DSELN_4XCF32I(t1, t0, _x2, _x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(t3, t2, _x3, _x1, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x1, _x0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x3, _x2, t3, t1, BBE_DSELI_INTERLEAVE_4);

        BBE_SVN_4XCF32_IP(_x0, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x1, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x2, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x3, p_dst, 2 * BBE_SIMD_WIDTH);

    }
    __Pragma("no_reorder");
    return 0;
}  /* first_stage_DFT4_FP */


/* First stage radix-4, with SVINT, best for N>=128 */
#ifdef IS_INV_FFTF
inline_ int first_stage_iDFT4_FP_v2
#else
inline_ int first_stage_DFT4_FP_v2
#endif
(const complex_float *tw, complex_float *x, complex_float *y, int N)
{
    int i;
    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(stride * 1 + (uintptr_t)x);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);
    valign uu0, uu1;
    xb_vecNx16 t0, t1, t2, t3;
    xb_vecN_4xcf32  tw1, tw2, tw3;


    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;

        BBE_LVN_4XCF32_XP(_x1, p_src, 2 * stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, -2 * stride);
        BBE_LVN_4XCF32_XP(_x0, p_src, stride + 2 * BBE_SIMD_WIDTH);

        BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

#ifdef IS_INV_FFTF
        IDFT4_FP(_x0, _x1, _x2, _x3);
        _x1 = BBE_MULJN_4XCF32(_x1, tw1);
        _x2 = BBE_MULJN_4XCF32(_x2, tw2);
        _x3 = BBE_MULJN_4XCF32(_x3, tw3);
#else
        DFT4_FP(_x0, _x1, _x2, _x3);
        _x1 = BBE_MULN_4XCF32(_x1, tw1);
        _x2 = BBE_MULN_4XCF32(_x2, tw2);
        _x3 = BBE_MULN_4XCF32(_x3, tw3);
#endif

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

    }
    /* 11 cycles per pipeline stage in steady state with unroll=1 */
    for (i = 0; i<count - 1; i++)
    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;

        BBE_LVN_4XCF32_XP(_x1, p_src, 2 * stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, -2 * stride);
        BBE_LVN_4XCF32_XP(_x0, p_src, stride + 2 * BBE_SIMD_WIDTH);

        BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

#ifdef IS_INV_FFTF
        IDFT4_FP(_x0, _x1, _x2, _x3);
        _x1 = BBE_MULJN_4XCF32(_x1, tw1);
        _x2 = BBE_MULJN_4XCF32(_x2, tw2);
        _x3 = BBE_MULJN_4XCF32(_x3, tw3);
#else
        DFT4_FP(_x0, _x1, _x2, _x3);
        _x1 = BBE_MULN_4XCF32(_x1, tw1);
        _x2 = BBE_MULN_4XCF32(_x2, tw2);
        _x3 = BBE_MULN_4XCF32(_x3, tw3);
#endif

        t0 = BBE_MOVNX16_FROMN_4XCF32(_x0);
        t1 = BBE_MOVNX16_FROMN_4XCF32(_x1);
        t2 = BBE_MOVNX16_FROMN_4XCF32(_x2);
        t3 = BBE_MOVNX16_FROMN_4XCF32(_x3);

        BBE_DSELNX16I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_4);

        BBE_SALIGNVRNX16_XP(t0, uu0, p_dst, 2 * 2 * BBE_SIMD_WIDTH);
        BBE_SALIGNVRNX16_XP(t1, uu1, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVINTLARNX16_XP(t2, uu0, p_dst, 2 * 2 * BBE_SIMD_WIDTH, 1);
        BBE_SVINTLARNX16_XP(t3, uu1, p_dst, -2 * BBE_SIMD_WIDTH, 1);
    }
    BBE_SALIGNVRNX16_XP(t3, uu0, p_dst, 2 * 2 * BBE_SIMD_WIDTH);
    BBE_SALIGNVRNX16_XP(t3, uu1, p_dst, 2 * BBE_SIMD_WIDTH);

    __Pragma("no_reorder");
    return 0;
} /* first_stage_DFT4_FP_v2 */


/*
    Inner stage radix-4  with merged loop. 
    format of the twiddles differ from other stage_inner*, since 
    used BBE_LVN_4XCF32_IC instead BBE_LAVN.
    Doesn't work for N > 8192
*/
#ifdef IS_INV_FFTF
inline_ int stage_inner_merged_iDFT4_FP
#else
inline_ int stage_inner_merged_DFT4_FP
#endif
    (const complex_float *_tw,
    complex_float *x   /*input*/,
    complex_float *y   /*output*/,
    const int N,
    const int v /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    )
{
    int i;
    const int vect_sz = (BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x));
    xb_vecN_4xcf32 *tw = (xb_vecN_4xcf32 *)_tw;

    const int stride_bytes = N / 4 * sizeof(*x);

    xb_vecN_4xcf32 *px0 = (xb_vecN_4xcf32 *)(0 * stride_bytes + (uintptr_t)x);
    xb_vecN_4xcf32 *px1 = (xb_vecN_4xcf32 *)(1 * stride_bytes + (uintptr_t)x);
    xb_vecN_4xcf32 *px2 = (xb_vecN_4xcf32 *)(2 * stride_bytes + (uintptr_t)x);
    xb_vecN_4xcf32 *px3 = (xb_vecN_4xcf32 *)(3 * stride_bytes + (uintptr_t)x);

    xb_vecN_4xcf32 *py0 = (xb_vecN_4xcf32 *)(0 * v * sizeof(*y) + (uintptr_t)y);
    xb_vecN_4xcf32 *py1 = (xb_vecN_4xcf32 *)(1 * v * sizeof(*y) + (uintptr_t)y);
    xb_vecN_4xcf32 *py2 = (xb_vecN_4xcf32 *)(2 * v * sizeof(*y) + (uintptr_t)y);
    xb_vecN_4xcf32 *py3 = (xb_vecN_4xcf32 *)(3 * v * sizeof(*y) + (uintptr_t)y);

    int M = v / vect_sz; // Number of iterations in the inner loop

    int offset_x = 0;
    int offset_y = 0;

    unsigned mod_add_x = (M * 2 * BBE_SIMD_WIDTH + ((N / 4 * sizeof(*x) - 2 * BBE_SIMD_WIDTH) << 16));
    unsigned mod_add_y = (4 * v * sizeof(*x) + ((N * sizeof(*x) - 2 * BBE_SIMD_WIDTH) << 16));

    xb_vecN_4xcf32  tw1, tw2, tw3;
    xb_vecN_4xcf32 x0, x1, x2, x3;
    xb_vecN_4xcf32 tmp2;

    NASSERT(M >= 4);
    WUR_CBEGIN((uintptr_t)tw);
    WUR_CEND(N / v*sizeof(*x) + (uintptr_t)tw);


    for (i = 0; i < N / (4 * vect_sz) - 1; i++) //for (i = 0; i < num_bfls; i++)
    {
        /* 9 cycles per pipeline stage in steady state with unroll=1 */

        BBE_LVN_4XCF32_IC(tmp2, tw);
        tw1 = BBE_REPN_4XCF32(tmp2, 0);
        tw2 = BBE_REPN_4XCF32(tmp2, 1);
        tw3 = BBE_REPN_4XCF32(tmp2, 2);

        x0 = BBE_LVN_4XCF32_X(px0, offset_x);
        x1 = BBE_LVN_4XCF32_X(px1, offset_x);
        x2 = BBE_LVN_4XCF32_X(px2, offset_x);
        x3 = BBE_LVN_4XCF32_X(px3, offset_x);

#ifdef IS_INV_FFTF
        IDFT4_FP(x0, x1, x2, x3);

        x1 = BBE_MULJN_4XCF32(x1, tw1);
        x2 = BBE_MULJN_4XCF32(x2, tw2);
        x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
        DFT4_FP(x0, x1, x2, x3);

        x1 = BBE_MULN_4XCF32(x1, tw1);
        x2 = BBE_MULN_4XCF32(x2, tw2);
        x3 = BBE_MULN_4XCF32(x3, tw3);
#endif

        BBE_SVN_4XCF32_X(x0, py0, offset_y);
        BBE_SVN_4XCF32_X(x1, py1, offset_y);
        BBE_SVN_4XCF32_X(x2, py2, offset_y);
        BBE_SVN_4XCF32_X(x3, py3, offset_y);

        offset_x = BBE_ADDMOD16U(offset_x, mod_add_x);
        offset_y = BBE_ADDMOD16U(offset_y, mod_add_y);
    }

    offset_x = N / 4 * sizeof(*x) - 2 * BBE_SIMD_WIDTH;

    BBE_LVN_4XCF32_IC(tmp2, tw);
    tw1 = BBE_REPN_4XCF32(tmp2, 0);
    tw2 = BBE_REPN_4XCF32(tmp2, 1);
    tw3 = BBE_REPN_4XCF32(tmp2, 2);

    x0 = BBE_LVN_4XCF32_X(px0, offset_x);
    x1 = BBE_LVN_4XCF32_X(px1, offset_x);
    x2 = BBE_LVN_4XCF32_X(px2, offset_x);
    x3 = BBE_LVN_4XCF32_X(px3, offset_x);

#ifdef IS_INV_FFTF
    IDFT4_FP(x0, x1, x2, x3);

    x1 = BBE_MULJN_4XCF32(x1, tw1);
    x2 = BBE_MULJN_4XCF32(x2, tw2);
    x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
    DFT4_FP(x0, x1, x2, x3);

    x1 = BBE_MULN_4XCF32(x1, tw1);
    x2 = BBE_MULN_4XCF32(x2, tw2);
    x3 = BBE_MULN_4XCF32(x3, tw3);
#endif

    BBE_SVN_4XCF32_X(x0, py0, offset_y);
    BBE_SVN_4XCF32_X(x1, py1, offset_y);
    BBE_SVN_4XCF32_X(x2, py2, offset_y);
    BBE_SVN_4XCF32_X(x3, py3, offset_y);

    __Pragma("no_reorder");
    return 0;
} //stage_inner_merged

/* Third stage radix-4, used in the cfftf16K, cfftf32K, 
   twiddles format same as stage_inner_merged */
#ifdef IS_INV_FFTF
inline_ int stage_inner_iDFT4_v16_FP
#else
inline_ int stage_inner_DFT4_v16_FP
#endif
   (const complex_float *_tw,
    complex_float *x   /*input*/,
    complex_float *y   /*output*/,
    const int N,
    const int v /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    )
{

    int i, j;
    xb_vecN_4xcf32 *tw = (xb_vecN_4xcf32 *)_tw;

    xb_vecN_4xcf32 *px = (xb_vecN_4xcf32 *)x;
    xb_vecN_4xcf32 *py = (xb_vecN_4xcf32 *)y;
    const int stride_bytes = N / 4 * sizeof(*x);
    const int num_bfls = N / (4 * v);
    int M = v / (BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x)); // Number of iterations in the inner loop

    WUR_CBEGIN((uintptr_t)tw);
    WUR_CEND(N*sizeof(*x) / v + (uintptr_t)tw);

    for (j = 0; j < M; j++)
    {
        py = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*j + (uintptr_t)y);
        px = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*j + stride_bytes + (uintptr_t)x);

        for (i = 0; i < num_bfls; i++)
        {
            xb_vecN_4xcf32  tw1, tw2, tw3;
            xb_vecN_4xcf32 x0, x1, x2, x3;
            xb_vecN_4xcf32 tmp2;
            /* 9 cycles per pipeline stage in steady state with unroll=1 */

            BBE_LVN_4XCF32_IC(tmp2, tw);

            tw1 = BBE_REPN_4XCF32(tmp2, 0);
            tw2 = BBE_REPN_4XCF32(tmp2, 1);
            tw3 = BBE_REPN_4XCF32(tmp2, 2);

            //  px = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*(j + i*M) + (uintptr_t)x);
            //  py = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*(j + i*4*M) + (uintptr_t)y);

            BBE_LVN_4XCF32_XP(x1, px, 2 * stride_bytes);
            BBE_LVN_4XCF32_XP(x3, px, -3 * stride_bytes);
            BBE_LVN_4XCF32_XP(x0, px, 2 * stride_bytes);
            BBE_LVN_4XCF32_XP(x2, px, -1 * stride_bytes + M * 2 * BBE_SIMD_WIDTH);
#ifdef IS_INV_FFTF
            IDFT4_FP(x0, x1, x2, x3);
            x1 = BBE_MULJN_4XCF32(x1, tw1);
            x2 = BBE_MULJN_4XCF32(x2, tw2);
            x3 = BBE_MULJN_4XCF32(x3, tw3);
#else
            DFT4_FP(x0, x1, x2, x3);
            x1 = BBE_MULN_4XCF32(x1, tw1);
            x2 = BBE_MULN_4XCF32(x2, tw2);
            x3 = BBE_MULN_4XCF32(x3, tw3);
#endif
            BBE_SVN_4XCF32_IP(x0, py, 16 * sizeof(*x));
            BBE_SVN_4XCF32_IP(x1, py, 16 * sizeof(*x));
            BBE_SVN_4XCF32_IP(x2, py, 16 * sizeof(*x));
            BBE_SVN_4XCF32_IP(x3, py, 16 * sizeof(*x));

        } //for (i = 0; i < num_bfls; i++)
    } //for (j = 0; j<M; j++)

    __Pragma("no_reorder");

    return 0;
} //stage_inner_DFT4_v16_FP

/* The in-place complex-to-real spectrum conversion */
inline_ void rifft_spec_conv_fp(complex_float *x,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N)
{
#if 1

    const int v_sz = (BBE_SIMD_WIDTH * 2 / sizeof(*x)); //Vector size (number of complex elements)
    xb_vecN_4xcf32 _a0, _a1, _b0, _b1, tw;

    xb_vecN_4xcf32 *ptw = (xb_vecN_4xcf32 *)(twiddle_table + N / 4 - BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(complex_float));
    xb_vecN_4xcf32 *px0 = (xb_vecN_4xcf32 *)(x + N / 4);
    xb_vecN_4xcf32 *py0 = (xb_vecN_4xcf32 *)(x + N / 4 - v_sz);
    xb_vecN_4xcf32 *px1 = (xb_vecN_4xcf32 *)(x + N / 4 + 1);
    xb_vecN_4xcf32 *py1 = (xb_vecN_4xcf32 *)(x + N / 4 + 1);
    valign vx1;
    valign vy1 = BBE_ZALIGN();
    xtcomplexfloat yN_4; // y[N/4]

    /* y[N / 4] = conj_fl32c(x[N / 4]); */
    BBE_LVN_4XCF32_XP(_a0, px0, -2 * BBE_SIMD_WIDTH);
    yN_4 = BBE_SELSN_4XCF32(_a0, 0);
    yN_4 = BBE_CONJCF32(yN_4);

    vx1 = BBE_LAN_4XCF32_PP(px1);

    int n;

    /* Calculation of a x[0] included into loop */
    //a0.s.re = x[0].s.re + x[N / 2].s.re;
    //a0.s.im = x[0].s.re - x[N / 2].s.re;

    /* Divide spectrum by 2 so that overall scale factor of RIFFT
    * is 1/2*2/N == N. */
    //  x[0].s.re = ldexpf(a0.s.re, -1);
    //  x[0].s.im = ldexpf(a0.s.im, -1);

    for (n = N / 4 - v_sz; n >= 0; n -= v_sz)
    {
        BBE_LVN_4XCF32_IP(_a0, px0, -2 * BBE_SIMD_WIDTH);
        BBE_LAN_4XCF32_IP(_a1, vx1, px1);
        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);

        _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);

        BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, BBE_CONJN_4XCF32(_a1));


        _b1 = BBE_MULJN_4XCF32(_b1, tw);
        _SCALE_N_4XCF32(_b0, 1.0f / 2);

        BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);
        _a1 = BBE_CONJN_4XCF32(_a1);
        _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
        BBE_SVN_4XCF32_IP(_a0, py0, -2 * BBE_SIMD_WIDTH);
        BBE_SAN_4XCF32_IP(_a1, vy1, py1);
    }
    //  Flush not needed, since no write to x[N/2] !
    //  BBE_SAN_4XCF32POS_FP(vy1, py1); 
    /* conj(x(N/4+1)) */
    //x[N / 4] = conj_fl32c(x[N / 4]);

    *((xtcomplexfloat*)x + N / 4) = yN_4;

#else
    // reference code
    complex_float a0, a1;
    complex_float b0, b1;

    int n;

    /* Calculation of a x[0] included into loop */
    //a0.s.re = x[0].s.re + x[N / 2].s.re;
    //a0.s.im = x[0].s.re - x[N / 2].s.re;

    /* Divide spectrum by 2 so that overall scale factor of RIFFT
    * is 1/2*2/N == N. */
    //  x[0].s.re = ldexpf(a0.s.re, -1);
    //  x[0].s.im = ldexpf(a0.s.im, -1);

    for (n = 0; n<N / 4; n++)
        //for (n = N / 4 - 1; n >= 0; n--) // works too
    {

        a0 = x[n];
        a1 = x[(N / 2 - n)];

        /* b0 =  a0+conj(a1); */
        b0.s.re = a0.s.re + a1.s.re;
        b0.s.im = a0.s.im - a1.s.im;
        /* b1 =  a0-conj(a1); */
        b1.s.re = a0.s.re - a1.s.re;
        b1.s.im = a0.s.im + a1.s.im;

        /* b1 <- b1*conj( -j *0.5*exp(-j*2*PI/N*n) ); */
        b1 = mul_fl32c(b1, conj_fl32c(twiddle_table[n]));

        b0.s.re = ldexpf(b0.s.re, -1);
        b0.s.im = ldexpf(b0.s.im, -1);


        /* a0 <- b0 + b1;             */
        a0.s.re = b0.s.re + b1.s.re;
        a0.s.im = b0.s.im + b1.s.im;
        /* a1 <- conj( b0 - conj(b1) );*/
        a1.s.re = b0.s.re - b1.s.re;
        a1.s.im = -b0.s.im + b1.s.im;

        x[n] = a0;

        if (n != 0)
            x[(N / 2 - n)] = a1;
    }

    /* conj(x(N/4+1)) */
    x[N / 4] = conj_fl32c(x[N / 4]);
#endif
}


/* The in-place real-to-complex spectrum conversion */
inline_ void rfft_spec_conv_fp(complex_float *y,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N)
{

    valign u, u1;
    int n;

    const int v_sz = (BBE_SIMD_WIDTH * 2 / sizeof(*y)); //Vector size (number of complex elements)
    xb_vecN_4xcf32 *px0 = (xb_vecN_4xcf32 *)(y + N / 4);
    xb_vecN_4xcf32 *px1 = (xb_vecN_4xcf32 *)(y + N / 4 + 1);
    xb_vecN_4xcf32 *py1 = (xb_vecN_4xcf32 *)(y + N / 4 + 1);
    xb_vecN_4xcf32 *py0 = (xb_vecN_4xcf32 *)(y + N / 4 - BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*y));
    xb_vecN_4xcf32 *ptw = (xb_vecN_4xcf32 *)(twiddle_table + N / 4 - BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(complex_float));

    u1 = BBE_ZALIGN();

    WUR_CBEGIN((uintptr_t)y);
    WUR_CEND((uintptr_t)(y + (N / 2)));
    BBE_LAN_4XCF32POS_PC(u, px1);

    NASSERT_ALIGN32(y);

    {

#if 1
        xb_vecN_4xcf32 _a0, _a1, _b0, _b1, tw;
        xtcomplexfloat yN_4; // y[N/4]

        /* y[N / 4] = conj_fl32c(y[N / 4]); */
        BBE_LVN_4XCF32_XP(_a0, px0, -2 * BBE_SIMD_WIDTH);
        yN_4 = BBE_SELSN_4XCF32(_a0, 0);
        yN_4 = BBE_CONJCF32(yN_4);
        /* 10 cycles per pipeline stage in steady state with unroll=2 */
        for (n = 0; n < N / (4 * v_sz); n++)
        {
            BBE_LVN_4XCF32_XP(_a0, px0, -2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IC(_a1, u, px1);
            BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw);
            xb_vecN_2xf32 tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);
            BBE_SVN_4XCF32_XP(_a0, py0, -2 * BBE_SIMD_WIDTH);
        }
        BBE_SAN_4XCF32POS_FP(u1, py1);

        *((xtcomplexfloat*)y + N / 4) = yN_4;

#else

        // Reference code
        complex_float a0, a1;
        complex_float b0, b1;

        a0.s.re = y[0].s.re + y[0].s.im; a0.s.im = 0;
        a1.s.re = y[0].s.re - y[0].s.im; a1.s.im = 0;

        y[0] = a0;
        y[N / 2] = conj_fl32c(a1);


        for (n = N / 4 - 1; n > 0; n--)
        {
            complex_float conj_a1;

            a0 = y[n];
            a1 = y[(N / 2 - n) & (N / 2 - 1)]; // Circular addresing

            /* b0 <- 1/2*(a0+conj(a1)); */
            /* b1 <- 1/2*(a0-conj(a1))*-1j */
            /* b1 <- b1*twd */

            conj_a1 = conj_fl32c(a1);

            b0.s.re = a0.s.re + conj_a1.s.re;
            b0.s.im = a0.s.im + conj_a1.s.im;

            b1.s.re = a0.s.re - conj_a1.s.re;
            b1.s.im = a0.s.im - conj_a1.s.im;

            b0.s.re *= 0.5f;
            b0.s.im *= 0.5f;
            b1 = mul_fl32c(twiddle_table[n], b1);

            /* a0 <- b0+b1 */
            a0.s.re = b0.s.re + b1.s.re;
            a0.s.im = b0.s.im + b1.s.im;
            /* a1 <- b0-b1 */
            a1.s.re = b0.s.re - b1.s.re;
            a1.s.im = b0.s.im - b1.s.im;

            y[n] = a0;
            y[N / 2 - n] = conj_fl32c(a1);
        }
        y[N / 4] = conj_fl32c(y[N / 4]);
#endif
    }//rfft_spec_conv_fp
}


/* first phase of the tfftN  */
inline_ int stage_first_tfft_DFT4_FP(const complex_float *tw,
    complex_float *x,
    complex_float *y,
    int N)
{

    int i;
    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);

    xb_vecN_4xcf32 t0, t1, t2, t3, tw1, tw2, tw3;
    xb_vecN_4xcf32 rotator;
    valign u = BBE_LAN_4XCF32_PP(p_tw);

    /* Load rotator exp(1j*2*pi/N *(0:v-1) , v is vector size 2*BBE_SIMD_WIDTH/sizeof(complex_float) */
    BBE_LAVN_4XCF32_XP(rotator, u, p_tw, 2 * BBE_SIMD_WIDTH);
    for (i = 0; i<count; i++)
    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
        BBE_LVN_4XCF32_XP(_x0, p_src, stride);
        BBE_LVN_4XCF32_XP(_x1, p_src, stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + 2 * BBE_SIMD_WIDTH);

        /* 35 cycles per pipeline stage in steady state with unroll = 2 */
        BBE_LAVN_4XCF32_XP(tw1, u, p_tw, sizeof(complex_float));

        tw1 = BBE_REPN_4XCF32(tw1, 0);
        tw1 = BBE_MULN_4XCF32(tw1, rotator);
        tw2 = BBE_MULN_4XCF32(tw1, tw1);
        tw3 = BBE_MULN_4XCF32(tw2, tw1);


        {
            xb_vecN_4xcf32 d0, d1, s0, s1;

            BBE_ADDSUBN_4XCF32(d0, s0, _x0, _x2);
            BBE_ADDSUBN_4XCF32(d1, s1, _x1, _x3);
            d1 = BBE_SHFLN_4XCF32I(d1, BBE_SHFLI_SWAP_2);
            d1 = BBE_CONJN_4XCF32(d1);
            BBE_ADDSUBN_4XCF32(_x2, _x0, s0, s1);
            BBE_ADDSUBN_4XCF32(_x3, _x1, d0, d1);
        }


        _x1 = BBE_MULN_4XCF32(_x1, tw1);
        _x2 = BBE_MULN_4XCF32(_x2, tw2);
        _x3 = BBE_MULN_4XCF32(_x3, tw3);

        BBE_DSELN_4XCF32I(t1, t0, _x2, _x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(t3, t2, _x3, _x1, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x1, _x0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x3, _x2, t3, t1, BBE_DSELI_INTERLEAVE_4);

        BBE_SVN_4XCF32_IP(_x0, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x1, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x2, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x3, p_dst, 2 * BBE_SIMD_WIDTH);

    }
    __Pragma("no_reorder");
    return 0;
}

/* First phase of the tifft  */
inline_ int stage_first_tifft_IDFT4_FP(const complex_float *tw, complex_float *x, complex_float *y, int N)
{
    int i;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);

    xb_vecN_4xcf32 t0, t1, t2, t3, tw1, tw2, tw3;
    xb_vecN_4xcf32 rotator;
    valign u = BBE_LAN_4XCF32_PP(p_tw);


    /* Load rotator exp(1j*2*pi/N *(0:v-1) , v is vector size 2*BBE_SIMD_WIDTH/sizeof(complex_float) */
    BBE_LAVN_4XCF32_XP(rotator, u, p_tw, 2 * BBE_SIMD_WIDTH);
    for (i = 0; i<count; i++)
    {
        xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
        BBE_LVN_4XCF32_XP(_x0, p_src, stride);
        BBE_LVN_4XCF32_XP(_x1, p_src, stride);
        BBE_LVN_4XCF32_XP(_x2, p_src, stride);
        BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + 2 * BBE_SIMD_WIDTH);

        BBE_LAVN_4XCF32_XP(tw1, u, p_tw, sizeof(complex_float));

        tw1 = BBE_REPN_4XCF32(tw1, 0);
        tw1 = BBE_MULN_4XCF32(tw1, rotator);
        tw2 = BBE_MULN_4XCF32(tw1, tw1);
        tw3 = BBE_MULN_4XCF32(tw2, tw1);

        IDFT4_FP(_x0, _x1, _x2, _x3);

        _x1 = BBE_MULJN_4XCF32(_x1, tw1);
        _x2 = BBE_MULJN_4XCF32(_x2, tw2);
        _x3 = BBE_MULJN_4XCF32(_x3, tw3);

        BBE_DSELN_4XCF32I(t1, t0, _x2, _x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(t3, t2, _x3, _x1, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x1, _x0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(_x3, _x2, t3, t1, BBE_DSELI_INTERLEAVE_4);

        BBE_SVN_4XCF32_IP(_x0, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x1, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x2, p_dst, 2 * BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(_x3, p_dst, 2 * BBE_SIMD_WIDTH);

    }

    return 0;
}



/* Last stage radix 8 used in the bcfft32,128 */
inline_ int blk_last_stage_DFT8_FP(complex_float *x,
    complex_float *y,
    int N, int L)
{


    int i, l;
    const int stride = N / 8 * sizeof(*x);
    int count = N / (8 * BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x));

    xb_vecN_4xcf32 * p_src0 = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_src1 = (xb_vecN_4xcf32 *)(x + stride / sizeof(*x));
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    xb_vecN_4xcf32 a0, a1, b0, b1, c0, c1, d0, d1;


    for (i = 0; i<count; i++)
    {

        p_src0 = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
        p_src1 = (xb_vecN_4xcf32*)(stride + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
        p_dst = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*i + (uintptr_t)y);
        /* 16 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {
            //    p_src0 = (xb_vecN_4xcf32*)(N*l*sizeof(*x) + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
            //    p_src1 = (xb_vecN_4xcf32*)(N*l*sizeof(*x) + stride + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
            //    p_dst = (xb_vecN_4xcf32*)(N*l *sizeof(*x) + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)y);

            BBE_LVN_4XCF32_XP(a0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(b0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(c0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(d0, p_src0, -6 * stride + N*sizeof(*x));

            BBE_LVN_4XCF32_XP(a1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(b1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(c1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(d1, p_src1, -6 * stride + N*sizeof(*x));

            DFT4_FP(a0, b0, c0, d0);
            {
                xb_vecN_4xcf32 _d0, _d1, _s0, _s1;
                xb_vecN_4xcf32 dd, sd;
                BBE_ADDSUBN_4XCF32(_d0, _s0, a1, c1);
                BBE_ADDSUBN_4XCF32(_d1, _s1, b1, d1);

                BBE_ADDSUBN_4XCF32(dd, sd, _d0, _d1);
                // sd *= -j
                sd = BBE_SHFLN_4XCF32I(sd, BBE_SHFLI_SWAP_2);
                sd = BBE_CONJN_4XCF32(sd);

                BBE_ADDSUBN_4XCF32(d1, b1, sd, dd);

                _SCALE_N_4XCF32(b1, 0.707106781186547f);
                _SCALE_N_4XCF32(d1, 0.707106781186547f);

                BBE_ADDSUBN_4XCF32(c1, a1, _s0, _s1);
            }

            // c1 *= 1j; 
            c1 = BBE_SHFLN_4XCF32I(c1, BBE_SHFLI_SWAP_2);
            c1 = BBE_CONJN_4XCF32(c1);


            BBE_ADDSUBN_4XCF32(a1, a0, a0, a1);
            BBE_ADDSUBN_4XCF32(b1, b0, b0, b1);
            BBE_ADDSUBN_4XCF32(c1, c0, c0, c1);
            BBE_ADDSUBN_4XCF32(d1, d0, d0, d1);



            BBE_SVN_4XCF32_XP(a0, p_dst, stride);
            BBE_SVN_4XCF32_XP(b0, p_dst, stride);
            BBE_SVN_4XCF32_XP(c0, p_dst, stride);
            BBE_SVN_4XCF32_XP(d0, p_dst, stride);

            BBE_SVN_4XCF32_XP(a1, p_dst, stride);
            BBE_SVN_4XCF32_XP(b1, p_dst, stride);
            BBE_SVN_4XCF32_XP(c1, p_dst, stride);
            BBE_SVN_4XCF32_XP(d1, p_dst, -7 * stride + N*sizeof(*x));

        } //for (l = 0; l < L; l++)
    }//for (i = 0; i<count; i++)

    __Pragma("no_reorder");
    return 0;
} //blk_last_stage_DFT8_FP

/* Last stage of the bcifftf32, 128 */
inline_ int blk_last_stage_iDFT8_FP(complex_float *x, complex_float *y, int N, int L)
{


    int i, l;
    const int stride = N / 8 * sizeof(*x);
    int count = N / (8 * BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x));

    xb_vecN_4xcf32 * p_src0 = (xb_vecN_4xcf32 *)(x);
    xb_vecN_4xcf32 * p_src1 = (xb_vecN_4xcf32 *)(x + stride / sizeof(*x));
    xb_vecN_4xcf32 * p_dst = (xb_vecN_4xcf32 *)(y);
    xb_vecN_4xcf32 a0, a1, b0, b1, c0, c1, d0, d1;



    for (i = 0; i<count; i++)
    {

        p_src0 = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
        p_src1 = (xb_vecN_4xcf32*)(stride + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
        p_dst = (xb_vecN_4xcf32*)(2 * BBE_SIMD_WIDTH*i + (uintptr_t)y);
        /* 16 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {
            //    p_src0 = (xb_vecN_4xcf32*)(N*l*sizeof(*x) + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
            //    p_src1 = (xb_vecN_4xcf32*)(N*l*sizeof(*x) + stride + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)x);
            //    p_dst = (xb_vecN_4xcf32*)(N*l *sizeof(*x) + 2 * BBE_SIMD_WIDTH*i + (uintptr_t)y);

            BBE_LVN_4XCF32_XP(a0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(b0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(c0, p_src0, 2 * stride);
            BBE_LVN_4XCF32_XP(d0, p_src0, -6 * stride + N*sizeof(*x));

            BBE_LVN_4XCF32_XP(a1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(b1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(c1, p_src1, 2 * stride);
            BBE_LVN_4XCF32_XP(d1, p_src1, -6 * stride + N*sizeof(*x));

            DFT4_FP(a0, b0, c0, d0);
            {
                xb_vecN_4xcf32 _d0, _d1, _s0, _s1;
                xb_vecN_4xcf32 dd, sd;
                BBE_ADDSUBN_4XCF32(_d0, _s0, a1, c1);
                BBE_ADDSUBN_4XCF32(_d1, _s1, b1, d1);

                BBE_ADDSUBN_4XCF32(dd, sd, _d0, _d1);
                // sd *= -j
                sd = BBE_SHFLN_4XCF32I(sd, BBE_SHFLI_SWAP_2);
                sd = BBE_CONJN_4XCF32(sd);

                BBE_ADDSUBN_4XCF32(d1, b1, sd, dd);

                _SCALE_N_4XCF32(b1, 0.707106781186547f / N);
                _SCALE_N_4XCF32(d1, 0.707106781186547f / N);


                BBE_ADDSUBN_4XCF32(c1, a1, _s0, _s1);
            }

            // c1 *= 1j; 
            c1 = BBE_SHFLN_4XCF32I(c1, BBE_SHFLI_SWAP_2);
            c1 = BBE_CONJN_4XCF32(c1);


            _SCALE_N_4XCF32(a0, 1.0f / N);
            _SCALE_N_4XCF32(d0, 1.0f / N);
            _SCALE_N_4XCF32(c0, 1.0f / N);

            _SCALE_N_4XCF32(a1, 1.0f / N);

            _SCALE_N_4XCF32(c1, 1.0f / N);
            _SCALE_N_4XCF32(b0, 1.0f / N);

            BBE_ADDSUBN_4XCF32(a1, a0, a0, a1);
            BBE_ADDSUBN_4XCF32(b1, b0, b0, b1);
            BBE_ADDSUBN_4XCF32(c1, c0, c0, c1);
            BBE_ADDSUBN_4XCF32(d1, d0, d0, d1);

            BBE_SVN_4XCF32_XP(a0, p_dst, stride);
            BBE_SVN_4XCF32_XP(d1, p_dst, stride);
            BBE_SVN_4XCF32_XP(c1, p_dst, stride);
            BBE_SVN_4XCF32_XP(b1, p_dst, stride);

            BBE_SVN_4XCF32_XP(a1, p_dst, stride);
            BBE_SVN_4XCF32_XP(d0, p_dst, stride);
            BBE_SVN_4XCF32_XP(c0, p_dst, stride);
            BBE_SVN_4XCF32_XP(b0, p_dst, -7 * stride + N*sizeof(*x));

        } //for (l = 0; l < L; l++)
    }//for (i = 0; i<count; i++)

    __Pragma("no_reorder");
    return 0;
} /* blk_last_stage_iDFT8_FP */

/* first stage radix-4 of bcfft  */
inline_ int blk_first_stage_DFT4_FP(
    const complex_float *tw,
    complex_float *x,
    complex_float *y,
    int N,
    int L)
{
    int i, l;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);
    xb_vecN_4xcf32  tw1, tw2, tw3;
    valign uu0, uu1;
    xb_vecNx16 t0, t1, t2, t3;

    for (i = 0; i<count; i++)
    {
        p_src = (xb_vecN_4xcf32 *)(x + 4 * i);
        p_dst = (xb_vecNx16 *)(y + 16 * i);

        BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

        for (l = 0; l < L; l++)
        {
            xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
            /* 11 cycles per pipeline stage in steady state with unroll=1 */
            BBE_LVN_4XCF32_XP(_x0, p_src, stride);
            BBE_LVN_4XCF32_XP(_x1, p_src, stride);
            BBE_LVN_4XCF32_XP(_x2, p_src, stride);
            BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + N*sizeof(*x));

            DFT4_FP(_x0, _x1, _x2, _x3);

            _x1 = BBE_MULN_4XCF32(_x1, tw1);
            _x2 = BBE_MULN_4XCF32(_x2, tw2);
            _x3 = BBE_MULN_4XCF32(_x3, tw3);

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
            BBE_SALIGNVRNX16_XP(t3, uu1, p_dst, -6 * BBE_SIMD_WIDTH + N*sizeof(*y));

        } //  for (l = 0; l < L; l++)
    } // for (i = 0; i<count; i++)

    __Pragma("no_reorder");
    return 0;
} //blk_first_stage_DFT4_FP

/* First stage IDFT4 for bcifftN */
inline_ int blk_first_stage_iDFT4_FP(
    const complex_float *tw,
    complex_float *x,
    complex_float *y,
    int N,
    int L)
{
    int i, l;

    int count = N / 4 / (BBE_SIMD_WIDTH * sizeof(int16_t) / sizeof(*x));
    int stride = N / 4 * sizeof(*x);
    xb_vecN_4xcf32 * p_tw = (xb_vecN_4xcf32 *)(tw);
    xb_vecN_4xcf32 * p_src = (xb_vecN_4xcf32 *)(x);
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);
    xb_vecN_4xcf32  tw1, tw2, tw3;
    valign uu0, uu1;
    xb_vecNx16 t0, t1, t2, t3;



    for (i = 0; i<count; i++)
    {
        p_src = (xb_vecN_4xcf32 *)(x + 4 * i);
        p_dst = (xb_vecNx16 *)(y + 16 * i);

        BBE_LVN_4XCF32_IP(tw1, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw2, p_tw, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(tw3, p_tw, 2 * BBE_SIMD_WIDTH);

        for (l = 0; l < L; l++)
        {
            xb_vecN_4xcf32 _x0, _x1, _x2, _x3;
            /* 11 cycles per pipeline stage in steady state with unroll=1 */
            BBE_LVN_4XCF32_XP(_x0, p_src, stride);
            BBE_LVN_4XCF32_XP(_x1, p_src, stride);
            BBE_LVN_4XCF32_XP(_x2, p_src, stride);
            BBE_LVN_4XCF32_XP(_x3, p_src, -3 * stride + N*sizeof(*x));

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
            BBE_SALIGNVRNX16_XP(t3, uu1, p_dst, -6 * BBE_SIMD_WIDTH + N*sizeof(*y));

        } //  for (l = 0; l < L; l++)
    } // for (i = 0; i<count; i++)

    __Pragma("no_reorder");
    return 0;
} /* blk_first_stage_iDFT4_FP */


/* Inner phase bcifftN  N=64, 128 */
inline_ void blk_inner_stage_iDFT4_FP(const complex_float *_tw,
    complex_float *x        /*input*/,
    complex_float *y        /*output*/,
    const int N,
    const int L             /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    )
{
    int i, l;
    const int v = 4;
    xb_vecN_4xcf32 *tw = (xb_vecN_4xcf32 *)_tw;

    xb_vecN_4xcf32 *px = (xb_vecN_4xcf32 *)x;
    xb_vecN_4xcf32 *py = (xb_vecN_4xcf32 *)y;

    valign u = BBE_LAVNX16_PP((xb_vecNx16*)tw);

    const int stride_bytes = N / 4 * sizeof(*x);
    const int num_bfls = N / (4 * v);

    NASSERT(BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x) == v);


    for (i = 0; i < num_bfls; i++)
    {
        xb_vecN_4xcf32 tw1, tw2, tw3;
        xb_vecN_4xcf32 tmp2;

        BBE_LAVN_4XCF32_XP(tmp2, u, tw, 3 * sizeof(*x));

        tw1 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_0X4);
        tw2 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_1X4);
        tw3 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_2X4);
        px = (xb_vecN_4xcf32 *)(x + i * 4);
        py = (xb_vecN_4xcf32 *)(y + i * 16);
        /* 8 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {

            xb_vecN_4xcf32 x0, x1, x2, x3;

            //px = (xb_vecN_4xcf32 *)(x + l*N + i * 4 );
            //py = (xb_vecN_4xcf32 *)(y + l*N + i * 16);

            BBE_LVN_4XCF32_XP(x0, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x1, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x2, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x3, px, stride_bytes);

            IDFT4_FP(x0, x1, x2, x3);

            x1 = BBE_MULJN_4XCF32(x1, tw1);
            x2 = BBE_MULJN_4XCF32(x2, tw2);
            x3 = BBE_MULJN_4XCF32(x3, tw3);

            BBE_SVN_4XCF32_IP(x0, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x1, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x2, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_XP(x3, py, N*sizeof(complex_float)-6 * BBE_SIMD_WIDTH);
        }

    }

    __Pragma("no_reorder");
} //blk_inner_stage_iDFT4_FP




/* Inner stage of the bcfftf64, bcfftf128 */
inline_ void blk_inner_stage_DFT4_FP(const complex_float *_tw,
    complex_float *x   /*input*/,
    complex_float *y   /*output*/,
    const int N,
    const int L /*vector length must be multiple of BBE_SIMD_WIDTH/4 */
    )
{
    int i, l;
    const int v = 4;
    xb_vecN_4xcf32 *tw = (xb_vecN_4xcf32 *)_tw;

    xb_vecN_4xcf32 *px = (xb_vecN_4xcf32 *)x;
    xb_vecN_4xcf32 *py = (xb_vecN_4xcf32 *)y;

    valign u = BBE_LAVNX16_PP((xb_vecNx16*)tw);

    const int stride_bytes = N / 4 * sizeof(*x);
    const int num_bfls = N / (4 * v);

    NASSERT(BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(*x) == v);


    for (i = 0; i < num_bfls; i++)
    {
        xb_vecN_4xcf32 tw1, tw2, tw3;
        xb_vecN_4xcf32 tmp2;

        BBE_LAVN_4XCF32_XP(tmp2, u, tw, 3 * sizeof(*x));

        tw1 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_0X4);
        tw2 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_1X4);
        tw3 = BBE_SHFLN_4XCF32I(tmp2, BBE_SHFLI_REP_2X4);
        px = (xb_vecN_4xcf32 *)(x + i * 4);
        py = (xb_vecN_4xcf32 *)(y + i * 16);
        /* 8 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {

            xb_vecN_4xcf32 x0, x1, x2, x3;

            //px = (xb_vecN_4xcf32 *)(x + l*N + i * 4 );
            //py = (xb_vecN_4xcf32 *)(y + l*N + i * 16);

            BBE_LVN_4XCF32_XP(x0, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x1, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x2, px, stride_bytes);
            BBE_LVN_4XCF32_XP(x3, px, stride_bytes);

            DFT4_FP(x0, x1, x2, x3);

            x1 = BBE_MULN_4XCF32(x1, tw1);
            x2 = BBE_MULN_4XCF32(x2, tw2);
            x3 = BBE_MULN_4XCF32(x3, tw3);

            BBE_SVN_4XCF32_IP(x0, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x1, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_IP(x2, py, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_XP(x3, py, N*sizeof(complex_float)-6 * BBE_SIMD_WIDTH);
        }

    }

    __Pragma("no_reorder");
} /* blk_inner_stage_DFT4_FP */

inline_ void blkrfft_spec_conv_packed_fp(
                            complex_float *_y,   /* in/out  complex array, size is N/2 */
                            const complex_float *twiddle_table,
                            int N, int L) ATTRIBUTE_ALWAYS_INLINE;

/* The in-place real-to-complex spectrum conversion, N>16 */
inline_ void blkrfft_spec_conv_packed_fp(
                complex_float *_y,   /* in/out  complex array, size is N/2 */
                const complex_float *twiddle_table,
                int N,
                int L)
{
    valign u, u1;
    int l, n;
    xb_vecN_4xcf32 *px0;
    xb_vecN_4xcf32 *px1;
    xb_vecN_4xcf32 *py1;
    xb_vecN_4xcf32 *py0;

    const int v_sz = (BBE_SIMD_WIDTH * 2 / sizeof(*_y)); //Vector size (number of complex elements)
    NASSERT_ALIGN32(_y);
    vboolN_4 mask = BBE_NOTBN_4(BBE_MOVN_4_FROMN(BBE_LTRNI(4)));

    xb_vecN_4xcf32 _a0, _a1, _b0, _b1, tw;
    xb_vecN_4xcf32 *ptw = (xb_vecN_4xcf32 *)(twiddle_table + N / 4 - v_sz);

    NASSERT(N>16); 

    for (n = 0; n < N / (4 * v_sz) - 2; n += 2)
    {
        xb_vecN_4xcf32 tw2;
        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(tw2, ptw, -2 * BBE_SIMD_WIDTH);

        px0 = (xb_vecN_4xcf32 *)(_y - v_sz*(n + 1) + N / 4);
        py0 = (xb_vecN_4xcf32 *)(_y - v_sz*(n + 1) + N / 4);
        px1 = (xb_vecN_4xcf32 *)(_y + v_sz*n + N / 4 + 1);
        py1 = (xb_vecN_4xcf32 *)(_y + v_sz*n + N / 4 + 1);

        /* 10 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l<L; l++)
        {  
             /*
             ptw = (xb_vecN_4xcf32 *)(twiddle_table - v_sz*n + N / 4 - BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(complex_float));
             px0 = (xb_vecN_4xcf32 *)((_y + N / 2 * l) - v_sz*(n+1) + N / 4);
             px1 = (xb_vecN_4xcf32 *)((_y + N / 2 * l) + v_sz*n + N / 4 + 1);
             py0 = (xb_vecN_4xcf32 *)((_y + N / 2 * l) - v_sz*(n + 1) + N / 4);
             py1 = (xb_vecN_4xcf32 *)((_y + N / 2 * l) + v_sz*n + N / 4 + 1);
            */
            u1 = BBE_ZALIGN();

            u = BBE_LAN_4XCF32_PP(px1);

            BBE_LVN_4XCF32_XP(_a0, px0, -2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, u, px1);


            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw);
            xb_vecN_2xf32 tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);
            BBE_SVN_4XCF32_XP(_a0, py0, -2 * BBE_SIMD_WIDTH);

            BBE_LVN_4XCF32_XP(_a0, px0, N / 2 * sizeof(*_y) + 2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, u, px1);


            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw2);

            tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);
            BBE_SVN_4XCF32_XP(_a0, py0, N / 2 * sizeof(*_y) + 2 * BBE_SIMD_WIDTH);

            BBE_SAN_4XCF32POS_FP(u1, py1);
            py1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 4 * BBE_SIMD_WIDTH + (uintptr_t)py1);
            px1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 4 * BBE_SIMD_WIDTH + (uintptr_t)px1);
        }// for (l=0....
    }
    {
        xb_vecN_4xcf32 tw2;
        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(tw2, ptw, -2 * BBE_SIMD_WIDTH);

        px0 = (xb_vecN_4xcf32 *)(_y - v_sz*(n + 1) + N / 4);
        py0 = (xb_vecN_4xcf32 *)(_y - v_sz*(n + 1) + N / 4);
        px1 = (xb_vecN_4xcf32 *)(_y + v_sz*n + N / 4 + 1);
        py1 = (xb_vecN_4xcf32 *)(_y + v_sz*n + N / 4 + 1);

        /* 10 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l<L; l++)
        {

            u1 = BBE_ZALIGN();
            u = BBE_LAN_4XCF32_PP(px1);

            BBE_LVN_4XCF32_XP(_a0, px0, -2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, u, px1);


            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw);
            xb_vecN_2xf32 tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);
            BBE_SVN_4XCF32_XP(_a0, py0, -2 * BBE_SIMD_WIDTH);

            BBE_LVN_4XCF32_XP(_a0, px0, N / 2 * sizeof(*_y) + 2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, u, px1);


            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw2);

            tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);
            BBE_SVN_4XCF32T_XP(_a0, py0, N / 2 * sizeof(*_y) + 2 * BBE_SIMD_WIDTH, mask);

            py1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 4 * BBE_SIMD_WIDTH + (uintptr_t)py1);
            px1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 4 * BBE_SIMD_WIDTH + (uintptr_t)px1);
        }// for (l=0....
    }


    __Pragma("no_reorder");

    px0 = (xb_vecN_4xcf32 *)(_y);
    py0 = (xb_vecN_4xcf32 *)(_y);
    /* 16 cycles per pipeline stage in steady state with unroll = 8 */
    for (l = 0; l<L; l++)
    {
        xb_vecN_4xcf32 a0, aN_4;

        BBE_LVN_4XCF32_XP(a0, px0, N / 4 * sizeof(*_y));
        BBE_LVN_4XCF32_XP(aN_4, px0, N / 4 * sizeof(*_y));
        aN_4 = BBE_CONJN_4XCF32(aN_4);
        a0 = BBE_MULJN_4XCF32(BBE_CONSTN_4XCF32(0x1), a0);
        BBE_SVN_4XCF32F_XP(a0, py0, N / 4 * sizeof(*_y), mask);
        BBE_SVN_4XCF32F_XP(aN_4, py0, N / 4 * sizeof(*_y), mask);
    }

} /* blkrfft_spec_conv_packed_fp */


inline_ void blkrifft_spec_conv_packed_fp(  complex_float *_x,   /* in/out  complex array, size is N/2+1 */
                                            const complex_float *twiddle_table,
                                            int N, int L)ATTRIBUTE_ALWAYS_INLINE;

/* The in-place complex-to-real spectrum conversion, N must be > 16 */
inline_ void blkrifft_spec_conv_packed_fp(complex_float *_x,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N, int L)
{
    int l;
    int n;
    valign vx1;
    valign vy1;
    vboolN_4 m = BBE_NOTBN_4(BBE_MOVN_4_FROMN(BBE_LTRNI(4)));
    // Vector size (number of complex elements)
    const int v_sz = (BBE_SIMD_WIDTH * 2 / sizeof(*_x)); 
    xb_vecN_4xcf32 _a0, _a1, _b0, _b1, tw, tw2;
    xb_vecN_4xcf32 *ptw = (xb_vecN_4xcf32 *)(twiddle_table + N / 4 - v_sz);
    xb_vecN_4xcf32 *px0;
    xb_vecN_4xcf32 *py0;
    xb_vecN_4xcf32 *px1;
    xb_vecN_4xcf32 *py1;

    NASSERT(N > 16);
    NASSERT((N&(N - 1)) == 0);

    for (n = 0; n < N / (4 * v_sz) - 2; n += 2)
    {
        px0 = (xb_vecN_4xcf32 *)(_x - v_sz*(n + 1) + N / 4);
        px1 = (xb_vecN_4xcf32 *)(_x + v_sz*n + N / 4 + 1);
        py0 = (xb_vecN_4xcf32 *)(_x - v_sz*(n + 1) + N / 4);
        py1 = (xb_vecN_4xcf32 *)(_x + v_sz*n + N / 4 + 1);

        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(tw2, ptw, -2 * BBE_SIMD_WIDTH);
        /* 10 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {
            /*
               ptw = (xb_vecN_4xcf32 *)(twiddle_table - v_sz*n + N / 4 - BBE_SIMD_WIDTH*sizeof(int16_t) / sizeof(complex_float));
                px0 = (xb_vecN_4xcf32 *)((_x + N / 2 * l) - v_sz*(n+1) + N / 4);
                px1 = (xb_vecN_4xcf32 *)((_x + N / 2 * l) + v_sz*n + N / 4 + 1);
              py0 = (xb_vecN_4xcf32 *)((_x + N / 2 * l) - v_sz*(n + 1) + N / 4);
                py1 = (xb_vecN_4xcf32 *)((_x + N / 2 * l) + v_sz*n + N / 4 + 1);
            */
            vy1 = BBE_ZALIGN();
            vx1 = BBE_LAN_4XCF32_PP(px1);
            /**************** step 1 ****************/
            BBE_LVN_4XCF32_IP(_a0, px0, -2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, vx1, px1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, BBE_CONJN_4XCF32(_a1));

            _b1 = BBE_MULJN_4XCF32(_b1, tw);
            _SCALE_N_4XCF32(_b0, 1.0f / 2);

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);
            _a1 = BBE_CONJN_4XCF32(_a1);
            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);

            BBE_SVN_4XCF32_IP(_a0, py0, -2 * BBE_SIMD_WIDTH);
            BBE_SAN_4XCF32_IP(_a1, vy1, py1);
            /**************** step 2 ****************/
            BBE_LVN_4XCF32_XP(_a0, px0, N / 2 * sizeof(*_x) + 2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, vx1, px1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, BBE_CONJN_4XCF32(_a1));

            _b1 = BBE_MULJN_4XCF32(_b1, tw2);
            _SCALE_N_4XCF32(_b0, 1.0f / 2);

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);
            _a1 = BBE_CONJN_4XCF32(_a1);
            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);

            BBE_SAN_4XCF32_IP(_a1, vy1, py1);
            BBE_SAN_4XCF32POS_FP(vy1, py1);
            BBE_SVN_4XCF32_XP(_a0, py0, N / 2 * sizeof(*_x) + 2 * BBE_SIMD_WIDTH);

            py1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_x) - 4 * BBE_SIMD_WIDTH + (uintptr_t)py1);
            px1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_x) - 4 * BBE_SIMD_WIDTH + (uintptr_t)px1);
        } // for (l = 0; l < L; l++) inner loop


    } // for (n = 0...
    {
        /* Last iteration, n = N / (4 * v_sz)-2 */
        px0 = (xb_vecN_4xcf32 *)(_x - v_sz*(n + 1) + N / 4);
        px1 = (xb_vecN_4xcf32 *)(_x + v_sz*n + N / 4 + 1);
        py0 = (xb_vecN_4xcf32 *)(_x - v_sz*(n + 1) + N / 4);
        py1 = (xb_vecN_4xcf32 *)(_x + v_sz*n + N / 4 + 1);

        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(tw2, ptw, -2 * BBE_SIMD_WIDTH);
        /* 10 cycles per pipeline stage in steady state with unroll=1 */
        for (l = 0; l < L; l++)
        {

            vy1 = BBE_ZALIGN();
            vx1 = BBE_LAN_4XCF32_PP(px1);
            /**************** step 1 ****************/
            BBE_LVN_4XCF32_IP(_a0, px0, -2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, vx1, px1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, BBE_CONJN_4XCF32(_a1));

            _b1 = BBE_MULJN_4XCF32(_b1, tw);
            _SCALE_N_4XCF32(_b0, 1.0f / 2);

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);
            _a1 = BBE_CONJN_4XCF32(_a1);
            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);

            BBE_SVN_4XCF32_IP(_a0, py0, -2 * BBE_SIMD_WIDTH);
            BBE_SAN_4XCF32_IP(_a1, vy1, py1);
            /**************** step 2 ****************/
            BBE_LVN_4XCF32_XP(_a0, px0, N / 2 * sizeof(*_x) + 2 * BBE_SIMD_WIDTH);
            BBE_LAN_4XCF32_IP(_a1, vx1, px1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, BBE_CONJN_4XCF32(_a1));

            _b1 = BBE_MULJN_4XCF32(_b1, tw2);
            _SCALE_N_4XCF32(_b0, 1.0f / 2);

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);
            _a1 = BBE_CONJN_4XCF32(_a1);
            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);

            //  Flush not needed, since no write to x[N/2] 
            BBE_SAN_4XCF32_IP(_a1, vy1, py1);
            BBE_SVN_4XCF32T_XP(_a0, py0, N / 2 * sizeof(*_x) + 2 * BBE_SIMD_WIDTH, m);

            py1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_x) - 4 * BBE_SIMD_WIDTH + (uintptr_t)py1);
            px1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_x) - 4 * BBE_SIMD_WIDTH + (uintptr_t)px1);
        } // for (l = 0; l < L; l++) inner loop
    }  /* Last iteration */

    px0 = (xb_vecN_4xcf32 *)(_x);
    py0 = (xb_vecN_4xcf32 *)(_x);
#if 1
    /* 16 cycles per pipeline stage in steady state with unroll=8 */
    for (l = 0; l<L; l++)
    {
        xb_vecN_4xcf32 a0, aN_4;

        BBE_LVN_4XCF32_XP(a0, px0, N / 4 * sizeof(*_x));
        BBE_LVN_4XCF32_XP(aN_4, px0, N / 4 * sizeof(*_x));
        aN_4 = BBE_CONJN_4XCF32(aN_4);
        a0 = BBE_MULJN_4XCF32(BBE_CONSTN_4XCF32(0x3), a0);
        BBE_SVN_4XCF32F_XP(a0, py0, N / 4 * sizeof(*_x), m);
        BBE_SVN_4XCF32F_XP(aN_4, py0, N / 4 * sizeof(*_x), m);
    }
#else
    for (l = 0; l < L; l++)
    {
        complex_float *x = _x + N / 2 * l;
        xtcomplexfloat yN_4; // y[N/4]
        xtcomplexfloat y0 = ((xtcomplexfloat*)x)[0];

        /* real y[0]
        y[0].s.re = 0.5*(a0.s.re + a0.s.im);
        real(y[N / 2])
        y[0].s.im = 0.5*(a0.s.re - a0.s.im); */
        y0 = BBE_CONJCF32(y0);
        y0 = BBE_MULCF32(y0, BBE_SELSN_4XCF32(BBE_CONSTN_4XCF32(0x3), 0)); /* (y0_re-1j*y0_im) * (0.5 + 0.5j) */


        /* y[N / 4] = conj_fl32c(x[N / 4]); */
        yN_4 = ((xtcomplexfloat*)x)[N / 4];
        yN_4 = BBE_CONJCF32(yN_4);

        ((xtcomplexfloat*)x)[0] = y0;
        ((xtcomplexfloat*)x)[N / 4] = yN_4;

    } //for (l = 0; l < L; l++)
#endif
} /* blkrifft_spec_conv_packed_fp */

#endif //#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU
#endif //#ifndef _FFT_FP_COMMON_H_

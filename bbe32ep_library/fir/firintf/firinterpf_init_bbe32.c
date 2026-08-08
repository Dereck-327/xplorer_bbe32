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

#include "firinterpf_common.h"

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

static void _alloc_dx(tFilterLayout* pFltr, int M, int D)
{
    pFltr->coefNum = M*D;
    pFltr->delLength = (M + (BBE_SIMD_WIDTH / 4 - 1))&(~(BBE_SIMD_WIDTH / 4 - 1));
}

static void _alloc_d_2_3_4_mx(tFilterLayout* pFltr, int M, int D)
{
    pFltr->coefNum = M*D;
    pFltr->delLength = (((M + (BBE_SIMD_WIDTH / 4 - 1))&(~(BBE_SIMD_WIDTH / 4 - 1))) + 4);
}

static void _alloc_d_6_12_mx(tFilterLayout* pFltr, int M, int D)
{
    pFltr->coefNum = M*D;
    pFltr->delLength = (((M + (BBE_SIMD_WIDTH / 4 - 1))&(~(BBE_SIMD_WIDTH / 4 - 1))) + 4) * D;
}

static void _alloc_dx_mx(tFilterLayout* pFltr, int M, int D)
{
    pFltr->coefNum = M*D;
    pFltr->delLength = ((M + (BBE_SIMD_WIDTH / 4 - 1))&(~(BBE_SIMD_WIDTH / 4 - 1))) * 2 + 4;
}

static void _init_dm(float32_t* restrict coef, const float32_t* restrict h, int M, int D)
{
    int m, d;
    for (d = 0; d < D; d++)
    {
        for (m = 0; m < M; m++)
        {
            coef[m*D + d] = h[D*(M - 1 - m) + d] * D;
        }
    }
    for (m = M * D; m < ((M + (BBE_SIMD_WIDTH / 2 - 1))&(~(BBE_SIMD_WIDTH / 2 - 1))); m++)
    {
        coef[m] = 0.f;
    }
}

static void _init_md(float32_t* restrict coef, const float32_t* restrict h, int M, int D)
{
    int m, d;
    for (d = 0; d < D; d++)
    {
        for (m = 0; m < M; m++)
        {
            coef[d*M + m] = h[D*(M - 1 - m) + d] * D;
        }
    }
    for (m = M * D; m < ((M + (BBE_SIMD_WIDTH / 2 - 1))&(~(BBE_SIMD_WIDTH / 2 - 1))); m++)
    {
        coef[m] = 0.f;
    }
}

const tFirInitAlloc_int firinterpf_dx         = { _alloc_dx,         _init_dm };
const tFirInitAlloc_int firinterpf_d_2_3_4_mx = { _alloc_d_2_3_4_mx, _init_dm };
const tFirInitAlloc_int firinterpf_d_6_12     = { _alloc_dx,         _init_md };
const tFirInitAlloc_int firinterpf_d_6_12_mx  = { _alloc_d_6_12_mx,  _init_md };
const tFirInitAlloc_int firinterpf_dx_mx      = { _alloc_dx_mx,      _init_md };

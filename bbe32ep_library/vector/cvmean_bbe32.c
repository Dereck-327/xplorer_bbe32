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
  NatureDSP_Baseband library. Vector Mathematics
    Mean of Vector Elements
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Mean of Vector Elements

Description: Compute the mean value over all elements of given vector (real 
or complex).

Representation:
rvmean,cvmean    16-bit signed fixed-point format
rvmeanf,cvmeanf  IEEE-754 Std. single precision floating-point format

Parameters:
Input:
x[N]    Input vector
N       Length of input vector, in real or complex samples
Output:
m[1]    Mean value, real or complex

Restrictions:
x       Aligned on 32-byte boundary
x,m     Must not overlap
N       Must be a multiple of either:
          4 (cvmeanf), or
          8 (cvmean, rvmeanf), or
          16 (rvmean)
-------------------------------------------------------------------------*/

void cvmean (complex_fract16 * m, const complex_fract16 * restrict x, int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    int n;
    xb_vecNx16 x0, z0, z1;
    xb_vecNx40 A;

    if (N <= 0) { m[0].s.re = m[0].s.im = 0; return; }
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % (BBE_SIMD_WIDTH / 2) == 0);

    A = 0;
    z0 = BBE_MOVVINT16(1);
    z1 = BBE_MOVVA16(N);
    for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_MULANX16(A, x0, z0);
    }
#if HAVE_DIV
    {
        xb_vecNx16 y0;
        xb_c40   a;
        a = BBE_RADDNX40C(A);
        A = BBE_MOVNX40_FROMC40(a);

        y0 = BBE_QUONX32(A, z1);
        BBE_SSNX16_I(y0, (int16_t *)m, 0);
        y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_1);
        BBE_SSNX16_I(y0, (int16_t *)m, 2);
    }
#else
    {
        int32_t s0, s1;
        xb_c40   a;
        a = BBE_RADDNX40C(A);
        A = BBE_MOVNX40_FROMC40(a);
        s0 = BBE_MOVAW32(A);
        A = BBE_SELNX40I(A, A, BBE_W_SELI_ROTATE_RIGHT_1);
        s1 = BBE_MOVAW32(A);

        s0 = (s0 + N / 2) / N;
        s1 = (s1 + N / 2) / N;
        s0 = XT_MAX(-32768, XT_MIN(32767, s0));
        s1 = XT_MAX(-32768, XT_MIN(32767, s1));
        m->s.re = (int16_t)s0;
        m->s.im = (int16_t)s1;
    }
#endif
} /* cvmean() */

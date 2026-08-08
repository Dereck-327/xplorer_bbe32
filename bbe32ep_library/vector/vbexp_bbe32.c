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
    Common Exponent
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
Common Block Exponent

Description: These functions compute base-2 exponent adjustment term needed
to normalize data in the input vector. Exact meaning of normalization depends
on the data format (see below).

Representation: 
vbexp,vbexp_fast        16-bit signed fixed-point format
                        For each input value functions count the number of
                        redundant sign bits (as if the value was loaded in
                        a 32-bit register) and return the minimum result over
                        the input data vector.
sbexp                   32-bit signed fixed-point format
                        Count the number of redundant sign bits and return 
                        the result.
vbexpf,sbexpf           IEEE-754 Std. single precision floating-point format
                        For each finite input value x, functions estimate the
                        integer E(x), such that 0.5 <= |x|*2^E(x) < 1 and
                        -128 <= E(x) <= 148. The minimum value of E(x) over
                        input data vector is the result.

Special cases:
   x    |  Result |    Extra Conditions    
--------+---------|---------------------------
0       |    0    |
+/-Inf  | -129    | floating-point functions
NaN     |    0    |
--------|---------|---------------------------
0       |   31    |
-32768  |   16    | fixed-point functions
32767   |   16    |

Parameters:
Input:
x[N]    Input data
N       Length of data vector
Returned Value: exponent adjustment term, or zero if N<=0
Restrictions:
vbexp(), vbexpf():
  No restrictions
vbexp_fast():
  x     Aligned on 32-byte boundary
  N     Multiple of 16
-------------------------------------------------------------------------*/

int vbexp ( const int16_t * restrict x, int N )
{
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    int n;
    xb_vecNx16  x0, x0_max, x0_min, x1, x1_max, x1_min, r0, t0;
    valign      x_align;
    xb_int16    r;
    vsaN        nsa;

    x_align = BBE_LANX16_PP(pX);
    x1_max = x0_max =
        x1_min = x0_min = 0;
    for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH + 1)); n++)
    {
        BBE_LANX16_IP(x0, x_align, pX);
        x0_max = BBE_MAXNX16(x0, x0_max);
        x0_min = BBE_MINNX16(x0, x0_min);
        BBE_LANX16_IP(x1, x_align, pX);
        x1_max = BBE_MAXNX16(x1, x1_max);
        x1_min = BBE_MINNX16(x1, x1_min);
    }
    x0_max = BBE_MAXNX16(x1_max, x0_max);
    x0_min = BBE_MINNX16(x1_min, x0_min);
    nsa = BBE_NSANX16(x0_min);
    r0 = BBE_MOVVVS(nsa);
    nsa = BBE_NSANX16(x0_max);
    t0 = BBE_MOVVVS(nsa);
    r0 = BBE_MINNX16(r0, t0);
    n = (N&(2 * BBE_SIMD_WIDTH - 1)) << 1;
    if (n)
    {
        BBE_LAVNX16_XP(x0, x_align, pX, n);
        nsa = BBE_NSANX16(x0);
        t0 = BBE_MOVVVS(nsa);
        r0 = BBE_MINNX16(r0, t0);
        n -= 2 * BBE_SIMD_WIDTH;
        BBE_LAVNX16_XP(x0, x_align, pX, n);
        nsa = BBE_NSANX16(x0);
        t0 = BBE_MOVVVS(nsa);
        r0 = BBE_MINNX16(r0, t0);
    }
    r = BBE_RMINNX16(r0);
    n = r;
    return n + 16;
} /* vbexp() */

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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
    Polynomial Fitting and Interpolation for Real Data
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"
#include "pfit_common.h"

/*-------------------------------------------------------------------------
Polynomial Fitting and Interpolation

Description: the pfit functions fit (in least squares sense) a degree N 
polynomial to input data sampled at a points grid of length M, and use
that polynomial to interpolate data at arbitrary query points. Namely,
the pfitN_grid() functions compute the Vandermonde matrix for the sample
points grid and perform the Cholesky decomposition of that matrix, the
pfitN_process() functions calculate the least squares solution for the
polynomial coefficients. Finally, the pfitN_eval() functions evaluate
the polynomial at query points.

Please refer to the NatureDSP Baseband Library Reference for full details
on these functions.

Representation:
pfit_gridN,      16-bit fixed-point data. Parameter specifications denote
pfit_processN,   fixed-point format for various data items
pfit_evalN       
pfitf_gridN,     IEEE-754 Std single precision floating-point data
pfitf_processN,  
pfitf_evalN     

Note:
Number of fractional bits specidied for various input/output arguments below apply
for the fixed-point variant

Parameters:
Input:
N                     Degree of polynomial, 1..6
M                     Number of sample points
P                     Number of query points
maxIter               Number of least squares solution enhancement iterations. Right 
                      choice depends on required accuracy, the ad-hoc value is (N+1)/2
x[M]                  Sample points grid, Q15 or floating point
y[M]                  Sampled data values, Q15 or floating point
xi[P]                 Query points, Q15 or floating point
M'=(M+7)&(~7), N'=8   for floating point API
M'=(M+15)&(~15),N'=16 for fixed-point API

Intermediate:
V[M'*8]               Vandermonde matrix, Q15 or floating point
R[N'*8]               Upper triangular Cholesky factor of matrix V, Q11 or floating point
Output:
yi[P]                 Data values interpolated at query points, Q15 or floating point
p[N+1]                Polynomial coefficients, Q8.23 or floating point
Temporary:
pScr                  Scratch memory area. To determine the scratch area size required by
                      a function pfitN_<fun>, use the respective helper function 
                      pfit_<fun>_getScratchSize(M,N)

Restrictions:
M>N                   The number of sample points must exceed the degree of polynomial
x,y,xi,yi,V,R,p,pScr  Must not overlap
V,R,pScr              Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
#if !(HAVE_VSAMATH && HAVE_RECIP)
DISCARD_FUN(void, pfit_process3, (
    void * restrict           pScr,
    int32_t * restrict        p,
    const int16_t * restrict  V,
    const int16_t * restrict  R,
    const int16_t * restrict  y,
    int                       M,
    int                       maxIter))
#else

static void pfit_discr4(int16_t* restrict f,
    const int16_t* restrict b,
    const int16_t* restrict A,
    const int32_t* restrict x,
    int M, int N
    )
{
    xb_vecNx40 wA, wB;
    xb_vecNx16 aa, xx, bb, hi, lo, ff, _128 = BBE_MOVVA16(128);
    const xb_vecNx16 *restrict px = (const xb_vecNx16 *)x;
    const xb_vecNx16 *restrict pb;
    const xb_vecNx16 *restrict pA;
    xb_vecNx16 *restrict pf = (xb_vecNx16 *)f;
    valign vb;
    vsaN _7;
    int m;
    NASSERT_ALIGN(A, BBE_SIMD_WIDTH * 2);
    NASSERT_ALIGN(f, BBE_SIMD_WIDTH * 2);

    pb = (const xb_vecNx16 *)b;
    vb = BBE_LA_PP(pb);
    _7 = BBE_MOVVSA32(7);
    pA = (const xb_vecNx16*)A;
    M = M + (-M & (BBE_SIMD_WIDTH - 1));
    
    BBE_LVNX16_XP(aa, pA, M << 1);
    for (m = 0; m<(M >> LOG2_BBE_SIMD_WIDTH); m++)
    {
        BBE_LANX16_IP(bb, vb, pb);
        wA = BBE_MULRNX16(bb, _128, _7);

        BBE_LPNX16_IP(xx, px, 4); lo = BBE_REPNX16(xx, 0); hi = BBE_REPNX16(xx, 1);
        wB = BBE_MULUSNX16(lo, aa);
        wB = BBE_SRAINX40(wB, 16);
        BBE_MULANX16(wB, hi, aa);
        wA = BBE_SUBNX40(wA, wB);

        BBE_LPNX16_IP(xx, px, 4); lo = BBE_REPNX16(xx, 0); hi = BBE_REPNX16(xx, 1);
        BBE_LVNX16_XP(aa, pA, M << 1);
        wB = BBE_MULUSNX16(lo, aa);
        wB = BBE_SRAINX40(wB, 16);
        BBE_MULANX16(wB, hi, aa);
        wA = BBE_SUBNX40(wA, wB);

        BBE_LPNX16_IP(xx, px, 4); lo = BBE_REPNX16(xx, 0); hi = BBE_REPNX16(xx, 1);
        BBE_LVNX16_XP(aa, pA, M << 1);
        wB = BBE_MULUSNX16(lo, aa);
        wB = BBE_SRAINX40(wB, 16);
        BBE_MULANX16(wB, hi, aa);
        wA = BBE_SUBNX40(wA, wB);

        BBE_LPNX16_XP(xx, px, -3 * 4); lo = BBE_REPNX16(xx, 0); hi = BBE_REPNX16(xx, 1);
        BBE_LVNX16_XP(aa, pA, 2 * BBE_SIMD_WIDTH - (4 - 1)*(M << 1));
        wB = BBE_MULUSNX16(lo, aa);
        wB = BBE_SRAINX40(wB, 16);
        BBE_MULANX16(wB, hi, aa);
        wA = BBE_SUBNX40(wA, wB);

        ff = BBE_PACKVNX40(wA, _7);
        BBE_SVNX16_XP(ff, pf, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(aa, pA, M << 1);
    }
}

void pfit_process3 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter )
{
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);

    pfit_process(pScr, p, V, R, y, M, 3, maxIter, pfit_discr4);
} /* pfit_process3() */

#endif

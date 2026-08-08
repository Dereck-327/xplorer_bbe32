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
#if !(HAVE_VSAMATH && HAVE_NSAENX40 && HAVE_RSQRT)
DISCARD_FUN(void, pfit_grid6, (
    void *              pScr,
    int16_t * restrict  V,
    int16_t * restrict  R,
    const int16_t *     x,
    int                 M))
#else

static void makeV6(int16_t* restrict V, const int16_t* restrict x, int M)
{
    const xb_vecNx16* restrict px = (const xb_vecNx16*)x;
    xb_vecNx16* restrict pV = (xb_vecNx16*)V;
    xb_vecNx16 X, Y, one, zero;
    valign vx;
    vboolN b;
    int m/*, _M*/;
    const int Vstride = M + (-M & (BBE_SIMD_WIDTH - 1));

    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    vx = BBE_LA_PP(px);
    one = BBE_MOVVA16(0x7FFF);
    zero = BBE_ZERONX16();

    for (m = 0; m<(M >> LOG2_BBE_SIMD_WIDTH); m++)
    {
        BBE_LANX16_IP(X, vx, px);
        BBE_SVNX16_X(one, pV, 6 * (Vstride << 1));
        BBE_SVNX16_X(X, pV, 5 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, X, X);
        BBE_SVNX16_X(Y, pV, 4 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 3 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 2 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 1 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_IP(Y, pV, 2 * BBE_SIMD_WIDTH);
    }
    M &= (BBE_SIMD_WIDTH - 1);
    if (M != 0)
    {
        b = BBE_LTRN(M);
        one = BBE_MOVNX16T(one, zero, b);
        BBE_LAVNX16_XP(X, vx, px, (M << 1));
        BBE_SVNX16_X(one, pV, 6 * (Vstride << 1));
        BBE_SVNX16_X(X, pV, 5 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, X, X);
        BBE_SVNX16_X(Y, pV, 4 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 3 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 2 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_X(Y, pV, 1 * (Vstride << 1));
        BBE_MULNX16PACKQ_SAT(Y, Y, X);
        BBE_SVNX16_IP(Y, pV, 2 * BBE_SIMD_WIDTH);
    }
    ////zeroing tail for clean ferret
    //_M = (8 - 7) * (Vstride >> LOG2_BBE_SIMD_WIDTH);
    //pV = (xb_vecNx16*)XT_ADD(12*Vstride, (uintptr_t)pV);
    //for (m = 0; m < _M; m++)
    //{
    //    BBE_SVNX16_IP(zero, pV, 2 * BBE_SIMD_WIDTH);
    //}
}

void pfit_grid6 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M )
{
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(M >= 7);

    makeV6(V, x, M); /* construct Vandermonde matrix */
    pfit_chol((int16_t *)pScr, R, V, 100, M, 7); /* make Cholesky decomposition */
} /* pfit_grid6() */

#endif

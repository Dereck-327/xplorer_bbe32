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

void pfit_eval1 (
  int16_t * restrict       yi, /* (out) Data values interpolated at query points */
  const int16_t * restrict xi, /* (in ) Query points                             */
  const int32_t * restrict p,  /* (in ) Polynomial coefficients                  */
  int P                        /* (in ) Number of query points                   */
)
{
    int n;
    xb_vecNx16 yp, xp, lo, hi;
    xb_vecNx40 Am, Ap;
    valign vx, vy;
    vsaN shr;
    int bytecount;
    int rbytecount;

    const xb_vecNx16 * restrict pxi = (const xb_vecNx16 *)xi;
    xb_vecNx16 * restrict pyi = (xb_vecNx16 *)yi;
    vx = BBE_LA_PP(pxi);
    shr = BBE_MOVVSA32(8);
    rbytecount = bytecount = P << 1;
    rbytecount += (P >> 31);  // add zero !!!! to disable compiler's remapping of bytecount to bytecount
    vy = BBE_ZALIGN();

    __Pragma("loop_count min=1");
    for (n = 0; n<((P + BBE_SIMD_WIDTH - 1) >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
        Am = BBE_MOVWA32(p[0]);
        TAKEHILO3(Am, hi, lo);
        Am = BBE_MULUSNX16(lo, xp);
        Am = BBE_SRAINX40(Am, 16);
        BBE_MULANX16(Am, hi, xp);
        Am = BBE_SLLINX40(Am, 1);
        Ap = BBE_MOVWA32(p[1]);
        Am = BBE_ADDNX40(Am, Ap);
        Am = BBE_RNDSADJNX40(Am, shr);
        yp = BBE_PACKVNX40(Am, shr);
        BBE_SAVNX16_XP(yp, vy, pyi, bytecount);

        rbytecount = XT_ADDX8(-(BBE_SIMD_WIDTH / 4), rbytecount);
        bytecount = XT_ADDX8(-(BBE_SIMD_WIDTH / 4), bytecount);
    }
    BBE_SAPOS_FP(vy, pyi);
} /* pfit_eval1() */

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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
  Polynomial Fitting and Interpolation for Real Data
  IntegrIT, 2006-2016
*/
#ifndef PFITF_COMMON_H__
#define PFITF_COMMON_H__
/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"

/*-----------------------------
    Cholesky forward recursion:
    y=R'\(A'*B)
    Input:
    A[M'*N]  
    B[M]     
    R[8*8]  
    M'       (M+7)&~7
    Output:
    Y[N]     
-----------------------------*/
void pfitf_cholfwd (
            float32_t* restrict y, 
            const float32_t* restrict R, 
            const float32_t* restrict D, 
            const float32_t* restrict A, 
            const float32_t* restrict B, 
            int M, int N);

/*-----------------------------
    Cholesky backward recursion:
    x=R\y
    Input:
    Rt[8*8]  transposed R
    y[N]    
    Output:
    x[N]    
-----------------------------*/
void  pfitf_cholbkw (
            float32_t* restrict x, 
            const float32_t* restrict R, 
            const float32_t* restrict D, 
            const float32_t* restrict y, 
            int N);

/*-----------------------------
    Cholesky 
    Input:
    A[M'*N]
    M'=(M+7)&~7
    Input/Output:
    R[N*8]  

    Scratch:
    Am[M]
    Rm[N]
-----------------------------*/
void pfitf_chol(       float32_t* t,
                       float32_t * restrict R, 
                const float32_t* restrict A, 
                float32_t sigma2,
                int M, int N);

typedef void (*fndiscr)(      float32_t* restrict f,
                  const float32_t* restrict b,
                  const float32_t* restrict A,
                  const float32_t* restrict x,
                  int M,int N);

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
void pfitf_process
(
  void * restrict             pScr,
  float32_t * restrict        p,
  const float32_t * restrict  V,
  const float32_t * restrict  R,
  const float32_t * restrict  y,
  int                         M,
  int                         N,
  int                         maxIter,
  fndiscr discr
);

#endif

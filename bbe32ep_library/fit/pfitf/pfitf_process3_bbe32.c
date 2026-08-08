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
    C code optimized for BBE32EP core with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"
#include "pfitf_common.h"

#if HAVE_VFPU

/*-----------------------------
    calculate discrepancy : 
    f=b-A*x
    Input:
    b[M]     
    A[M'][N] 
    x[N]     
    M'       closer biggest multiple of BBE_SIMD_WIDTH/2
    Output:
    f[M]     
-----------------------------*/
static void pfitf_discr4(      float32_t* restrict f,
                  const float32_t* restrict b,
                  const float32_t* restrict A,
                  const float32_t* restrict x,
                  int M,int N
                  )
#if 0
{
  int m;
  const int Astride = (M+BBE_SIMD_WIDTH/2-1)&~(BBE_SIMD_WIDTH/2-1);
  NASSERT(N==4);
  NASSERT_ALIGN(A, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(f, 2*BBE_SIMD_WIDTH);
  for(m=0; m<M; m++)
  {
    float32_t Acc=b[m];
    Acc -= x[0]*A[0*Astride+m];
    Acc -= x[1]*A[1*Astride+m];
    Acc -= x[2]*A[2*Astride+m];
    Acc -= x[3]*A[3*Astride+m];
    f[m] = Acc;
  }
}
#else
{
  const int Astride = (M+BBE_SIMD_WIDTH/2-1)&~(BBE_SIMD_WIDTH/2-1);
    const xb_vecN_2xf32 * restrict pb=(const xb_vecN_2xf32 *)b;
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)(A+3*Astride);
          xb_vecN_2xf32 * restrict pf=(      xb_vecN_2xf32 *)f;
    xb_vecN_2xf32 x0,x1,x2,x3;

    int m,nbytes;
    valign ab;
    NASSERT(N==4);
    NASSERT_ALIGN(A, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(f, 2*BBE_SIMD_WIDTH);
    nbytes=M*sizeof(float32_t);
    ab=BBE_LAN_2XF32_PP(pb);
    x0=BBE_LSN_2XF32_I((const xtfloat *)x,0);
    x1=BBE_LSN_2XF32_I((const xtfloat *)x,1*sizeof(float32_t));
    x2=BBE_LSN_2XF32_I((const xtfloat *)x,2*sizeof(float32_t));
    x3=BBE_LSN_2XF32_I((const xtfloat *)x,3*sizeof(float32_t));
    x0=BBE_REPN_2XF32(x0,0);
    x1=BBE_REPN_2XF32(x1,0);
    x2=BBE_REPN_2XF32(x2,0);
    x3=BBE_REPN_2XF32(x3,0);
    for(m=0; m<(Astride>>(LOG2_BBE_SIMD_WIDTH-1)); m++)
    {
        xb_vecN_2xf32 Acc,a0,a1,a2,a3;

        BBE_LAVN_2XF32_XP(Acc,ab,pb,nbytes); nbytes-=sizeof(xb_vecN_2xf32);
        BBE_LVN_2XF32_XP(a3,pA,-Astride*sizeof(float32_t));
        BBE_LVN_2XF32_XP(a2,pA,-Astride*sizeof(float32_t));
        BBE_LVN_2XF32_XP(a1,pA,-Astride*sizeof(float32_t));
        BBE_LVN_2XF32_XP(a0,pA,3*Astride*sizeof(float32_t)+sizeof(xb_vecN_2xf32));

        BBE_MULSN_2XF32(Acc,x3,a3);
        BBE_MULSN_2XF32(Acc,x2,a2);
        BBE_MULSN_2XF32(Acc,x1,a1);
        BBE_MULSN_2XF32(Acc,x0,a0);
        BBE_SVN_2XF32_IP(Acc,pf,sizeof(xb_vecN_2xf32));
    }
}
#endif

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

void pfitf_process3
(
  void * restrict             pScr,
  float32_t * restrict        p,
  const float32_t * restrict  V,
  const float32_t * restrict  R,
  const float32_t * restrict  y,
  int                         M,
  int                         maxIter
)
{
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R   ,2*BBE_SIMD_WIDTH);
    NASSERT(M>=3);
    pfitf_process(pScr,p,V,R,y,M,3,maxIter,pfitf_discr4);
}

#else
DISCARD_FUN(void,pfitf_process3,
(
  void * restrict             pScr,
  float32_t * restrict        p,
  const float32_t * restrict  V,
  const float32_t * restrict  R,
  const float32_t * restrict  y,
  int                         M,
  int                         maxIter
));
#endif

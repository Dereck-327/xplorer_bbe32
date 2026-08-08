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

/*-----------------------------------------------
    construct Vandermonde matrix
    Input:
    x[M]       
    Output:
    V[M'][N+1]
    Temporary:
    t[M]

    M'=((M + 7)& ~7)
-----------------------------------------------*/
static void makeVf2(float32_t* restrict t,float32_t* restrict V,const float32_t* restrict x,int M,int N)
#if 0
{
  int m;
  const int Vstride = ((M + 7)& ~7);

  NASSERT(N==2);
  NASSERT_ALIGN(t,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);

  (void)t;
  for (m=0; m<M; m++)
  {
    V[2*Vstride+m]=1.f;
    V[1*Vstride+m]=x[m];
    V[          m]=x[m]*x[m];
  }
  for (;m<Vstride;m++)
  {
    V[2*Vstride+m]=0.f;
    V[1*Vstride+m]=0.f;
    V[          m]=0.f;
  }
}
#elif 0
{
    int nBytes;
    const int Vstride = ((M + 7)& ~7);
    const   xb_vecN_2xf32 *restrict px=(const xb_vecN_2xf32 *)x;
            xb_vecN_2xf32 *restrict pV=(      xb_vecN_2xf32 *)(V+2*Vstride);
    valign ax;
    int m;

    NASSERT(N==2);
    NASSERT_ALIGN(t,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);

    (void)t;

    ax=BBE_LAN_2XF32_PP(px);
    nBytes=M*sizeof(float32_t);
    for (m=0; m<((M + 7)& ~7)>>(LOG2_BBE_SIMD_WIDTH-1); m++)
    {
      xb_vecN_2xf32 X,X2;
      BBE_LAVN_2XF32_XP(X,ax,px,nBytes);
      X2=BBE_MULN_2XF32(X,X);
      BBE_SVN_2XF32_XP(BBE_CONSTN_2XF32(1),pV,-Vstride*sizeof(float32_t));
      nBytes-=sizeof(X);
      BBE_SVN_2XF32_XP(X  ,pV,-Vstride*sizeof(float32_t));
      BBE_SVN_2XF32_XP(X2 ,pV,2*Vstride*sizeof(float32_t)+sizeof(xb_vecN_2xf32));
    }
}
#else
{
    int nBytes;
    const int Vstride = ((M + 7)& ~7);
    const   xb_vecN_2xf32 *restrict px=(const xb_vecN_2xf32 *)x;
            xb_vecN_2xf32 *restrict pV0=(      xb_vecN_2xf32 *)(V+0*Vstride);
            xb_vecN_2xf32 *restrict pV1=(      xb_vecN_2xf32 *)(V+1*Vstride);
            xb_vecN_2xf32 *restrict pV2=(      xb_vecN_2xf32 *)(V+2*Vstride);
    valign ax;
    int m;

    NASSERT(N==2);
    NASSERT_ALIGN(t,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);

    (void)t;

    ax=BBE_LAN_2XF32_PP(px);
    nBytes=M*sizeof(float32_t);
    for (m=0; m<((M + 7)& ~7)>>(LOG2_BBE_SIMD_WIDTH-1); m++)
    {
      xb_vecN_2xf32 X,X2;
      BBE_LAVN_2XF32_XP(X,ax,px,nBytes);
      X2=BBE_MULN_2XF32(X,X);
      BBE_SVN_2XF32_IP(BBE_CONSTN_2XF32(1),pV2,sizeof(xb_vecN_2xf32));
      BBE_SVN_2XF32_IP(X ,pV1,sizeof(xb_vecN_2xf32));
      BBE_SVN_2XF32_IP(X2,pV0,sizeof(xb_vecN_2xf32));
      nBytes-=sizeof(X);
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
void pfitf_grid2
(
  void *                pScr,
  float32_t * restrict  V, 
  float32_t * restrict  R, 
  const float32_t *     x,
  int                   M
)
{
  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
  NASSERT(M>=2+1);

  makeVf2((float32_t *)pScr,V,x,M,2); /* construct Vandermonde matrix */
  pfitf_chol((float32_t *)pScr,R,V,(float32_t)1.e-6,M,2+1); /* make Cholesky decomposition */
}
#else
DISCARD_FUN(void,pfitf_grid2,
(
  void *                pScr,
  float32_t * restrict  V, 
  float32_t * restrict  R, 
  const float32_t *     x,
  int                   M
));
#endif

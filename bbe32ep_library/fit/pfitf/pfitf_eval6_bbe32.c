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

#if HAVE_VFPU
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
void pfitf_eval6
(
  float32_t * restrict        yi,
  const float32_t * restrict  xi,
  const float32_t * restrict  poly,
  int P
)
#if 0
{
    float32_t p0=poly[0],p1=poly[1],p2=poly[2],p3=poly[3],p4=poly[4],p5=poly[5],p6=poly[6];
    int p;
    for (p=0; p<P; p++)
    {
        float32_t xp=xi[p];
        float32_t A;
        A =        p0;
        A = A*xp + p1;
        A = A*xp + p2;
        A = A*xp + p3;
        A = A*xp + p4;
        A = A*xp + p5;
        A = A*xp + p6;
        yi[p]=A;
    }
}
#else
{
    xb_vecN_2xf32 pp,y,x,t,p0,p1,p2,p3,p4,p5,p6;
    const xb_vecN_2xf32 *pX=(const xb_vecN_2xf32*)xi;
          xb_vecN_2xf32 *pY=(      xb_vecN_2xf32*)yi;
    const xb_vecN_2xf32 *pP=(const xb_vecN_2xf32*)poly;
    valign aX,aY,aP;
    int p,nbytes;
    if (P<=0) return;
    aP=BBE_LAN_2XF32_PP(pP);
    aX=BBE_LAN_2XF32_PP(pX);
    aY=BBE_ZALIGN();
    BBE_LAVN_2XF32_XP(pp,aP,pP,7*sizeof(float32_t));
    p0=BBE_REPN_2XF32(pp,0);
    p1=BBE_REPN_2XF32(pp,1);
    p2=BBE_REPN_2XF32(pp,2);
    p3=BBE_REPN_2XF32(pp,3);
    p4=BBE_REPN_2XF32(pp,4);
    p5=BBE_REPN_2XF32(pp,5);
    p6=BBE_REPN_2XF32(pp,6);
    for (p=0; p<(P>>(LOG2_BBE_SIMD_WIDTH-1)); p++)
    {
        BBE_LAN_2XF32_IP(x,aX,pX);
        y=p0;
        t=p1; BBE_MULAN_2XF32(t,y,x); y=t;
        t=p2; BBE_MULAN_2XF32(t,y,x); y=t;
        t=p3; BBE_MULAN_2XF32(t,y,x); y=t;
        t=p4; BBE_MULAN_2XF32(t,y,x); y=t;
        t=p5; BBE_MULAN_2XF32(t,y,x); y=t;
        t=p6; BBE_MULAN_2XF32(t,y,x); y=t;
        BBE_SAN_2XF32_IP(y,aY,pY);
    }
    nbytes=(P*sizeof(float32_t))&(2*BBE_SIMD_WIDTH-1);
    BBE_LAVN_2XF32_XP(x,aX,pX,nbytes);
    y=p0;
    t=p1; BBE_MULAN_2XF32(t,y,x); y=t;
    t=p2; BBE_MULAN_2XF32(t,y,x); y=t;
    t=p3; BBE_MULAN_2XF32(t,y,x); y=t;
    t=p4; BBE_MULAN_2XF32(t,y,x); y=t;
    t=p5; BBE_MULAN_2XF32(t,y,x); y=t;
    t=p6; BBE_MULAN_2XF32(t,y,x); y=t;
    BBE_SAVN_2XF32_XP(y,aY,pY,nbytes);
    BBE_SAN_2XF32POS_FP(aY,pY);
}
#endif
#else
DISCARD_FUN(void,pfitf_eval6,(float32_t * yi,const float32_t * xi,const float32_t * poly,int P))
#endif

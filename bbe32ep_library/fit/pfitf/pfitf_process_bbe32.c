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

#define RSTRIDE 8

/*
  calculate main diagonal 
  Input:
  N        matrix size
  R[8*?]   R matrix
  Output:   
  D[?]     reciprocal of main diagonals
*/
void pfitf_calcDf(float32_t * restrict D, const float32_t* R, int N)
#if 0
{
  int m;
  NASSERT(N>=2 && N<=7);
  /* TBD: optimize */
  for (m=0; m<N; m++)
  {
    D[m]=R[m*RSTRIDE+m];
  }
  for (; m<BBE_SIMD_WIDTH/2; m++)D[m]=1.f;

  for (m=0; m<BBE_SIMD_WIDTH/2; m++)
  {
    D[m]=1.f/D[m];
  }
}
#else
{
    int m;
    const xtfloat* pR=(const xtfloat*)(&R[(N-1)*(RSTRIDE+1)]);
    xb_vecN_2xf32  r,d;
    d=BBE_CONSTN_2XF32(1);

    for (m=0; m<N; m++)
    {
        BBE_LSN_2XF32_XP(r,pR,-(RSTRIDE+1)*(int)sizeof(float32_t));
        d=BBE_SELN_2XF32I(d,r,BBE_SELI_PACK_2);
    }
    d=BBE_RECIPN_2XF32(d);   
    BBE_SVN_2XF32_I(d,(xb_vecN_2xf32*)D,0);
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
)
{
  float32_t *yy; // [N+1]
  float32_t *xx; // [N+1]
  float32_t *D;  // [N+1]
  float32_t *b;  // [M]
  float32_t *Rt;  // [8x8]
  int iter;

  NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(R   ,2*BBE_SIMD_WIDTH);
  NASSERT(N>=1 && N<=6);
  NASSERT(M>=N);

  { /* allocate scratch */
    uintptr_t a=(uintptr_t)pScr;
    size_t xx_sz,yy_sz,D_sz,b_sz,Rt_sz;
    xx_sz=(N+1)*sizeof(float32_t); xx_sz = (xx_sz+ 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    yy_sz=(N+1)*sizeof(float32_t); yy_sz = (yy_sz+ 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    D_sz =8*sizeof(float32_t);     D_sz  = (D_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    b_sz =M*sizeof(float32_t);     b_sz  = (b_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Rt_sz =8*8*sizeof(float32_t);  Rt_sz  = (Rt_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);

    yy=(float32_t*)a; a+=xx_sz;
    xx=(float32_t*)a; a+=yy_sz;
    D =(float32_t*)a; a+=D_sz;
    b =(float32_t*)a; a+=b_sz;
    Rt=(float32_t*)a; a+=Rt_sz;
  }
  N=N+1;

  NASSERT_ALIGN(xx, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(yy, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(b , 2*BBE_SIMD_WIDTH);

  /* transpose R */
#if 0
  {
      int m,n;
      for (m=0; m<8; m++)
      for (n=0; n<8; n++) Rt[8*m+n]=R[m+8*n];
  }
#else
  {
        xb_vecN_2xf32 r0,r1,r2,r3,r4,r5,r6,r7;
        r0=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,0*sizeof(xb_vecN_2xf32));
        r1=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,1*sizeof(xb_vecN_2xf32));
        r2=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,2*sizeof(xb_vecN_2xf32));
        r3=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,3*sizeof(xb_vecN_2xf32));
        r4=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,4*sizeof(xb_vecN_2xf32));
        r5=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,5*sizeof(xb_vecN_2xf32));
        r6=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,6*sizeof(xb_vecN_2xf32));
        r7=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)R,7*sizeof(xb_vecN_2xf32));
        BBE_DSELN_2XF32I(r1, r0, r1, r0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r3, r2, r3, r2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r5, r4, r5, r4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r7, r6, r7, r6, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r2, r0, r2, r0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r3, r1, r3, r1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r6, r4, r6, r4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r7, r5, r7, r5, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r4, r0, r4, r0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r5, r1, r5, r1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r6, r2, r6, r2, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(r7, r3, r7, r3, BBE_DSELI_DEINTERLEAVE_2);
        BBE_SVN_2XF32_I(r0,(xb_vecN_2xf32*)Rt,0*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r1,(xb_vecN_2xf32*)Rt,1*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r2,(xb_vecN_2xf32*)Rt,2*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r3,(xb_vecN_2xf32*)Rt,3*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r4,(xb_vecN_2xf32*)Rt,4*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r5,(xb_vecN_2xf32*)Rt,5*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r6,(xb_vecN_2xf32*)Rt,6*sizeof(xb_vecN_2xf32));
        BBE_SVN_2XF32_I(r7,(xb_vecN_2xf32*)Rt,7*sizeof(xb_vecN_2xf32));
    }
#endif
    pfitf_calcDf(D,R,N);
    pfitf_cholfwd (yy, R, D,V,  y, M, N);
    pfitf_cholbkw (p,  Rt, D,yy, N);

    /* calculate discrepancy : v=y-V*p */
    for (iter=0; iter<maxIter; iter++)
    {
        discr(b,y,V,p,M,N);
        pfitf_cholfwd(yy, R, D, V, b, M, N);
        pfitf_cholbkw(xx, Rt, D, yy, N);
#if 0
        { /* auto-vectorizable */
            int n;
            for (n=0; n<N; n++)
            {
            p[n]+=xx[n];
            }
        }
#else
        {
            valign ap;
            xb_vecN_2xf32 XX,PP;
            xb_vecN_2xf32* pP;
            pP=(xb_vecN_2xf32*)p;
            XX=BBE_LVN_2XF32_I((const xb_vecN_2xf32*)xx,0);
            ap=BBE_LAN_2XF32_PP(pP);
            BBE_LAVN_2XF32_XP(PP,ap,pP,N*sizeof(float32_t));
            PP=BBE_ADDN_2XF32(PP,XX);
            ap=BBE_ZALIGN();
            pP=(xb_vecN_2xf32*)p;
            BBE_SAVN_2XF32_XP(PP,ap,pP,N*sizeof(float32_t));
            BBE_SAN_2XF32POS_FP(ap,pP);
        }
#endif
    }
}

/* Return the scratch area size, in bytes. */
size_t pfitf_process_getScratchSize( int M, int N )
{
    size_t xx_sz,yy_sz,D_sz,b_sz,Rt_sz;
    xx_sz=(N+1)*sizeof(float32_t); xx_sz = (xx_sz+ 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    yy_sz=(N+1)*sizeof(float32_t); yy_sz = (yy_sz+ 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    D_sz =8*sizeof(float32_t);     D_sz  = (D_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    b_sz =M*sizeof(float32_t);     b_sz  = (b_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Rt_sz =8*8*sizeof(float32_t);  Rt_sz  = (Rt_sz + 2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    return xx_sz+yy_sz+D_sz+b_sz+Rt_sz;
}


#else

size_t pfitf_process_getScratchSize( int M, int N )
{
  (void)M,(void)N;
  return 0;
}

#endif

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
 * NatureDSP_Baseband Library API
 * Matrix Decomposition and Inversion Functions

    Compute Determinant from LU decomposition for real matrices (stream ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#if HAVE_VFPU
#if 0
#include "NatureDSP_Math.h"
#include <math.h>
#include <float.h>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define VECLEN (BBE_SIMD_WIDTH/2)
#include <complex.h>
static complex_float mulc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)*crealf(y) - cimagf(x)*cimagf(y);
    z.s.im=crealf(x)*cimagf(y) + cimagf(x)*crealf(y);
    return z.u;
}
#endif

/*-------------------------------------------------------------------------
Determinant For Stream Ordered Matrices

Description: compute determinant of a real/complex matrix from its LU 
decomposition (see lu<size>sf() and clu<size>sf() functions) by multiplying
together diagonal elements of the upper triangular factor U. This operation
is accomplished for a sequence of LU matrices stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Input:
  N           Matrix size
  L           Number of matrices
Input/Output:
  LU[N*N][L]  Packed L and U factors computed by [c]lu<size>sf()
Output:
  D[L]        Determinant values
Restrictions:
  D,LU        Must not overlap and must be aligned on 32-byte boundary 
  N           Must be greater than 1
  L           Must be a multiple of 8 for real-valued functions, or a mutiple
              of 4 for complex-valued functions.
---------------------------------------------------------------------------*/
void cdetnxnsf ( complex_float * restrict D, const complex_float * restrict LU, int N, int L )
#if 0
{
    int l,k;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 );
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0 || N<=1) return;

    for (l=0; l<L; l++)
    {
        complex_float d=LU[l];
        for (k=1; k<N; k++) d=mulc(d,LU[l+L*(k+k*N)]);
        D[l]=d;
    }
}
#elif 1
{
    const xb_vecN_4xcf32 * restrict pLU;
          xb_vecN_4xcf32 * restrict pD ;
    int n,l;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0 ) return;

    pD =(      xb_vecN_4xcf32 *)D;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 d,a;
        pLU=((const xb_vecN_4xcf32 *)LU)+l;
        BBE_LVN_4XCF32_XP(d,pLU,L*(N+1)*sizeof(complex_float));
        for (n=1; n<N; n++)
        {
            BBE_LVN_4XCF32_XP(a,pLU,L*(N+1)*sizeof(complex_float));
            d=BBE_MULN_4XCF32(a,d);
        }
        BBE_SVN_4XCF32_IP(d,pD ,2*BBE_SIMD_WIDTH);
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pLU;
          xb_vecN_4xcf32 * restrict pD ;
    const xb_vecN_4xcf32 * restrict pDr;
    int n,l;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0 ) return;

    pLU=(const xb_vecN_4xcf32 *)LU;
    pD =(      xb_vecN_4xcf32 *)D;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 d;
        BBE_LVN_4XCF32_IP(d,pLU,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(d,pD ,2*BBE_SIMD_WIDTH);
    }
    for (n=1; n<N; n++)
    {
        __Pragma("no_reorder")
        pLU=(const xb_vecN_4xcf32 *)(LU+(n+n*N)*L);
        pD =(      xb_vecN_4xcf32 *)D;
        pDr=(const xb_vecN_4xcf32 *)D;
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
        {
            xb_vecN_4xcf32 d,a;
            BBE_LVN_4XCF32_IP(a,pLU,2*BBE_SIMD_WIDTH);
            BBE_LVN_4XCF32_IP(d,pDr,2*BBE_SIMD_WIDTH);
            d=BBE_MULN_4XCF32(a,d);
            BBE_SVN_4XCF32_IP(d,pD ,2*BBE_SIMD_WIDTH);
        }
    }
}
#endif

#else
DISCARD_FUN(void, cdetnxnsf, ( complex_float * restrict D, const complex_float * restrict LU, int N, int L ))
#endif

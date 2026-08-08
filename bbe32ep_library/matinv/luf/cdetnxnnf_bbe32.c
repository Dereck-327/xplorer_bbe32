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

    Compute Determinant from LU decomposition for complex matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#if HAVE_VFPU
// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    m=XT_MIN(m,LOG2_BBE_SIMD_WIDTH-1);
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}
#if 0
#include <math.h>
#include <float.h>
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
Determinant For Block Ordered Matrices

Description: compute determinant of a real/complex matrix from its LU 
decomposition (see lu<size>nf() and clu<size>nf() functions) by multiplying
together diagonal elements of the upper triangular factor U. This operation
is accomplished for a sequence of LU matrices stored in block order.

Storage size SLU denotes the number of data elements required for each NxN LU 
matrix when stored in block order. If matrix size N*N is less than the SIMD
vector size, then the storage size equals N*N rounded up to the next power of
two, otherwise storage size is N*N rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

Data format: IEEE-754 Std single precision floating-point

Input:
  N           Matrix size
  L           Number of matrices
Input/Output:
  LU[L][SLU]  Packed L and U factors computed by [c]lu<size>nf()
Output:
  D[L]        Determinant values
Restrictions:
  D,LU        Must not overlap and must be aligned on 32-byte boundary 
  N           Must be a positive multiple of 4
---------------------------------------------------------------------------*/
void cdetnxnnf   ( complex_float * restrict D, const complex_float * restrict LU, int N, int L )
#if 0
{
    int l,k,SA;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0 || N<=1) return;

    SA=getSpace((N*N)<<1)>>1;
    for (l=0; l<L; l++,LU+=SA)
    {
        complex_float d=LU[0];
        for (k=1; k<N; k++) d=mulc(d,LU[k+k*N]);
        D[l]=d;
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pDrd;
          xb_vecN_4xcf32 * restrict pDwr;
    const xtcomplexfloat * restrict pLU;
    const xtcomplexfloat * restrict pDsrd;
          xtcomplexfloat * restrict pDswr;
    int l,k,K,SA,L0;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    if (L<=0 ) return;

    SA=getSpace((N*N)<<1)>>1;
    L0=L&~(BBE_SIMD_WIDTH/4-1);
    K=L0>>(LOG2_BBE_SIMD_WIDTH-2);
    if (K)
    {
        xb_vecN_4xcf32 t0,t1,t2,t3,T,X;
        xb_vecN_2xf32 tmp;
        pDwr=(      xb_vecN_4xcf32 *)D; 
        tmp=BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0),BBE_CONSTN_2XF32(1),BBE_SELI_INTERLEAVE_2_LO);
        X=BBE_MOVN_4XCF32_FROMN_2XF32(tmp);
        for (l=0; l<K; l++) 
        {
            BBE_SVN_4XCF32_IP(X,pDwr,2*BBE_SIMD_WIDTH);
        }
        for (k=0; k<N; k++) 
        {
            pDrd=(const xb_vecN_4xcf32 *)D; 
            pDwr=(      xb_vecN_4xcf32 *)D; 
            pLU=(const xtcomplexfloat *)&LU[k+k*N];
            for (l=0; l<K; l++)
            {
                BBE_LSN_4XCF32_XP(t0,pLU,SA*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pLU,SA*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pLU,SA*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pLU,SA*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T =BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
                BBE_LVN_4XCF32_IP(X,pDrd,2*BBE_SIMD_WIDTH);
                X=BBE_MULN_4XCF32(X,T);
                BBE_SVN_4XCF32_IP(X,pDwr,2*BBE_SIMD_WIDTH);
            }
            __Pragma("no_reorder");
        }
    }

    D+=L0; LU+=SA*L0; L-=L0;
    if (L>0)
    {
        xb_vecN_4xcf32 X,T;
        xb_vecN_2xf32 tmp;
        pDswr=(      xtcomplexfloat *)D; 
        tmp=BBE_SELN_2XF32I(BBE_CONSTN_2XF32(0),BBE_CONSTN_2XF32(1),BBE_SELI_INTERLEAVE_2_LO);
        X=BBE_MOVN_4XCF32_FROMN_2XF32(tmp);
        for (l=0; l<L; l++) 
        {
            BBE_SSN_4XCF32_IP(X,pDswr,sizeof(complex_float));
        }
        for (k=0; k<N; k++) 
        {
            pDsrd=(const xtcomplexfloat *)D; 
            pDswr=(      xtcomplexfloat *)D; 
            pLU  =(const xtcomplexfloat *)&LU[k+k*N];
            for (l=0; l<L; l++)
            {
                BBE_LSN_4XCF32_XP(T,pLU,SA*sizeof(complex_float));
                BBE_LSN_4XCF32_IP(X,pDsrd,sizeof(complex_float));
                X=BBE_MULN_4XCF32(X,T);
                BBE_SSN_4XCF32_IP(X,pDswr,sizeof(complex_float));
            }
        }
    }
}
#endif

#else
DISCARD_FUN(void, cdetnxnnf ,  ( complex_float * restrict D, const complex_float * restrict LU, int N, int L ))
#endif

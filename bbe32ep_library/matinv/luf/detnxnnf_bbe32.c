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

    Compute Determinant from LU decomposition for real matrices (block ordered)
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
void detnxnnf   ( float32_t * restrict D, const float32_t * restrict LU, int N, int L )
#if 0
{
    int l,k,SA;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(N>0 && N%4==0);
    if (L<=0 || N<=1) return;

    SA=getSpace(N*N);
    for (l=0; l<L; l++,LU+=SA)
    {
        float32_t d=1.0f;
        for (k=0; k<N; k++) d*=LU[k+k*N];
        D[l]=d;
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pDrd;
          xb_vecN_2xf32 * restrict pDwr;
          const xtfloat * restrict pLU;
    const xtfloat * restrict pDsrd;
          xtfloat * restrict pDswr;
    int l,k,K,SA,L0;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    if (L<=0 ) return;

    SA=getSpace(N*N);
    L0=L&~(BBE_SIMD_WIDTH/2-1);
    K=L0>>(LOG2_BBE_SIMD_WIDTH-1);
    if (K)
    {
        xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T,X;
        pDwr=(      xb_vecN_2xf32 *)D; 
        X=BBE_CONSTN_2XF32(1);
        for (l=0; l<K; l++) 
        {
            BBE_SVN_2XF32_IP(X,pDwr,2*BBE_SIMD_WIDTH);
        }
        for (k=0; k<N; k++) 
        {
            pDrd=(const xb_vecN_2xf32 *)D; 
            pDwr=(      xb_vecN_2xf32 *)D; 
            pLU=(const xtfloat *)&LU[k+k*N];
            for (l=0; l<K; l++)
            {
                BBE_LSN_2XF32_XP(t0,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t1,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t2,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t3,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t4,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t5,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t6,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_XP(t7,pLU,SA*sizeof(float32_t));
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
                t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
                T=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
                BBE_LVN_2XF32_IP(X,pDrd,2*BBE_SIMD_WIDTH);
                X=BBE_MULN_2XF32(X,T);
                BBE_SVN_2XF32_IP(X,pDwr,2*BBE_SIMD_WIDTH);
            }
            __Pragma("no_reorder");
        }
    }

    D+=L0; LU+=SA*L0; L-=L0;
    if (L>0)
    {
        pDswr=(      xtfloat *)D; 
        for (l=0; l<L; l++) 
        {
            BBE_SSN_2XF32_IP(BBE_CONSTN_2XF32(1),pDswr,sizeof(float32_t));
        }
        for (k=0; k<N; k++) 
        {
            pDsrd=(const xtfloat *)D; 
            pDswr=(      xtfloat *)D; 
            pLU=(const xtfloat *)&LU[k+k*N];
            for (l=0; l<L; l++)
            {
                xb_vecN_2xf32 X,T;
                BBE_LSN_2XF32_XP(T,pLU,SA*sizeof(float32_t));
                BBE_LSN_2XF32_IP(X,pDsrd,sizeof(float32_t));
                X=BBE_MULN_2XF32(X,T);
                BBE_SSN_2XF32_IP(X,pDswr,sizeof(float32_t));
            }
        }
    }
}
#endif
#else
DISCARD_FUN(void, detnxnnf,   ( float32_t * restrict D, const float32_t * restrict LU, int N, int L ))
#endif

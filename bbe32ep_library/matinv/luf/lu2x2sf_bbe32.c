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

    LU decomposition for real matrices (stream ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_Math.h"
#include "common.h"
#include <math.h>
#include <float.h>

#if HAVE_VFPU

/*-------------------------------------------------------------------------
LU Decomposition For Stream Ordered Matrices

Description: compute LU decomposition of a square matrix using partial pivoting
with row interchanges: P*A = L*U, where P is the permutation matrix, L is the 
lower triangular factor (with diagnoal elements equal to 1), R is the upper 
triangular factor.

Algorithm is applied in-place to a sequence of real/complex matrices stored
in stream order. For each input matrix A[k], the resutling factor U[k] replaces
the upper triangle of input matrix A[k], and subdiagonal elements of factor L[k]
replace the lower triangle of matrix A[k]. Ones on the main diagonal of resulting
factor L[k] are discarded.

Row ordering used by LU decomposition algorithm for input matrix A[k] is stored
to the vector of permutation indices P[0..N-1][k]: i-th row of the matrix was 
interchanged with row P[i][k], i=0..N-1.

Decomposition result is not defined for a close to singular input matrix.

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [c]lu<size>sf_getScratchSize(N,L)
Input:
  N          Matrix size
  L          Number of matrices
Input/Output:
  A[N*N][L]  Input matrices, packed L and U factors on output
Output:
  P[N][L]    Permutation index vectors
Restrictions:
  pScr,A,P   Must not overlap and must be aligned on 32-byte boundary 
  N          Must be greater than 1
  L          Must be a multiple of 8 for real-valued functions, or a multiple
             of 4 for complex-valued functions
---------------------------------------------------------------------------*/
void lu2x2sf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L )
#if 0
{
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
    {
        P[l+L*0]=0;P[l+L*1]=1;
    }

    for (l=0; l<L; l++)
    {
        float32_t maxVal, s, t;
        int pivot;
        float32_t norm;

        /* find those raw permutation which maximizes A(i,i) after the row update */
        maxVal = fmaxf(fabsf(A[l]),FLT_MIN); pivot = 0;
        s = fabsf(A[l+L*2]);
        if (s>maxVal) pivot = 1;
        /* permute i-th and pivot-th rows */
        if (pivot)
        {
            t = A[l+L*(pivot*2 + 0)]; A[l+L*(pivot*2 + 0)] = A[l+L*0]; A[l+L*0] = t;
            t = A[l+L*(pivot*2 + 1)]; A[l+L*(pivot*2 + 1)] = A[l+L*1]; A[l+L*1] = t;
            { int16_t t; t=P[l+L]; P[l+L] = P[l]; P[l] = t; }
        }
        norm=1.0f/A[l];
        A[l+L*2] *= norm;
        A[l+L*3] -= A[l+L*2] * A[l+L*1];
    }
}
#else
{
    valign aP0,aP1;
    const xb_vecN_2xf32 * restrict pAr=(const xb_vecN_2xf32 *)A;
          xb_vecN_2xf32 * restrict pAw=(      xb_vecN_2xf32 *)A;
          xb_vecNx16    * restrict pP0=(      xb_vecNx16    *)P;
          xb_vecNx16    * restrict pP1=(      xb_vecNx16    *)(P+L);
    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    aP0=aP1=BBE_ZALIGN();
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecN_2xf32 A0,A1,A2,A3,MAXVAL,S,T;
        vboolN_2 bPIVOT;
        xb_vecNx16 PIVOT;
        A1=BBE_LVN_2XF32_X (pAr,L*1*sizeof(float32_t));
        A2=BBE_LVN_2XF32_X (pAr,L*2*sizeof(float32_t));
        A3=BBE_LVN_2XF32_X (pAr,L*3*sizeof(float32_t));
        BBE_LVN_2XF32_IP(A0,pAr,2*BBE_SIMD_WIDTH);
        MAXVAL=BBE_MAXN_2XF32(BBE_ABSN_2XF32(A0),FLT_MIN);
        S=BBE_ABSN_2XF32(A2);
        bPIVOT=BBE_OGTN_2XF32(S,MAXVAL);
        T=A0;  A0=BBE_MOVN_2XF32T(A2,A0,bPIVOT); A2=BBE_MOVN_2XF32T(T,A2,bPIVOT);
        T=A1;  A1=BBE_MOVN_2XF32T(A3,A1,bPIVOT); A3=BBE_MOVN_2XF32T(T,A3,bPIVOT);
        PIVOT=BBE_MOVNX16T(1,0,BBE_JOINB(~bPIVOT,bPIVOT));
        BBE_SAVNX16_XP(PIVOT,aP0,pP0,BBE_SIMD_WIDTH/2*sizeof(int16_t));
        BBE_SAPOS_FP(aP0,pP0);
        BBE_SAVNX16_XP(BBE_SELNX16I(PIVOT,PIVOT,BBE_SELI_ROTATE_RIGHT_8),aP1,pP1,BBE_SIMD_WIDTH/2*sizeof(int16_t));
        BBE_SAPOS_FP(aP1,pP1);
        A2=BBE_MULN_2XF32(A2,BBE_RECIPN_2XF32(A0));
        BBE_MULSN_2XF32(A3,A2,A1);
        BBE_SVN_2XF32_X (A1,pAw,L*1*sizeof(float32_t));
        BBE_SVN_2XF32_X (A2,pAw,L*2*sizeof(float32_t));
        BBE_SVN_2XF32_X (A3,pAw,L*3*sizeof(float32_t));
        BBE_SVN_2XF32_IP(A0,pAw,2*BBE_SIMD_WIDTH);
    }
}
#endif

size_t lu2x2sf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==2);
    return 0;
}
#else
DISCARD_FUN(void, lu2x2sf, ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L ))

size_t lu2x2sf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(N==2);
    return 0;
}
#endif

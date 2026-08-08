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
void lu4x4sf ( 
            void *pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int L )
#if 0
{
    float32_t maxVal, s, t;
    float32_t norm;
    int pivot;
    int l,j;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
        for (j=0; j<4; j++) P[l+L*j]=j;

    for (l=0; l<L; l++)
    {
        /* find those raw permutation which maximizes A(i,i) after the row update */
        maxVal = fmaxf(FLT_MIN,fabsf(A[l+L*0])); pivot=0;
        s=fabsf(A[l+L* 4]); if (s>maxVal) { maxVal = s; pivot = 1; }
        s=fabsf(A[l+L* 8]); if (s>maxVal) { maxVal = s; pivot = 2; }
        s=fabsf(A[l+L*12]); if (s>maxVal) { maxVal = s; pivot = 3; }
        /* permute i-th and pivot-th rows */
        t = A[l+L*(pivot*4 + 0)]; A[l+L*(pivot*4 + 0)] = A[l+L*0]; A[l+L*0] = t;
        t = A[l+L*(pivot*4 + 1)]; A[l+L*(pivot*4 + 1)] = A[l+L*1]; A[l+L*1] = t;
        t = A[l+L*(pivot*4 + 2)]; A[l+L*(pivot*4 + 2)] = A[l+L*2]; A[l+L*2] = t;
        t = A[l+L*(pivot*4 + 3)]; A[l+L*(pivot*4 + 3)] = A[l+L*3]; A[l+L*3] = t;
        { int16_t t; t=P[l+L*pivot]; P[l+L*pivot] = P[l+L*0]; P[l+L*0] = t; }
    }

    for (l=0; l<L; l++)
    {
        norm=1.0f/A[l+L* 0];
        A[l+L* 4] *= norm;
        A[l+L* 5] -= A[l+L* 4] * A[l+L*1];
        A[l+L* 6] -= A[l+L* 4] * A[l+L*2];
        A[l+L* 7] -= A[l+L* 4] * A[l+L*3];

        A[l+L* 8] *= norm;
        A[l+L* 9] -= A[l+L* 8] * A[l+L*1];
        A[l+L*10] -= A[l+L* 8] * A[l+L*2];
        A[l+L*11] -= A[l+L* 8] * A[l+L*3];

        A[l+L*12] *= norm;
        A[l+L*13] -= A[l+L*12] * A[l+L*1];
        A[l+L*14] -= A[l+L*12] * A[l+L*2];
        A[l+L*15] -= A[l+L*12] * A[l+L*3];
    }
    for (l=0; l<L; l++)
    {
        /* find those raw permutation which maximizes A(i,i) after the row update */
        maxVal = fmaxf(FLT_MIN,fabsf(A[l+L*5])); pivot=1;
        s=fabsf(A[l+L* 9]); if (s>maxVal) { maxVal = s; pivot = 2; }
        s=fabsf(A[l+L*13]); if (s>maxVal) { maxVal = s; pivot = 3; }

        /* permute i-th and pivot-th rows */
        /* permute i-th and pivot-th rows */
        t = A[l+L*(pivot*4 + 0)]; A[l+L*(pivot*4 + 0)] = A[l+L* 4]; A[l+L* 4] = t;
        t = A[l+L*(pivot*4 + 1)]; A[l+L*(pivot*4 + 1)] = A[l+L* 5]; A[l+L* 5] = t;
        t = A[l+L*(pivot*4 + 2)]; A[l+L*(pivot*4 + 2)] = A[l+L* 6]; A[l+L* 6] = t;
        t = A[l+L*(pivot*4 + 3)]; A[l+L*(pivot*4 + 3)] = A[l+L* 7]; A[l+L* 7] = t;
        { int16_t t; t=P[l+L*pivot]; P[l+L*pivot] = P[l+L*1]; P[l+L*1] = t; }
    }
    for (l=0; l<L; l++)
    {
        norm=1.0f/A[l+L*5];
        A[l+L* 9] *= norm;
        A[l+L*10] -= A[l+L* 9] * A[l+L*6];
        A[l+L*11] -= A[l+L* 9] * A[l+L*7];
        A[l+L*13] *= norm;
        A[l+L*14] -= A[l+L*13] * A[l+L*6];
        A[l+L*15] -= A[l+L*13] * A[l+L*7];
    }
    for (l=0; l<L; l++)
    {
        /* find those raw permutation which maximizes A(i,i) after the row update */
        maxVal = fmaxf(FLT_MIN,fabs(A[l+L*10])); pivot = 2;
        s = fabsf(A[l+L*14]); if (s>maxVal) { maxVal = s; pivot = 3; }
        /* permute i-th and pivot-th rows */
        if (2 != pivot)
        {
            t = A[l+L*12]; A[l+L*12] = A[l+L* 8]; A[l+L* 8] = t;
            t = A[l+L*13]; A[l+L*13] = A[l+L* 9]; A[l+L* 9] = t;
            t = A[l+L*14]; A[l+L*14] = A[l+L*10]; A[l+L*10] = t;
            t = A[l+L*15]; A[l+L*15] = A[l+L*11]; A[l+L*11] = t;
            { int16_t t; t=P[l+L*3]; P[l+L*3] = P[l+L*2]; P[l+L*2] = t; }
        }

        norm=1.0f/A[l+L*10];
        A[l+L*14] *= norm;
        A[l+L*15] -= A[l+L*14] * A[l+L*11];
    }
}
#else
{
    int l;
    const xb_vecN_2xf32 * restrict pAr;
          xb_vecN_2xf32 * restrict pAw;
    xb_vecNx16 * restrict pP0;
    xb_vecNx16 * restrict pP1;
    xb_vecNx16 * restrict pP2;
    xb_vecNx16 * restrict pP3;
    valign aP0,aP1,aP2,aP3;
    xb_vecN_2xf32 A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,NORM,T,MAXVAL;
    vboolN_2 bPIVOT0,bPIVOT1,bPIVOT2,bPIVOT3,bPIVOT4,bPIVOT5;
    vboolN_2 * restrict pbPivot;
    const vboolN_2 * restrict pbPivot0;
    const vboolN_2 * restrict pbPivot1;
    const vboolN_2 * restrict pbPivot2;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    pbPivot=(vboolN_2 *)pScr;
    pAr=(const xb_vecN_2xf32 * )(A+0*L);
    pAw=(      xb_vecN_2xf32 * )(A+0*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A0 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A4 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A8 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,-12*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        MAXVAL=BBE_MAXN_2XF32(FLT_MIN,BBE_ABSN_2XF32(A0));
        T=BBE_ABSN_2XF32(A4 ); bPIVOT0=BBE_OGTN_2XF32(T,MAXVAL); MAXVAL=BBE_MAXN_2XF32(T,MAXVAL);
        T=BBE_ABSN_2XF32(A8 ); bPIVOT1=BBE_OGTN_2XF32(T,MAXVAL); MAXVAL=BBE_MAXN_2XF32(T,MAXVAL);
        T=BBE_ABSN_2XF32(A12); bPIVOT2=BBE_OGTN_2XF32(T,MAXVAL); 
        bPIVOT0 &= ~bPIVOT1;
        bPIVOT0 &= ~bPIVOT2;
        bPIVOT1 &= ~bPIVOT2;
        BBE_SBN_2_IP(bPIVOT0,pbPivot,sizeof(vboolN_2));
        BBE_SBN_2_IP(bPIVOT1,pbPivot,sizeof(vboolN_2));
        BBE_SBN_2_IP(bPIVOT2,pbPivot,sizeof(vboolN_2));
        T=A0;  A0=BBE_MOVN_2XF32T(A4 ,A0,bPIVOT0); A4 =BBE_MOVN_2XF32T(T,A4 ,bPIVOT0);
        T=A0;  A0=BBE_MOVN_2XF32T(A8 ,A0,bPIVOT1); A8 =BBE_MOVN_2XF32T(T,A8 ,bPIVOT1);
        T=A0;  A0=BBE_MOVN_2XF32T(A12,A0,bPIVOT2); A12=BBE_MOVN_2XF32T(T,A12,bPIVOT2);
        BBE_SVN_2XF32_XP(A0 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A4 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A8 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,-12*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pbPivot0=(vboolN_2 *)pScr;
    pAr=(const xb_vecN_2xf32 * )(A+1*L);
    pAw=(      xb_vecN_2xf32 * )(A+1*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        /* permute i-th and pivot-th rows */
        BBE_LBN_2_IP(bPIVOT0,pbPivot0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT1,pbPivot0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT2,pbPivot0,sizeof(vboolN_2));

        BBE_LVN_2XF32_XP(A1 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A5 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,-11*L*sizeof(float32_t));
        T=A1;  A1=BBE_MOVN_2XF32T(A5 ,A1,bPIVOT0); A5 =BBE_MOVN_2XF32T(T,A5 ,bPIVOT0);
        T=A1;  A1=BBE_MOVN_2XF32T(A9 ,A1,bPIVOT1); A9 =BBE_MOVN_2XF32T(T,A9 ,bPIVOT1);
        T=A1;  A1=BBE_MOVN_2XF32T(A13,A1,bPIVOT2); A13=BBE_MOVN_2XF32T(T,A13,bPIVOT2);
        BBE_SVN_2XF32_XP(A1 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A5 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A9 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,-11*L*sizeof(float32_t));

        BBE_LVN_2XF32_XP(A2 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A6 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,-11*L*sizeof(float32_t));
        T=A2;  A2=BBE_MOVN_2XF32T(A6 ,A2,bPIVOT0); A6 =BBE_MOVN_2XF32T(T,A6 ,bPIVOT0);
        T=A2;  A2=BBE_MOVN_2XF32T(A10,A2,bPIVOT1); A10=BBE_MOVN_2XF32T(T,A10,bPIVOT1);
        T=A2;  A2=BBE_MOVN_2XF32T(A14,A2,bPIVOT2); A14=BBE_MOVN_2XF32T(T,A14,bPIVOT2);
        BBE_SVN_2XF32_XP(A2 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A6 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,-11*L*sizeof(float32_t));

        BBE_LVN_2XF32_XP(A3 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A7 ,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr,4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-14*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        T=A3;  A3=BBE_MOVN_2XF32T(A7 ,A3,bPIVOT0); A7 =BBE_MOVN_2XF32T(T,A7 ,bPIVOT0);
        T=A3;  A3=BBE_MOVN_2XF32T(A11,A3,bPIVOT1); A11=BBE_MOVN_2XF32T(T,A11,bPIVOT1);
        T=A3;  A3=BBE_MOVN_2XF32T(A15,A3,bPIVOT2); A15=BBE_MOVN_2XF32T(T,A15,bPIVOT2);
        BBE_SVN_2XF32_XP(A3 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A7 ,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw,4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-14*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }

    __Pragma("no_reorder")
    pAr=(const xb_vecN_2xf32 * )A;
    pAw=(      xb_vecN_2xf32 * )(A+4*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A0 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A1 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A2 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A3 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A4 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A5 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A6 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A7 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A8 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-15*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        NORM = BBE_RECIPN_2XF32(A0);
        A4=BBE_MULN_2XF32(A4,NORM);
        BBE_MULSN_2XF32(A5,A4,A1);
        BBE_MULSN_2XF32(A6,A4,A2);
        BBE_MULSN_2XF32(A7,A4,A3);

        A8=BBE_MULN_2XF32(A8,NORM);
        BBE_MULSN_2XF32(A9 ,A8,A1);
        BBE_MULSN_2XF32(A10,A8,A2);
        BBE_MULSN_2XF32(A11,A8,A3);

        A12=BBE_MULN_2XF32(A12,NORM);
        BBE_MULSN_2XF32(A13,A12,A1);
        BBE_MULSN_2XF32(A14,A12,A2);
        BBE_MULSN_2XF32(A15,A12,A3);

        BBE_SVN_2XF32_XP(A4 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A5 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A6 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A7 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A8 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A9 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-11*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }
#if 0
    pAr=(const xb_vecN_2xf32 * )(A+4*L);
    pAw=(      xb_vecN_2xf32 * )(A+4*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A4 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A5 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A6 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A7 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A8 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-11*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        MAXVAL=BBE_MAXN_2XF32(FLT_MIN,BBE_ABSN_2XF32(A5));
        T=BBE_ABSN_2XF32(A9 ); bPIVOT3=BBE_OGTN_2XF32(T,MAXVAL); MAXVAL=BBE_MAXN_2XF32(T,MAXVAL);
        T=BBE_ABSN_2XF32(A13); bPIVOT4=BBE_OGTN_2XF32(T,MAXVAL); 
        bPIVOT3 &= ~bPIVOT4;
        BBE_SBN_2_IP(bPIVOT3,pbPivot,sizeof(vboolN_2));
        BBE_SBN_2_IP(bPIVOT4,pbPivot,sizeof(vboolN_2));
        /* permute i-th and pivot-th rows */
        T=A4;  A4=BBE_MOVN_2XF32T(A8 ,A4,bPIVOT3); A8 =BBE_MOVN_2XF32T(T,A8 ,bPIVOT3);
        T=A4;  A4=BBE_MOVN_2XF32T(A12,A4,bPIVOT4); A12=BBE_MOVN_2XF32T(T,A12,bPIVOT4);
        T=A5;  A5=BBE_MOVN_2XF32T(A9 ,A5,bPIVOT3); A9 =BBE_MOVN_2XF32T(T,A9 ,bPIVOT3);
        T=A5;  A5=BBE_MOVN_2XF32T(A13,A5,bPIVOT4); A13=BBE_MOVN_2XF32T(T,A13,bPIVOT4);
        T=A6;  A6=BBE_MOVN_2XF32T(A10,A6,bPIVOT3); A10=BBE_MOVN_2XF32T(T,A10,bPIVOT3);
        T=A6;  A6=BBE_MOVN_2XF32T(A14,A6,bPIVOT4); A14=BBE_MOVN_2XF32T(T,A14,bPIVOT4);
        T=A7;  A7=BBE_MOVN_2XF32T(A11,A7,bPIVOT3); A11=BBE_MOVN_2XF32T(T,A11,bPIVOT3);
        T=A7;  A7=BBE_MOVN_2XF32T(A15,A7,bPIVOT4); A15=BBE_MOVN_2XF32T(T,A15,bPIVOT4);

        BBE_SVN_2XF32_XP(A4 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A5 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A6 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A7 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A8 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A9 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-11*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }
#else
    __Pragma("no_reorder")
    pAr=(const xb_vecN_2xf32 * )(A+5*L);
    pAw=(      xb_vecN_2xf32 * )(A+5*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A5 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,-9*L*sizeof(float32_t));
        MAXVAL=BBE_MAXN_2XF32(FLT_MIN,BBE_ABSN_2XF32(A5));
        T=BBE_ABSN_2XF32(A9 ); bPIVOT3=BBE_OGTN_2XF32(T,MAXVAL); MAXVAL=BBE_MAXN_2XF32(T,MAXVAL);
        T=BBE_ABSN_2XF32(A13); bPIVOT4=BBE_OGTN_2XF32(T,MAXVAL); 
        bPIVOT3 &= ~bPIVOT4;
        BBE_SBN_2_IP(bPIVOT3,pbPivot,sizeof(vboolN_2));
        BBE_SBN_2_IP(bPIVOT4,pbPivot,sizeof(vboolN_2));
        /* permute i-th and pivot-th rows */
        T=A5;  A5=BBE_MOVN_2XF32T(A9 ,A5,bPIVOT3); A9 =BBE_MOVN_2XF32T(T,A9 ,bPIVOT3);
        T=A5;  A5=BBE_MOVN_2XF32T(A13,A5,bPIVOT4); A13=BBE_MOVN_2XF32T(T,A13,bPIVOT4);
        BBE_SVN_2XF32_XP(A5 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A9 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,-9*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A4 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A8 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,-6*L*sizeof(float32_t));
        T=A4;  A4=BBE_MOVN_2XF32T(A8 ,A4,bPIVOT3); A8 =BBE_MOVN_2XF32T(T,A8 ,bPIVOT3);
        T=A4;  A4=BBE_MOVN_2XF32T(A12,A4,bPIVOT4); A12=BBE_MOVN_2XF32T(T,A12,bPIVOT4);
        BBE_SVN_2XF32_XP(A4 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A8 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,-6*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A6 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,-7*L*sizeof(float32_t));
        T=A6;  A6=BBE_MOVN_2XF32T(A10,A6,bPIVOT3); A10=BBE_MOVN_2XF32T(T,A10,bPIVOT3);
        T=A6;  A6=BBE_MOVN_2XF32T(A14,A6,bPIVOT4); A14=BBE_MOVN_2XF32T(T,A14,bPIVOT4);
        BBE_SVN_2XF32_XP(A6 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,-7*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A7 ,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr, 4*L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-10*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        T=A7;  A7=BBE_MOVN_2XF32T(A11,A7,bPIVOT3); A11=BBE_MOVN_2XF32T(T,A11,bPIVOT3);
        T=A7;  A7=BBE_MOVN_2XF32T(A15,A7,bPIVOT4); A15=BBE_MOVN_2XF32T(T,A15,bPIVOT4);
        BBE_SVN_2XF32_XP(A7 ,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-10*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }
#endif
    __Pragma("no_reorder")
    pAr=(const xb_vecN_2xf32 * )(A+5*L);
    pAw=(      xb_vecN_2xf32 * )(A+9*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A5 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A6 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A7 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A8 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-10*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
        NORM = BBE_RECIPN_2XF32(A5);
        A9=BBE_MULN_2XF32(A9,NORM);
        BBE_MULSN_2XF32(A10,A9,A6);
        BBE_MULSN_2XF32(A11,A9,A7);

        A13=BBE_MULN_2XF32(A13,NORM);
        BBE_MULSN_2XF32(A14,A13,A6);
        BBE_MULSN_2XF32(A15,A13,A7);

        BBE_SVN_2XF32_XP(A9 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-6*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pAr=(const xb_vecN_2xf32 * )(A+8*L);
    pAw=(      xb_vecN_2xf32 * )(A+8*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        BBE_LVN_2XF32_XP(A8 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A9 ,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A10,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A11,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A12,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A13,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A14,pAr,L*sizeof(float32_t));
        BBE_LVN_2XF32_XP(A15,pAr,-7*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);

        bPIVOT5=BBE_OGTN_2XF32(BBE_ABSN_2XF32(A14),BBE_MAXN_2XF32(FLT_MIN,BBE_ABSN_2XF32(A10)));
        BBE_SBN_2_IP(bPIVOT5,pbPivot,sizeof(vboolN_2));
        /* permute i-th and pivot-th rows */
        T=A8 ;  A8 =BBE_MOVN_2XF32T(A12,A8 ,bPIVOT5); A12=BBE_MOVN_2XF32T(T,A12,bPIVOT5);
        T=A9 ;  A9 =BBE_MOVN_2XF32T(A13,A9 ,bPIVOT5); A13=BBE_MOVN_2XF32T(T,A13,bPIVOT5);
        T=A10;  A10=BBE_MOVN_2XF32T(A14,A10,bPIVOT5); A14=BBE_MOVN_2XF32T(T,A14,bPIVOT5);
        T=A11;  A11=BBE_MOVN_2XF32T(A15,A11,bPIVOT5); A15=BBE_MOVN_2XF32T(T,A15,bPIVOT5);
        NORM = BBE_RECIPN_2XF32(A10);
        A14=BBE_MULN_2XF32(A14,NORM);
        BBE_MULSN_2XF32(A15,A14,A11);

        BBE_SVN_2XF32_XP(A8 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A9 ,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A10,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A11,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A12,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A13,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A14,pAw,L*sizeof(float32_t));
        BBE_SVN_2XF32_XP(A15,pAw,-7*L*sizeof(float32_t)+2*BBE_SIMD_WIDTH);
    }

    /* forming permutation vector from bPIVOTxx from previous stages */
    __Pragma("no_reorder")
    pbPivot0=(vboolN_2*)pScr;
    pbPivot1=pbPivot0+3*(L/(BBE_SIMD_WIDTH/2));
    pbPivot2=pbPivot0+5*(L/(BBE_SIMD_WIDTH/2));
    pP0=(xb_vecNx16 *)(P+0*L);
    pP1=(xb_vecNx16 *)(P+1*L);
    pP2=(xb_vecNx16 *)(P+2*L);
    pP3=(xb_vecNx16 *)(P+3*L);
    aP0=aP1=aP2=aP3=BBE_ZALIGN();
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecN_2xc16 p0,p1,p2,p3,t;
        xb_vecNx16 x0,x1,x2,x3;
        BBE_LBN_2_IP(bPIVOT0,pbPivot0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT1,pbPivot0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT2,pbPivot0,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT3,pbPivot1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT4,pbPivot1,sizeof(vboolN_2));
        BBE_LBN_2_IP(bPIVOT5,pbPivot2,sizeof(vboolN_2));
        //    pivot0 : P0 <-> P1
        //    pivot1 : P0 <-> P2
        //    pivot2 : P0 <-> P3
        p0=BBE_MOVN_2XC16T(1, 0,bPIVOT0);
        p1=BBE_MOVN_2XC16T(0, 1,bPIVOT0);
        p0=BBE_MOVN_2XC16T(2,p0,bPIVOT1);
        p2=BBE_MOVN_2XC16T(0, 2,bPIVOT1);
        p0=BBE_MOVN_2XC16T(3,p0,bPIVOT2);
        p3=BBE_MOVN_2XC16T(0, 3,bPIVOT2);
        //    pivot3: P1 <-> P2
        //    pivot4: P1 <-> P3
        t=p1; p1=BBE_MOVN_2XC16T(p2,p1,bPIVOT3);
              p2=BBE_MOVN_2XC16T(t ,p2,bPIVOT3);
        t=p1; p1=BBE_MOVN_2XC16T(p3,p1,bPIVOT4);
              p3=BBE_MOVN_2XC16T(t ,p3,bPIVOT4);
        //    pivot5: P2 <-> P3
        t=p2; p2=BBE_MOVN_2XC16T(p3,p2,bPIVOT5);
              p3=BBE_MOVN_2XC16T(t ,p3,bPIVOT5);
        x0=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p0,p0,BBE_SELI_EXTRACT_1_OF_2_OFF_0));
        x1=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p1,p1,BBE_SELI_EXTRACT_1_OF_2_OFF_0));
        x2=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p2,p2,BBE_SELI_EXTRACT_1_OF_2_OFF_0));
        x3=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p3,p3,BBE_SELI_EXTRACT_1_OF_2_OFF_0));
        BBE_SAVNX16_XP(x0,aP0,pP0,BBE_SIMD_WIDTH);
        BBE_SAVNX16_XP(x1,aP1,pP1,BBE_SIMD_WIDTH);
        BBE_SAVNX16_XP(x2,aP2,pP2,BBE_SIMD_WIDTH);
        BBE_SAVNX16_XP(x3,aP3,pP3,BBE_SIMD_WIDTH);
    }
    BBE_SAPOS_FP(aP0,pP0);
    BBE_SAPOS_FP(aP1,pP1);
    BBE_SAPOS_FP(aP2,pP2);
    BBE_SAPOS_FP(aP3,pP3);
}
#endif
size_t lu4x4sf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==4);
    return 6*(L/(BBE_SIMD_WIDTH/2))*sizeof(vboolN_2);
}
#else
DISCARD_FUN(void, lu4x4sf, ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L ))

size_t lu4x4sf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(N==4);
    return 0;
}
#endif

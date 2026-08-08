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

    LU decomposition for complex matrices (stream ordered)
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
#if 0
#include <complex.h>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define VECLEN (BBE_SIMD_WIDTH/4)

static complex_float subc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)-crealf(y);
    z.s.im=cimagf(x)-cimagf(y);
    return z.u;
}

static complex_float mulc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)*crealf(y) - cimagf(x)*cimagf(y);
    z.s.im=crealf(x)*cimagf(y) + cimagf(x)*crealf(y);
    return z.u;
}

static complex_float recipc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    float32_t d;
    d=crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
    d=1.0f/d;
    z.s.re= crealf(x)*d;
    z.s.im=-cimagf(x)*d;
    return z.u;
}
static float32_t sqrc(complex_float x)
{
    return crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
}
#endif
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
#include <stdio.h>
void clu3x3sf ( 
            void *pScr,
            complex_float * restrict A,
            int16_t   * restrict P,
            int L )
#if 0
{
    int l;
    float32_t maxVal, s;
    complex_float norm, t;
    int pivot;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
    {
        P[l+L*0]=0; P[l+L*1]=1; P[l+L*2]=2;
    }
    for (l=0; l<L; l++)
    {
        /* find those raw permutation which maximizes A(i,i) after the row update */
        maxVal = fmaxf(FLT_MIN,sqrc(A[l+L*0])); pivot = 0;
        s = sqrc(A[l+L*3]); if (s>maxVal) { maxVal = s; pivot = 1; }
        s = sqrc(A[l+L*6]); if (s>maxVal) { maxVal = s; pivot = 2; }
        /* permute i-th and pivot-th rows */
        t = A[l+L*(pivot*3 + 0)]; A[l+L*(pivot*3 + 0)] = A[l+L*0]; A[l+L*0] = t;
        t = A[l+L*(pivot*3 + 1)]; A[l+L*(pivot*3 + 1)] = A[l+L*1]; A[l+L*1] = t;
        t = A[l+L*(pivot*3 + 2)]; A[l+L*(pivot*3 + 2)] = A[l+L*2]; A[l+L*2] = t;
        P[l+L*pivot] = 0; P[l+L*0] = pivot;
    }
    for (l=0; l<L; l++)
    {
        norm=recipc(A[l+L*0]);
        A[l+L*3] =mulc(A[l+L*3 ], norm);
        A[l+L*4] =subc(A[l+L*4],mulc( A[l+L*3],A[l+L*1]));
        A[l+L*5] =subc(A[l+L*5],mulc( A[l+L*3],A[l+L*2]));
        A[l+L*6] = mulc(A[l+L*6],norm);
        A[l+L*7] = subc(A[l+L*7],mulc(A[l+L*6],A[l+L*1]));
        A[l+L*8] = subc(A[l+L*8],mulc(A[l+L*6],A[l+L*2]));
    }
    // 2-nd iteration
    for (l=0; l<L; l++)
    {
        pivot = sqrc(A[l+L*7])>fmaxf(FLT_MIN,sqrc(A[l+L*4]));
        /* permute i-th and pivot-th rows */
        if (pivot)
        {
            t = A[l+L*6]; A[l+L*6] = A[l+L*3]; A[l+L*3] = t;
            t = A[l+L*7]; A[l+L*7] = A[l+L*4]; A[l+L*4] = t;
            t = A[l+L*8]; A[l+L*8] = A[l+L*5]; A[l+L*5] = t;
            { int16_t t; t=P[l+L*2]; P[l+L*2] = P[l+L*1]; P[l+L*1] = t; }
        }

        norm=recipc(A[l+L*4]);
        A[l+L*7] = mulc(norm,A[l+L*7]);
        A[l+L*8] = subc(A[l+L*8],mulc(A[l+L*7],A[l+L*5]));
    }
}
#else
{
    vboolN_4 * restrict pbPivot;
    const vboolN_4 * restrict pbPivot0;
    const vboolN_4 * restrict pbPivot1;
    const xb_vecN_4xcf32 * restrict pAr;
          xb_vecN_4xcf32 * restrict pAw;
    xb_vecNx16 * restrict pP0;
    xb_vecNx16 * restrict pP1;
    xb_vecNx16 * restrict pP2;
    valign aP0,aP1,aP2;
    xb_vecN_4xcf32 A0,A1,A2,A3,A4,A5,A6,A7,A8,NORM,T;
    xb_vecN_2xf32 S,MAXVAL;
    vboolN_4 bPIVOT0,bPIVOT1,bPIVOT2;

    int l;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    pbPivot = (vboolN_4 *)pScr;
    if (L<=0) return;
    pAr=(const xb_vecN_4xcf32 * )A;
    pAw=(      xb_vecN_4xcf32 * )A;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        BBE_LVN_4XCF32_XP(A0,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A1,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A2,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A3,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A4,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A5,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A6,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A7,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A8,pAr,-8*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        T=BBE_MULMN_4XCF32(A0,A0,0,0); 
        BBE_MULMASN_4XCF32(T,A0,A0,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        MAXVAL=BBE_MAXN_2XF32(S,FLT_MIN);
        T=BBE_MULMN_4XCF32(A3,A3,0,0); 
        BBE_MULMASN_4XCF32(T,A3,A3,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        bPIVOT0=BBE_MOVN_4_FROMN(BBE_MOVN_FROMN_2(BBE_OGTN_2XF32(S,MAXVAL))); 
        MAXVAL=BBE_MAXN_2XF32(S,MAXVAL);
        T=BBE_MULMN_4XCF32(A6,A6,0,0); 
        BBE_MULMASN_4XCF32(T,A6,A6,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        bPIVOT1=BBE_MOVN_4_FROMN(BBE_MOVN_FROMN_2(BBE_OGTN_2XF32(S,MAXVAL))); 
        bPIVOT0 &= ~bPIVOT1;
        BBE_SBN_4_IP(bPIVOT0,pbPivot,sizeof(vboolN_4));
        BBE_SBN_4_IP(bPIVOT1,pbPivot,sizeof(vboolN_4));
        /* permute i-th and pivot-th rows */
        T=A0;  A0=BBE_MOVN_4XCF32T(A3,A0,bPIVOT0); A3=BBE_MOVN_4XCF32T(T,A3,bPIVOT0);
        T=A0;  A0=BBE_MOVN_4XCF32T(A6,A0,bPIVOT1); A6=BBE_MOVN_4XCF32T(T,A6,bPIVOT1);
        T=A1;  A1=BBE_MOVN_4XCF32T(A4,A1,bPIVOT0); A4=BBE_MOVN_4XCF32T(T,A4,bPIVOT0);
        T=A1;  A1=BBE_MOVN_4XCF32T(A7,A1,bPIVOT1); A7=BBE_MOVN_4XCF32T(T,A7,bPIVOT1);
        T=A2;  A2=BBE_MOVN_4XCF32T(A5,A2,bPIVOT0); A5=BBE_MOVN_4XCF32T(T,A5,bPIVOT0);
        T=A2;  A2=BBE_MOVN_4XCF32T(A8,A2,bPIVOT1); A8=BBE_MOVN_4XCF32T(T,A8,bPIVOT1);

        BBE_SVN_4XCF32_XP(A0,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A1,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A2,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A3,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A4,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A5,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A6,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A7,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A8,pAw,-8*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
    }

    __Pragma("no_reorder")
    pAr=(const xb_vecN_4xcf32 * )A;
    pAw=(      xb_vecN_4xcf32 * )(A+3*L);

    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        BBE_LVN_4XCF32_XP(A0,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A1,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A2,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A3,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A4,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A5,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A6,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A7,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A8,pAr,-8*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        NORM = BBE_RECIPN_4XCF32(A0);
        A3=BBE_MULN_4XCF32(A3,NORM);
        BBE_MULSN_4XCF32(A4,A3,A1);
        BBE_MULSN_4XCF32(A5,A3,A2);
        A6=BBE_MULN_4XCF32(A6,NORM);
        BBE_MULSN_4XCF32(A7,A6,A1);
        BBE_MULSN_4XCF32(A8,A6,A2);

        BBE_SVN_4XCF32_XP(A3,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A4,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A5,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A6,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A7,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A8,pAw,-5*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
    }
    // 2-nd iteration
    __Pragma("no_reorder")
    pAr=(const xb_vecN_4xcf32 * )(A+3*L);
    pAw=(      xb_vecN_4xcf32 * )(A+3*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        BBE_LVN_4XCF32_XP(A3,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A4,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A5,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A6,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A7,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A8,pAr,-5*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        T=BBE_MULMN_4XCF32(A4,A4,0,0); 
        BBE_MULMASN_4XCF32(T,A4,A4,0,15); 
        MAXVAL=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        T=BBE_MULMN_4XCF32(A7,A7,0,0); 
        BBE_MULMASN_4XCF32(T,A7,A7,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        bPIVOT2=BBE_MOVN_4_FROMN(BBE_MOVN_FROMN_2(BBE_OGTN_2XF32(S,MAXVAL))); 
        BBE_SBN_4_IP(bPIVOT2,pbPivot,sizeof(vboolN_4));
        /* permute i-th and pivot-th rows */
        T=A3;  A3=BBE_MOVN_4XCF32T(A6,A3,bPIVOT2); A6=BBE_MOVN_4XCF32T(T,A6,bPIVOT2);
        T=A4;  A4=BBE_MOVN_4XCF32T(A7,A4,bPIVOT2); A7=BBE_MOVN_4XCF32T(T,A7,bPIVOT2);
        T=A5;  A5=BBE_MOVN_4XCF32T(A8,A5,bPIVOT2); A8=BBE_MOVN_4XCF32T(T,A8,bPIVOT2);
        BBE_SVN_4XCF32_XP(A3,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A4,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A5,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A6,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A7,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A8,pAw,-5*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pAr=(const xb_vecN_4xcf32 * )(A+4*L);
    pAw=(      xb_vecN_4xcf32 * )(A+8*L);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        BBE_LVN_4XCF32_XP(A4,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A5,pAr,2*L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A7,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A8,pAr,-4*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        NORM = BBE_RECIPN_4XCF32(A4);
        A7=BBE_MULN_4XCF32(A7,NORM);
        BBE_MULSN_4XCF32(A8,A7,A5);
        BBE_SVN_4XCF32_X (A7,pAw,-L*sizeof(complex_float));
        BBE_SVN_4XCF32_IP(A8,pAw,2*BBE_SIMD_WIDTH);
    }

    /* forming permutation vector from bPIVOTxx from previous stages */
    __Pragma("no_reorder")
    pbPivot0=(vboolN_4*)pScr;
    pbPivot1=pbPivot0+2*(L/(BBE_SIMD_WIDTH/4));
    pP0=(xb_vecNx16 *)(P+0*L);
    pP1=(xb_vecNx16 *)(P+1*L);
    pP2=(xb_vecNx16 *)(P+2*L);
    aP0=aP1=aP2=BBE_ZALIGN();
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        xb_vecN_2xc16 p0,p1,p2,t;
        xb_vecNx16 x0,x1,x2;
        BBE_LBN_4_IP(bPIVOT0,pbPivot0,sizeof(vboolN_4));
        BBE_LBN_4_IP(bPIVOT1,pbPivot0,sizeof(vboolN_4));
        BBE_LBN_4_IP(bPIVOT2,pbPivot1,sizeof(vboolN_4));
        p0=BBE_MOVN_2XC16T(1, 0,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT0)));
        p0=BBE_MOVN_2XC16T(2,p0,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT1)));
        p1=BBE_MOVN_2XC16T(0, 1,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT0)));
        p2=BBE_MOVN_2XC16T(0, 2,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT1)));
        t=p1;
        p1=BBE_MOVN_2XC16T(p2,p1,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT2)));
        p2=BBE_MOVN_2XC16T(t ,p2,BBE_MOVN_2_FROMN(BBE_MOVN_FROMN_4(bPIVOT2)));
        x0=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p0,p0,BBE_SELI_EXTRACT_1_OF_4_OFF_0));
        x1=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p1,p1,BBE_SELI_EXTRACT_1_OF_4_OFF_0));
        x2=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(p2,p2,BBE_SELI_EXTRACT_1_OF_4_OFF_0));
        BBE_SAVNX16_XP(x0,aP0,pP0,BBE_SIMD_WIDTH/2);
        BBE_SAVNX16_XP(x1,aP1,pP1,BBE_SIMD_WIDTH/2);
        BBE_SAVNX16_XP(x2,aP2,pP2,BBE_SIMD_WIDTH/2);
    }
    BBE_SAPOS_FP(aP0,pP0);
    BBE_SAPOS_FP(aP1,pP1);
    BBE_SAPOS_FP(aP2,pP2);
}
#endif

size_t clu3x3sf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(N==3);
    return 3*(L/(BBE_SIMD_WIDTH/4))*sizeof(vboolN_4);
}
#else
DISCARD_FUN(void, clu3x3sf, ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int L ))
size_t clu3x3sf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(N==3);
    return 0;
}
#endif

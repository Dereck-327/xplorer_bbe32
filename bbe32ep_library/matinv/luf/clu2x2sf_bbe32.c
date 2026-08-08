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
void clu2x2sf ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int L )
#if 0
{
    int l,p;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
    {
        P[l+L*0]=0;P[l+L*1]=1;
    }

    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        float32_t s[BBE_SIMD_WIDTH/4],maxVal[BBE_SIMD_WIDTH/4];
        complex_float t[BBE_SIMD_WIDTH/4];
        int pivot[BBE_SIMD_WIDTH/4];
        complex_float norm[BBE_SIMD_WIDTH/4];

        /* find those raw permutation which maximizes A(i,i) after the row update */
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) { maxVal[p] = fmaxf(sqrc(A[l+p]),FLT_MIN); pivot[p] = 0;}
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) { s[p] = sqrc(A[l+p+L*2]);                              }
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) { if (s[p]>maxVal[p]) pivot[p] = 1;                     }
        /* permute i-th and pivot-th rows */
        for (p=0;p<BBE_SIMD_WIDTH/4;p++)if (pivot[p])
        {
            t[p] = A[l+p+L*(pivot[p]*2 + 0)]; A[l+p+L*(pivot[p]*2 + 0)] = A[l+p+L*0]; A[l+p+L*0] = t[p];
            t[p] = A[l+p+L*(pivot[p]*2 + 1)]; A[l+p+L*(pivot[p]*2 + 1)] = A[l+p+L*1]; A[l+p+L*1] = t[p];
            { int16_t t; t=P[l+p+L]; P[l+p+L] = P[l+p]; P[l+p] = t; }
        }
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) norm[p]=recipc(A[l+p]);
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) A[l+p+L*2] = mulc(A[l+p+L*2],norm[p]);
        for (p=0;p<BBE_SIMD_WIDTH/4;p++) A[l+p+L*3] = subc(A[l+p+L*3], mulc( A[l+p+L*2] , A[l+p+L*1]));
    }
}
#else
{
    xb_vecNx16    * restrict pP0=(      xb_vecNx16    *)P;
    xb_vecNx16    * restrict pP1=(      xb_vecNx16    *)(P+L);
    int l;
    valign aP0,aP1;

    const xb_vecN_4xcf32 *restrict pAr;
          xb_vecN_4xcf32 *restrict pAw;
    xb_vecN_4xcf32 A0,A1,A2,A3,NORM,T;
    xb_vecN_2xf32 MAXVAL,S;
    vboolN_4 bPIVOT;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0) return;
    aP0=aP1=BBE_ZALIGN();
    pAr=(const xb_vecN_4xcf32 *)A;
    pAw=(      xb_vecN_4xcf32 *)A;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        xb_vecNx16 PIVOT;
        BBE_LVN_4XCF32_XP(A0,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A1,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A2,pAr,L*sizeof(complex_float));
        BBE_LVN_4XCF32_XP(A3,pAr,-3*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        T=BBE_MULMN_4XCF32(A0,A0,0,0); 
        BBE_MULMASN_4XCF32(T,A0,A0,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        MAXVAL=BBE_MAXN_2XF32(S,FLT_MIN);
        /* find those raw permutation which maximizes A(i,i) after the row update */
        T=BBE_MULMN_4XCF32(A2,A2,0,0); 
        BBE_MULMASN_4XCF32(T,A2,A2,0,15); 
        S=BBE_MOVN_2XF32_FROMN_4XCF32(T);
        /* permute i-th and pivot-th rows */
        bPIVOT=BBE_MOVN_4_FROMN(BBE_MOVN_FROMN_2( BBE_OGTN_2XF32(S,MAXVAL)));
        T=A0;  A0=BBE_MOVN_4XCF32T(A2,A0,bPIVOT); A2=BBE_MOVN_4XCF32T(T,A2,bPIVOT);
        T=A1;  A1=BBE_MOVN_4XCF32T(A3,A1,bPIVOT); A3=BBE_MOVN_4XCF32T(T,A3,bPIVOT);
        PIVOT=BBE_MOVNX16T(1,0,BBE_MOVN_FROMN_4(bPIVOT));
        PIVOT=BBE_SELNX16I(PIVOT,PIVOT,BBE_SELI_EXTRACT_1_OF_4_OFF_0);
        BBE_SAVNX16_XP(PIVOT,aP0,pP0,(BBE_SIMD_WIDTH/4)*sizeof(int16_t));
        PIVOT=BBE_MOVNX16T(0,1,BBE_MOVN_FROMN_4(bPIVOT));
        PIVOT=BBE_SELNX16I(PIVOT,PIVOT,BBE_SELI_EXTRACT_1_OF_4_OFF_0);
        BBE_SAVNX16_XP(PIVOT,aP1,pP1,(BBE_SIMD_WIDTH/4)*sizeof(int16_t));

        NORM=BBE_RECIPN_4XCF32(A0);
        A2=BBE_MULN_4XCF32(A2,NORM);
        BBE_MULSN_4XCF32(A3,A2,A1);
        BBE_SVN_4XCF32_XP(A0,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A1,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A2,pAw,L*sizeof(complex_float));
        BBE_SVN_4XCF32_XP(A3,pAw,-3*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
    }
    BBE_SAPOS_FP(aP0,pP0);
    BBE_SAPOS_FP(aP1,pP1);
}
#endif

size_t clu2x2sf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(N==2);
    return 0;
}
#else
DISCARD_FUN(void, clu2x2sf, ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int L ))

size_t clu2x2sf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    NASSERT(N==2);
    return 0;
}
#endif

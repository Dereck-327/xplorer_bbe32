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
void lunxnsf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int N, int L )
#if 0
{
    int16_t* pivots=(int16_t*)pScr;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>1);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
    {
        for (i=0; i<N; i++) P[l+L*i]=i;
    }

    for (i=0; i<N-1; i++)
    {
        for (l=0; l<L; l++)
        {
            float32_t maxVal, s;
            int pivot;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            maxVal = FLT_MIN; pivot = i;
            for (j = i; j<N; j++)
            {
                s = A[l+L*(j*N + i)];
                s = fabs(s);
                if (s>maxVal) { maxVal = s; pivot = j; }
            }
            pivots[l]=pivot;
        }
            /* permute i-th and pivot-th rows */
        for (l=0; l<L; l++)
        {
            float32_t t;
            int pivot=pivots[l];
            for (k = 0; k<N; k++)
            {
                t = A[l+L*(pivot*N + k)]; A[l+L*(pivot*N + k)] = A[l+L*(i*N + k)]; A[l+L*(i*N + k)] = t;
            }
            { int16_t t; t=P[l+L*pivot]; P[l+L*pivot] = P[l+L*i]; P[l+L*i] = t; }
        }
        if (i==N-2) // last element to be updated
        {
            for (l=0; l<L; l++)
            {
                float32_t norm;
                norm=1.0f/A[l+L*((N-2)*(N+1))];
                A[l+L*((N-1)*N + (N-2))] *= norm;
                A[l+L*((N-1)*N + (N-1))] -= A[l+L*((N-1)*N + (N-2))] * A[l+L*(((N-2)*N + (N-1)))];
            }
        }
        else
        {
            for (l=0; l<L; l++)
            {
                float32_t norm;
                norm=1.0f/A[l+L*(i*N + i)];
                for (j = i + 1; j < N; j++)
                {
                    A[l+L*(j*N + i)] *= norm;
                    for (k = i + 1; k < N; k++)
                    {
                        A[l+L*(j*N + k)] -= A[l+L*(j*N + i)] * A[l+L*(i*N + k)];
                    }
                }
            }
        }
    }
}
#else
{
    xb_vecNx16* restrict pP;
    valign aP;
          xb_vecN_2xf32 * restrict pAw;
    const xb_vecN_2xf32 * restrict pAr;
    const xb_vecN_2xf32 * restrict pAi;
    const xb_vecN_2xf32 * restrict pAj;
          xtfloat *       restrict pA0;
          xtfloat *       restrict pA1;
          xtfloat *       restrict pA2;
          xtfloat *       restrict pA3;
          xtfloat *       restrict pA4;
          xtfloat *       restrict pA5;
          xtfloat *       restrict pA6;
          xtfloat *       restrict pA7;
    int16_t* pivots=(int16_t*)pScr;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    NASSERT(N>1);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    for (l=0; l<L; l++)
    {
        for (i=0; i<N; i++) P[l+L*i]=i;
    }

    for (i=0; i<N-1; i++)
    {
        /* find those raw permutation which maximizes A(i,i) after the row update */
        aP=BBE_ZALIGN();
        pP=(xb_vecNx16*)pivots;
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            xb_vecNx16 x0;
            vboolN_2 b;
            xb_vecN_2xc16 PIVOT,IDX;
            xb_vecN_2xf32 MAXVAL,Aij,S;
            pAj=(const xb_vecN_2xf32*)&A[l+L*(i*N + i)];
            BBE_LVN_2XF32_XP(Aij,pAj,(L*N)*sizeof(float32_t));
            PIVOT=IDX=i;
            MAXVAL=BBE_MAXN_2XF32(FLT_MIN,BBE_ABSN_2XF32(Aij));
            for (j=i+1; j<N; j++)
            {
                IDX=BBE_ADDN_2XC16(IDX,1);
                BBE_LVN_2XF32_XP(Aij,pAj,L*N*sizeof(float32_t));
                S=BBE_ABSN_2XF32(Aij);
                b=BBE_OGTN_2XF32(S,MAXVAL);
                MAXVAL=BBE_MAXN_2XF32(MAXVAL,S);
                PIVOT=BBE_MOVN_2XC16T(IDX,PIVOT,b);
            }
            x0=BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16I(PIVOT,PIVOT,BBE_SELI_EXTRACT_1_OF_2_OFF_0));
            BBE_SAVNX16_XP(x0,aP,pP,BBE_SIMD_WIDTH); BBE_SAPOS_FP(aP,pP);
        }
        __Pragma("no_reorder")
        __Pragma("no_unroll")
        for (l=0; l<L; l++)
        {
            int pivot=pivots[l];
            int16_t t; t=P[l+L*pivot]; P[l+L*pivot] = P[l+L*i]; P[l+L*i] = t;
        }
        pP=(xb_vecNx16*)pivots;
        aP=BBE_LANX16_PP(pP);
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            xb_vecNx16 PIVOTS,al,ah,addr;
            xb_vecNx40 w;

            BBE_LAVNX16_XP(PIVOTS,aP,pP,BBE_SIMD_WIDTH/2*sizeof(int16_t));
            
            w=BBE_SEQNX40();
            BBE_MULANX16(w,PIVOTS,L*N);
            w=BBE_SLLINX40(w,2);
            w=BBE_ADDNX40(w,(int32_t)(A+l));
            al=BBE_PACKLNX40(w);
            ah=BBE_PACKVNX40(w,16);
            addr=BBE_SELNX16I(ah,al,BBE_SELI_INTERLEAVE_1_LO);
            pA0=(xtfloat *)BBE_EXTRNX16C(addr,0); /* pA0=(xtfloat *)&A[l+0+L*(pivots[l+0]*N)]; */
            pA1=(xtfloat *)BBE_EXTRNX16C(addr,1); /* pA1=(xtfloat *)&A[l+1+L*(pivots[l+1]*N)]; */
            pA2=(xtfloat *)BBE_EXTRNX16C(addr,2); /* pA2=(xtfloat *)&A[l+2+L*(pivots[l+2]*N)]; */
            pA3=(xtfloat *)BBE_EXTRNX16C(addr,3); /* pA3=(xtfloat *)&A[l+3+L*(pivots[l+3]*N)]; */
            pA4=(xtfloat *)BBE_EXTRNX16C(addr,4); /* pA4=(xtfloat *)&A[l+4+L*(pivots[l+4]*N)]; */
            pA5=(xtfloat *)BBE_EXTRNX16C(addr,5); /* pA5=(xtfloat *)&A[l+5+L*(pivots[l+5]*N)]; */
            pA6=(xtfloat *)BBE_EXTRNX16C(addr,6); /* pA6=(xtfloat *)&A[l+6+L*(pivots[l+6]*N)]; */
            pA7=(xtfloat *)BBE_EXTRNX16C(addr,7); /* pA7=(xtfloat *)&A[l+7+L*(pivots[l+7]*N)]; */

            pAw=(xb_vecN_2xf32*)&A[l+L*i*N];
            for (k = 0; k<N; k++)
            {
                xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T0,T1;
                T1=BBE_LVN_2XF32_I(pAw,0);
                t0=BBE_LSN_2XF32_I(pA0,0);
                t1=BBE_LSN_2XF32_I(pA1,0);
                t2=BBE_LSN_2XF32_I(pA2,0);
                t3=BBE_LSN_2XF32_I(pA3,0);
                t4=BBE_LSN_2XF32_I(pA4,0);
                t5=BBE_LSN_2XF32_I(pA5,0);
                t6=BBE_LSN_2XF32_I(pA6,0);
                t7=BBE_LSN_2XF32_I(pA7,0);
                t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
                t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
                t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
                t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
                t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
                t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
                T0=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
                BBE_SVN_2XF32_XP(T0,pAw,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(T1,pA0,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_2 ),pA1,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_4 ),pA2,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_6 ),pA3,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_8 ),pA4,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_10),pA5,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_12),pA6,L*sizeof(float32_t));
                BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_14),pA7,L*sizeof(float32_t));
            }
        }
        __Pragma("no_reorder")
        /* update of (current and below) */
        if (i==N-2) // last element to be updated
        {
            pAi=(const xb_vecN_2xf32 *)&A[L*(i*N + i)];
            pAj=(const xb_vecN_2xf32 *)&A[L*(i*N + i+N)];
            pAw=(      xb_vecN_2xf32 *)pAj;
            for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
            {
                xb_vecN_2xf32 Aii,Aij,Aki,Akj,NORM;
                Aki=BBE_LVN_2XF32_X (pAi,L*sizeof(float32_t));
                BBE_LVN_2XF32_IP(Aii,pAi,2*BBE_SIMD_WIDTH);
                NORM=BBE_RECIPN_2XF32(Aii);
                Akj=BBE_LVN_2XF32_X (pAj,L*sizeof(float32_t));
                BBE_LVN_2XF32_IP(Aij,pAj,2*BBE_SIMD_WIDTH);
                Aij=BBE_MULN_2XF32(Aij,NORM);
                BBE_MULSN_2XF32(Akj,Aij,Aki);
                BBE_SVN_2XF32_X (Akj,pAw,L*sizeof(float32_t));
                BBE_SVN_2XF32_IP(Aij,pAw,2*BBE_SIMD_WIDTH);
            }
        }
        else
        {
            for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
            {
                xb_vecN_2xf32 Aii,Aij,Aki,Akj,NORM;
                pAi=(const xb_vecN_2xf32 *)&A[l+L*(i*N + i)];
                BBE_LVN_2XF32_XP(Aii,pAi,L*sizeof(float32_t));
                pAr=pAi;
                NORM=BBE_RECIPN_2XF32(Aii);
                for (j = i + 1; j < N; j++)
                {
                    pAi=pAr;
                    pAj=(const xb_vecN_2xf32 *)&A[l+L*(j*N + i)];
                    pAw=(xb_vecN_2xf32*)pAj;
                    BBE_LVN_2XF32_XP(Aij,pAj,L*sizeof(float32_t));
                    Aij=BBE_MULN_2XF32(Aij,NORM);
                    BBE_SVN_2XF32_XP(Aij,pAw,L*sizeof(float32_t));
                    __Pragma("loop_count min=2")
                    for (k = i + 1; k < N; k++)
                    {
                        BBE_LVN_2XF32_XP(Akj,pAj,L*sizeof(float32_t));
                        BBE_LVN_2XF32_XP(Aki,pAi,L*sizeof(float32_t));
                        BBE_MULSN_2XF32(Akj,Aij,Aki);
                        BBE_SVN_2XF32_XP(Akj,pAw,L*sizeof(float32_t));
                    }
                }
            }
        }
    }
}
#endif

/* Return the scratch area size, in bytes. */
size_t lunxnsf_getScratchSize   ( int N, int L )
{
    NASSERT(N>1);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    (void)N, (void)L;
    return sizeof(int16_t)*XT_MAX(L,0);
}
#else
DISCARD_FUN(void, lunxnsf ,( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int N, int L ))

size_t lunxnsf_getScratchSize   ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N>1);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0 || N<=1) return 0;
    return 0;
}
#endif

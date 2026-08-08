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

    LU decomposition for complex matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "NatureDSP_Math.h"
#include "clunf_common.h"
#include "common.h"
#include <math.h>
#include <float.h>

#if HAVE_VFPU 
#if 0
#include <complex.h>

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
#endif
typedef void (*fnupdate)(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
/*-------------------------------------------------------------------------
LU Decomposition For Block Ordered Matrices

Description: compute LU decomposition of a square matrix using partial pivoting
with row interchanges: P*A = L*U, where P is the permutation matrix, L is the 
lower triangular factor (with diagnoal elements equal to 1), U is the upper 
triangular factor.

Algorithm is applied in-place to a sequence of real/complex matrices stored
in block order. For each input matrix A[k], the resutling factor U[k] replaces
the upper triangle of input matrix A[k], and subdiagonal elements of factor L[k]
replace the lower triangle of matrix A[k]. Ones on the main diagonal of resulting
factor L[k] are discarded.

Row ordering used by LU decomposition algorithm for input matrix A[k] is stored
to the vector of permutation indices P[k][0..N-1]: i-th row of the matrix was 
interchanged with row P[k][i], i=0..N-1.

Decomposition result is not defined for a close to singular input matrix.

Storage sizes SA,SP denote the number of data elements required to store a
matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)
SP = storage_size(N)

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by 
            functions [c]lu<size>nf_getScratchSize(N,L)
Input:
  N         Matrix size
  L         Number of matrices
Input/Output:
  A[L][SA]  Input matrices, packed L and U factors on output
Output:
  P[L][SP]  Permutation index vectors
Restrictions:
  pScr,A,P  Must not overlap and must be aligned on 32-byte boundary 
  N         Must be a positive multiple of 4
---------------------------------------------------------------------------*/
void clu16x16nf ( 
            void * pScr,
            complex_float * restrict A,
            int16_t   * restrict P,
            int L )
#if 0
{
    complex_float* restrict pNorm=(complex_float*)pScr;
    const int N=16;
    const int SA=256, SP=16;
    int l,i,j,k;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    if (L<=0) return;

    /*
      reference Matlab code
      function [P,A]=lu_fullpivot_2(A)
      sz=size(A);
      N=sz(1);
      assert(sz(1)==sz(2),'for square matrices only');
      P=feval(class(A),eye(N,N));
      for i=1:N
          % find those raw permutation which maximizes A(i,i) 
          % after the row update
          s=A(i:N,i);
          [amax,imax]=max(abs(s));
          imax=imax+i-1;
          % permutation
          t=A(imax,:);    A(imax,:)=A(i,:);    A(i,:)=t;
          t=P(imax,:);    P(imax,:)=P(i,:);    P(i,:)=t;
          for j=(i+1):N
              A(j,i)=A(j,i)/A(i,i);
              k=(i+1:N);
              A(j,k) = A(j,k) - A(j,i)*A(i,k);
          end
      end
      % L=tril(A,-1)+eye(N,N);
      % U=triu(A);
    */
    for (l=0; l<L; l++)
    {
        for (i=0; i<N; i++) P[l*SP+i]=i;
    }

    for (i=0; i<N-1; i++)
    {
        for (l=0; l<L; l++)
        {
            float32_t maxVal, s;
            complex_float t;
            int pivot;

            /* find those raw permutation which maximizes A(i,i) after the row update */
            maxVal = FLT_MIN; pivot = 0;
            for (j = i; j<N; j++)
            {
                t = A[l*SA+j*N + i];
                s = fabs(crealf(t)*crealf(t)+cimagf(t)*cimagf(t));
                if (s>maxVal) { maxVal = s; pivot = j; }
            }
            /* permute i-th and pivot-th rows */
            for (k = 0; k<N; k++)
            {
                t = A[l*SA+pivot*N + k]; A[l*SA+pivot*N + k] = A[l*SA+i*N + k]; A[l*SA+i*N + k] = t;
            }
            { int16_t t; t=P[l*SP+pivot]; P[l*SP+pivot] = P[l*SP+i]; P[l*SP+i] = t; }
        }
        pNorm=(complex_float*)pScr;
        for (l=0; l<L; l++)
        {
            pNorm[l]=recipc(A[l*SA+i*N + i]);
        }
        for (l=0; l<L; l++)
        {
            for (j = i + 1; j < N; j++)
            {
                complex_float norm;
                norm=pNorm[l];
                A[l*SA+j*N + i] = mulc(A[l*SA+j*N + i],norm);
                for (k = i + 1; k < N; k++)
                {
                    A[l*SA+j*N + k] = subc(A[l*SA+j*N + k],mulc(A[l*SA+j*N + i],A[l*SA+i*N + k]));
                }
            }
        }
    }
}
#else
{
    static const fnupdate upd[]=
    {
        clunf_update16,
        clunf_update15,
        clunf_update14,
        clunf_update13,
        clunf_update12,
        clunf_update11,
        clunf_update10,
        clunf_update9,
        clunf_update8,
        clunf_update7,
        clunf_update6,
        clunf_update5,
        clunf_update4,
        clunf_update3,
        clunf_update2
    };
    const xtcomplexfloat  * restrict pAii;
            xb_vecN_4xcf32* restrict pNormwr;
            xb_vecN_4xcf32* restrict pAw;
    xb_vecN_4xcf32 t0,t1,t2,t3;
    complex_float* restrict pNorm=(complex_float*)pScr;
    int16_t      * restrict pPivot=(int16_t*)(pNorm+L);
    xb_vecNx16* restrict pP;
    int l,i;
    int K,nbytes;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P,2*BBE_SIMD_WIDTH);
    if (L<=0) return;

    pP=(xb_vecNx16*)P;
    for (l=0; l<L; l++)
    {
        BBE_SVNX16_IP(BBE_SEQNX16(),pP,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    for (i=0; i<16-1; i++)
    {
        pAii=(xtcomplexfloat*)&A[i];
        /* find those raw permutation which maximizes A(i,i) after the row update */
        if (i<8)
        {
            for (l=0; l<L; l++)
            {
                int pivot,dummy;
                xb_vecN_4xcf32 T0,T1,T2,T3,A0,A1,A2,A3;
                xb_vecN_2xf32 TT0,TT1,maxT;
                vboolN_2 b0,b1;
                vboolN b;
                vselN vb;
                (void)dummy;

                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,16*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T0=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,16*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T1=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,16*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T2=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,16*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T3=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);

                A0=BBE_MULMN_4XCF32(T0,T0,0,0); 
                BBE_MULMASN_4XCF32(A0,T0,T0,0,15); 
                A1=BBE_MULMN_4XCF32(T1,T1,0,0); 
                BBE_MULMASN_4XCF32(A1,T1,T1,0,15); 
                A2=BBE_MULMN_4XCF32(T2,T2,0,0); 
                BBE_MULMASN_4XCF32(A2,T2,T2,0,15); 
                A3=BBE_MULMN_4XCF32(T3,T3,0,0); 
                BBE_MULMASN_4XCF32(A3,T3,T3,0,15); 
                T0=BBE_SELN_4XCF32I(A1,A0,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                T1=BBE_SELN_4XCF32I(A3,A2,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                TT0=BBE_MOVN_2XF32_FROMN_4XCF32(T0);
                TT1=BBE_MOVN_2XF32_FROMN_4XCF32(T1);

                BBE_CONSTN_2XF32T(TT0,0,BBE_LTRN_2(i));
                maxT=BBE_RMAXNUMN_2XF32(BBE_MAXNUMN_2XF32(TT0,TT1));
                b0=BBE_OEQN_2XF32(maxT,TT0);
                b1=BBE_OEQN_2XF32(maxT,TT1);
                b=BBE_JOINBN_2(b1,b0);
                BBE_SQZN(vb,dummy,b);
                pivot=BBE_EXTRNX16(BBE_MOVVVS(vb),0);
                pPivot[l]=pivot;
            }
        }
        else // i>=8
        {
            pAii+=16*8;
            for (l=0; l<L; l++)
            {
                int pivot,dummy;
                xb_vecN_4xcf32 T2,T3,A2,A3;
                xb_vecN_2xf32 TT1;
                xtfloat maxT;
                vboolN_2 b0,b1;
                vboolN b;
                vselN vb;
                (void)dummy;

                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,16*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T2=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
                BBE_LSN_4XCF32_XP(t0,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t1,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t2,pAii,16*sizeof(complex_float));
                BBE_LSN_4XCF32_XP(t3,pAii,(16*8+16)*sizeof(complex_float));
                t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
                t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
                T3=BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);

                A2=BBE_MULMN_4XCF32(T2,T2,0,0); 
                BBE_MULMASN_4XCF32(A2,T2,T2,0,15); 
                A3=BBE_MULMN_4XCF32(T3,T3,0,0); 
                BBE_MULMASN_4XCF32(A3,T3,T3,0,15); 
                T3=BBE_SELN_4XCF32I(A3,A2,BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                TT1=BBE_MOVN_2XF32_FROMN_4XCF32(T3);

                BBE_CONSTN_2XF32T(TT1,0,BBE_LTRN_2(i));
                BBE_RBMAXNUMN_2XF32(b1,maxT,TT1);
                b0=BBE_LTRN_2I(0);
                b=BBE_JOINBN_2(b1,b0);
                BBE_SQZN(vb,dummy,b);
                pivot=BBE_EXTRNX16(BBE_MOVVVS(vb),0);
                pPivot[l]=pivot;
            }
        }
        __Pragma("no_reorder")
        pAw=(xb_vecN_4xcf32*)&A[i*16];
        for (l=0; l<L; l++)
        {
            xb_vecN_4xcf32 Ai0,Ai1,Ai2,Ai3;
            xb_vecN_4xcf32 Aj0,Aj1,Aj2,Aj3;
            int pivot;
            pivot=pPivot[l];
            /* permute i-th and pivot-th rows */
            { int16_t t; t=P[l*16+pivot]; P[l*16+pivot] = P[l*16+i]; P[l*16+i] = t; }
            pivot-=i;
            Ai0=BBE_LVN_4XCF32_I(pAw,0);
            Ai1=BBE_LVN_4XCF32_I(pAw,1*2*BBE_SIMD_WIDTH);
            Ai2=BBE_LVN_4XCF32_I(pAw,2*2*BBE_SIMD_WIDTH);
            Ai3=BBE_LVN_4XCF32_I(pAw,3*2*BBE_SIMD_WIDTH);
            Aj0=BBE_LVN_4XCF32_I((pAw+4*pivot),0);
            Aj1=BBE_LVN_4XCF32_I((pAw+4*pivot),1*2*BBE_SIMD_WIDTH);
            Aj2=BBE_LVN_4XCF32_I((pAw+4*pivot),2*2*BBE_SIMD_WIDTH);
            Aj3=BBE_LVN_4XCF32_I((pAw+4*pivot),3*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Ai0,(pAw+4*pivot),0);
            BBE_SVN_4XCF32_I(Ai1,(pAw+4*pivot),1*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Ai2,(pAw+4*pivot),2*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Ai3,(pAw+4*pivot),3*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Aj1,pAw,1*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Aj2,pAw,2*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I(Aj3,pAw,3*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_XP(Aj0,pAw,256*sizeof(complex_float));
        }
        __Pragma("no_reorder")
        // compute normalization factor
        K=L>>(LOG2_BBE_SIMD_WIDTH-2);
        nbytes=(L&(BBE_SIMD_WIDTH/4-1))*sizeof(complex_float);
        pAii=(const xtcomplexfloat*)&A[i*16 + i];
        pNormwr=(xb_vecN_4xcf32*)pScr;
        for (l=0; l<K; l++)
        {
            xb_vecN_4xcf32 T;
            BBE_LSN_4XCF32_XP(t0,pAii,256*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t1,pAii,256*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t2,pAii,256*sizeof(complex_float));
            BBE_LSN_4XCF32_XP(t3,pAii,256*sizeof(complex_float));
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T =BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            T =BBE_RECIPN_4XCF32(T);
            BBE_SVN_4XCF32_IP(T,pNormwr,2*BBE_SIMD_WIDTH);
        }
        if (nbytes)
        {
            int Alast=(int)(uintptr_t)(A+L*256-1);
            xb_vecN_4xcf32 T;
            valign aN=BBE_ZALIGN();
            BBE_LSN_4XCF32_XP(t0,pAii,256*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t1,pAii,256*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t2,pAii,256*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            BBE_LSN_4XCF32_XP(t3,pAii,256*sizeof(complex_float)); pAii=(xtcomplexfloat*)XT_MIN((int)pAii,Alast);
            t0=BBE_SELN_4XCF32I(t1,t0,BBE_SELI_PACK_4);
            t2=BBE_SELN_4XCF32I(t3,t2,BBE_SELI_PACK_4);
            T =BBE_SELN_4XCF32I(t2,t0,BBE_SELI_PACK_8);
            T =BBE_RECIPN_4XCF32(T);
            BBE_SAVN_4XCF32_XP(T,aN,pNormwr,nbytes);
            BBE_SAN_4XCF32POS_FP(aN,pNormwr);
        }
        upd[i](A+i*16,(complex_float*)pScr,i,16,L);
    }
}
#endif
size_t clu16x16nf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==16);
    return L*sizeof(complex_float)+L*sizeof(int16_t); // for 1/Aii
}

#else
DISCARD_FUN(void, clu16x16nf, ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t   * restrict P,
            int L ))
size_t clu16x16nf_getScratchSize ( int N, int L )
{
    (void)N, (void)L;
    NASSERT(N==16);
    return 0;
}

#endif

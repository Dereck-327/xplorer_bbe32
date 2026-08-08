/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
    Cholesky forward recursion, floating point real data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "rcholnf_common.h"

#if (HAVE_VFPU)

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    MxNxP=8x8x1
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void rcomputeAB8x8x1f(float32_t* Z,const float32_t* A,const float32_t* B, int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pA = (const xb_vecN_2xf32 *)(A);
    const xb_vecN_2xf32 * restrict pB = (const xb_vecN_2xf32 *)(B);
          xb_vecN_2xf32 * restrict pZ = (      xb_vecN_2xf32 *)(Z);

    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, A0, B0;

    for (l = 0; l < L; l++)
    {
        xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;

        BBE_LVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
        t0 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_0X4);
        t2 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_1X4);
        BBE_DSELN_2XF32I(t1, t0, t0, t0, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(t3, t2, t2, t2, BBE_DSELI_DEINTERLEAVE_2);
        t4 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_2X4);
        t6 = BBE_SHFLN_2XF32I(B0, BBE_SHFLI_REP_3X4);
        BBE_DSELN_2XF32I(t5, t4, t4, t4, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(t7, t6, t6, t6, BBE_DSELI_DEINTERLEAVE_2);

        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(A0, t0);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        Acc1 = BBE_MULN_2XF32(A0, t1);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        Acc2 = BBE_MULN_2XF32(A0, t2);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        Acc3 = BBE_MULN_2XF32(A0, t3);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc, A0, t4);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc1, A0, t5);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc2, A0, t6);
        BBE_LVN_2XF32_IP(A0, pA, 2 * BBE_SIMD_WIDTH);
        BBE_MULAN_2XF32(Acc3, A0, t7);

        Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
        Acc = BBE_ADDN_2XF32(Acc, Acc1);
        Acc = BBE_ADDN_2XF32(Acc, Acc2);

        BBE_SVN_2XF32_IP(Acc, pZ, 2 * BBE_SIMD_WIDTH);
    }
}

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SA,SR,SD,SB,SY denote the number of data elements required to store
a matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next power
of two, otherwise it is matrix_size rounded up to the next multiple of the SIMD 
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SB = storage_size(M*P)
SY = storage_size(N*P)

Scratch size in bytes is defined by [r]cholfwd<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 M         Matrix dimension (number of rows in matrices A)
 N         Matrix dimension (number of columns and rows in 
           matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 A[L][SA]  Sequence of L complex matrices A
 B[L][SB]  Sequence of original right-side matrices B
 D[L][SD]  Reciprocal of main diagonal
Output:
 y[L][SY]   Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4 
4. M>=N
---------------------------------------------------------------------------*/
void rcholfwd8x8x1nf(
            void *pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
            int L)
{
    float32_t* Z=(float32_t* )pScr;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);

    // compute A'*B
    rcomputeAB8x8x1f(Z,A,B, L);
    rcholnfFwdrec8(y,R,D,Z,8,L,40,8,8,8);
}

size_t rcholfwd8x8x1nf_getScratchSize  (int M,int N, int P,int L) 
{ 
    size_t Z_size;
    NASSERT(M==8 && N==8 && P==1);
    M=XT_MAX(0,M);
    N=XT_MAX(0,N);
    P=XT_MAX(0,P);
    L=XT_MAX(0,L);
    Z_size= (L*N*P*M*sizeof(float32_t));
    return (Z_size);
}

#else
DISCARD_FUN(void, rcholfwd8x8x1nf,(
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L ))

size_t rcholfwd8x8x1nf_getScratchSize(int M,int N, int P,int L)
{
  (void)M; (void)N; (void)P; (void)L;
  return 0;
}

#endif

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
    Cholesky backward recursion, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
    */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cholnf_common.h"

#if (HAVE_VFPU)

static void cholnfTransformR8(float32_t* Rt, const float32_t* R, int N, int L)
{
    //const int SR = N*(N + 1);
    int l;
    //const long long * restrict pR;
          long long * restrict pRw = (      long long *)Rt;
    //const long long * restrict pR2 = (const long long *)R;
    xb_vecN_4x64 vTmp;
    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7, t8;
    const xb_vecN_2xf32 * restrict pR_r = (const xb_vecN_2xf32 *)R;
          //xb_vecN_2xf32 * restrict pR_w = (      xb_vecN_2xf32 *)Rt;

    __Pragma("loop_count min=1");
    for (l = 0; l < L; l++)
    {
        BBE_LVN_2XF32_IP(t0, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t1, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t2, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t3, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t4, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t5, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t6, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t7, pR_r, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(t8, pR_r, 2 * BBE_SIMD_WIDTH);

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t8, BBE_SHFLI_REP_2X4)));//34
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t6, BBE_SHFLI_REP_2X4)));//26
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t8, BBE_SHFLI_REP_1X4)));//33
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t4, BBE_SHFLI_REP_3X4)));//19
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t6, BBE_SHFLI_REP_1X4)));//25
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t8, BBE_SHFLI_REP_0X4)));//32
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t3, BBE_SHFLI_REP_1X4)));//13
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t4, BBE_SHFLI_REP_2X4)));//18
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t6, BBE_SHFLI_REP_0X4)));//24
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t7, BBE_SHFLI_REP_3X4)));//31
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t2, BBE_SHFLI_REP_0X4)));//8
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t3, BBE_SHFLI_REP_0X4)));//12
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t4, BBE_SHFLI_REP_1X4)));//17
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t5, BBE_SHFLI_REP_3X4)));//23
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t7, BBE_SHFLI_REP_2X4)));//30
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t1, BBE_SHFLI_REP_0X4)));//4
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t1, BBE_SHFLI_REP_3X4)));//7
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t2, BBE_SHFLI_REP_3X4)));//11
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t4, BBE_SHFLI_REP_0X4)));//16
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t5, BBE_SHFLI_REP_2X4)));//22
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t7, BBE_SHFLI_REP_1X4)));//29
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t0, BBE_SHFLI_REP_1X4)));//1
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t0, BBE_SHFLI_REP_3X4)));//3
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t1, BBE_SHFLI_REP_2X4)));//6
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));

        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t2, BBE_SHFLI_REP_2X4)));//10
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t3, BBE_SHFLI_REP_3X4)));//15
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t5, BBE_SHFLI_REP_1X4)));//21
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));
        vTmp = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(t7, BBE_SHFLI_REP_0X4)));//28
        BBE_SSN_4X64_IP(vTmp, pRw, 2 * sizeof(float32_t));


    }
}

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SR,SD,SY,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SY = storage_size(N*P)
SX = storage_size(N*P)

Scratch size in bytes is defined by [r]cholbkw<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 N         Matrix dimension (number of columns and rows in matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal
 y[L][SY]  Sequence of intermediate decision matrices y
Output:         
 x[L][SX]  Sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4
---------------------------------------------------------------------------*/
void cholbkw8x1nf(
            void * pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
    int L)
{
          float32_t* restrict x=(      float32_t*)_x;
    const float32_t* restrict R=(const float32_t*)_R;
    const float32_t* restrict D=(const float32_t*)_D;
    const float32_t* restrict y=(const float32_t*)_y;
    float32_t* Rt = (float32_t*)pScr;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);

    cholnfTransformR8(Rt, R, 8, L);
    cholnfBkwnx1(x, Rt,D,y, 8,L);
}

/*
    return Scratch size
    Input:
    N:      matrix size
    P:      ignored
    L       number of matrices
    */
size_t cholbkw8x1nf_getScratchSize(int N, int P, int L) 
{ 
    NASSERT(N==8 && P==1);
    L=XT_MAX(0,L);
    return (L*N*(N - 1)*sizeof(float32_t));
}

#else
DISCARD_FUN(void, cholbkw8x1nf,(
            void * pScr,
            complex_float * restrict _x, 
      const complex_float * restrict _R,
      const complex_float * restrict _D,
      const complex_float * restrict _y, 
      int L ))

size_t cholbkw8x1nf_getScratchSize(int N, int P, int L)
{
  (void)N; (void)P; (void)L;
  return 0;
}

#endif

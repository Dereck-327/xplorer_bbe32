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
    Cholesky decomposition, floating point complex data, block format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cholnf_common.h"

#if (HAVE_VFPU)

#define VECLEN (BBE_SIMD_WIDTH/4)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-2)) m=LOG2_BBE_SIMD_WIDTH-2;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/* compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   sigma2[L]    regularization terms for each matrix
   n            number of column
   Output:
   Z[L][n+1][2] results
*/
static void conv(float32_t* Z,const float32_t* A,const float32_t * restrict sigma2,int n,int M,int N,int L,int SA)
{
#if 0
    int l, k, m;
    float32_t B_re, B_im;
    float32_t amk_re, amk_im, amn_re, amn_im;
    for (l = 0; l < L; l++)
    {
        for (k = 0; k < n + 1; k++)
        {
            B_re = sigma2[l];
            B_im = 0.f;
            for (m = 0; m < M; m++)
            {
                amk_re = A[2 * (m*N + k) + 0]; amk_im = A[2 * (m*N + k) + 1];
                amn_re = A[2 * (m*N + n) + 0]; amn_im = A[2 * (m*N + n) + 1];
                B_re += amk_re*amn_re + amk_im*amn_im;
                B_im += amk_re*amn_im - amk_im*amn_re;
            }
            Z[0] = B_re; Z[1] = B_im;
            Z += 2;
        }
        A += SA;
    }
#endif // 0
    
    int l, k, m;
          xb_vecN_2xf32 * restrict pZ;
    const xb_vecN_2xf32 * restrict pAk;
    const xb_vecN_2xf32 * restrict pA;
    const xtfloat       * restrict pAn;
    const xtfloat       * restrict pS;
    int delta = 0;

    valign vz;
    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, sigma;
    xb_vecN_2xf32 A0, A1;

    pS = (const xtfloat *)sigma2;
    vz = BBE_ZALIGN();
    pZ = (xb_vecN_2xf32 *)(Z);
    pA = (const xb_vecN_2xf32 *)(A);

    for (l = 0; l < L; l++)
    {
        delta = XT_XOR(delta, delta);
        BBE_LSN_2XF32_IP(sigma, pS, 4);
        sigma = BBE_SHFLN_2XF32I(sigma, BBE_SHFLI_REP_0X4);
        for (k = n + 1; k > 0; k -= VECLEN)
        {
            Acc = sigma;
            Acc1 = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            Acc3 = BBE_ZERON_2XF32();
            pAk = (const xb_vecN_2xf32*)XT_ADDX8(delta, (uintptr_t)pA);
            delta = XT_ADDI_N(delta, VECLEN);
            pAn = (const xtfloat *)XT_ADDX8(n, (uintptr_t)pA);

            __Pragma("loop_count min=2");
            for (m = 0; m < (M >> 1); m++)
            {
                A1 = BBE_LSN_2XF32_I(pAn, 4);
                BBE_LSN_2XF32_XP(A0, pAn, 2 * N * sizeof(float32_t));
                A1 = BBE_SELN_2XF32I(A1, A0, BBE_SELI_PACK_2);
                A1 = BBE_SHFLN_2XF32I(A1, BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(A0, pAk, 2 * N * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, A0, A1, 0, 4);
                BBE_MULMASN_2XF32(Acc1, A0, A1, 2, 11);

                A1 = BBE_LSN_2XF32_I(pAn, 4);
                BBE_LSN_2XF32_XP(A0, pAn, 2 * N * sizeof(float32_t));
                A1 = BBE_SELN_2XF32I(A1, A0, BBE_SELI_PACK_2);
                A1 = BBE_SHFLN_2XF32I(A1, BBE_SHFLI_REP_0X4);
                BBE_LVN_2XF32_XP(A0, pAk, 2 * N * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc2, A0, A1, 0, 4);
                BBE_MULMASN_2XF32(Acc3, A0, A1, 2, 11);
            }
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);
            BBE_SAVN_2XF32_XP(Acc, vz, pZ, 2 * k * sizeof(float32_t));
        }
        pA = (const xb_vecN_2xf32 *)XT_ADDX4(SA, (uintptr_t)pA);
    }
    BBE_SAN_2XF32POS_FP(vz, pZ);
}

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex- or real-valued least squares problem: A*X=B, 
where A is an MxN coefficient matrix with M >= N; X is an NxP matrix of
unknowns; and B is an MxP right hand matrix.

The decomposition results in an upper triangular NxN matrix R with real and
positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the (conjugate) transpose of a matrix, and sigma2*I is
an NxN identity matrix multiplied by the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Matrix R is stored in special format: only upper-diagonal elements are 
stored and they are written column by column. So, total number of elements
in one matrix R is the sum of arithmetic progression 1,2...N == ((N+1)*N)/2

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see [r]cholfwdmxnxpsf() and [r]cholbkwnxpsf(), 
respectively.

Storage sizes SA,SR,SD denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next
power of two, otherwise it is matrix_size rounded up to the next multiple of
the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)

Scratch size in bytes is defined by scratch allocation functions

Data format: IEEE-754 Std. single precision floating-point

Input:
 M, N      Dimensional parameters
 L         Number of matrices
 A[L][SA]  Sequence of L complex matrices A
 sigma2[L] regularization term
Output:
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. M and N must be positive multiples of 4
3. Number of columns for input matrices A must not exceed the number
   of rows: N<=M.
---------------------------------------------------------------------------*/
void cholmxnnf(
            void *pScr,
            complex_float * restrict _R, 
            complex_float * restrict _D,
      const complex_float * restrict _A, 
      const float32_t     * restrict sigma2,
            int M,int N,
            int L)
{
    float32_t *         restrict R=(float32_t *      )_R;
    float32_t *         restrict D=(float32_t *      )_D;
    const float32_t * restrict   A=(const float32_t *)_A;
    static const tCholnfIteration it[]=
    {
        {cholnfConv0_3  ,cholnfFwdrec0 ,cholnfDiagUpd4 },
        {cholnfConv0_3  ,cholnfFwdrec4 ,cholnfDiagUpd4 },
        {cholnfConv0_3  ,cholnfFwdrec4 ,cholnfDiagUpd4 },
        {cholnfConv0_3  ,cholnfFwdrec4 ,cholnfDiagUpd4 },
        {cholnfConv4_7,  cholnfFwdrec4 ,cholnfDiagUpd8 },
        {cholnfConv4_7,  cholnfFwdrec8 ,cholnfDiagUpd8 },
        {cholnfConv4_7,  cholnfFwdrec8 ,cholnfDiagUpd8 },
        {cholnfConv4_7,  cholnfFwdrec8 ,cholnfDiagUpd8 },
        {cholnfConv8_11 ,cholnfFwdrec8 ,cholnfDiagUpd12},
        {cholnfConv8_11 ,cholnfFwdrec12,cholnfDiagUpd12},
        {cholnfConv8_11 ,cholnfFwdrec12,cholnfDiagUpd12},
        {cholnfConv8_11 ,cholnfFwdrec12,cholnfDiagUpd12},
        {cholnfConv12_15,cholnfFwdrec12,cholnfDiagUpd16},
        {cholnfConv12_15,cholnfFwdrec16,cholnfDiagUpd16},
        {cholnfConv12_15,cholnfFwdrec16,cholnfDiagUpd16},
        {cholnfConv12_15,cholnfFwdrec16,cholnfDiagUpd16}
    };
    float32_t* Z=(float32_t*)pScr; // L columns of A'*A
    float32_t* y; // pointer to the new column in R
    int n;
    int SA=2*getSpace(M*N);
    int SR=2*getSpace((N*(N+1))>>1);
    int SD=2*getSpace(N);

    if (M<=0 || N<=0 || L<=0) return;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    y=R;

    for (n=0; n<XT_MIN(N,16); n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        it[n].conv(Z,A,sigma2,n,M,N,L,SA);
        // make forward recursion to update n new column elements
        it[n].fwdrec(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        it[n].diagUpd(y,D+2*n,Z,n,L,SR,SD);
    }

    for (; n<N; n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        conv(Z,A,sigma2,n,M,N,L,SA);
        // make forward recursion to update n new column elements
        cholnfFwdrec(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        cholnfDiagUpd(y,D+2*n,Z,n,L,SR,SD);
    }
} /* cholmxnnf() */

#define CHOLN_SCRATCH(M,N,L) (2*(N)*(L)*sizeof(float32_t))

size_t cholmxnnf_getScratchSize  (int M,int N, int L)
{
    M=XT_MAX(M,0);
    N=XT_MAX(N,0);
    L=XT_MAX(L,0);
    return CHOLN_SCRATCH(M,N,L);
}

#else
DISCARD_FUN(void, cholmxnnf,(
            void *pScr,
            complex_float * restrict _R, 
            complex_float * restrict _D,
      const complex_float * restrict _A, 
      const float32_t     * restrict sigma2,
            int M,int N,
            int L))

size_t cholmxnnf_getScratchSize  (int M,int N, int L)
{
  (void)M; (void)N; (void)L;
  return 0;
}

#endif

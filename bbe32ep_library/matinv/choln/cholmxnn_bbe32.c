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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Cholesky decomposition for block ordered matrices:
    Apply the Cholesky decomposition to the matrix of normal equations system
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "choln_common.h"

#if !HAVE_CHOLN
DISCARD_FUN(void,cholmxnn,( void    * pScr,
                              complex_fract16 * restrict _R, 
                              complex_fract16 * restrict _D,
                        const complex_fract16 * restrict _A, 
                            int M,
                            int N,
                            int L))
#else

/* compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   n            number of column
   Output:
   Z[L][n+1][2] results
*/
static void conv(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA)
#if 0
{
    int l,k,m;
    int64_t B_re,B_im;
    int16_t amk_re,amk_im,amn_re,amn_im;
    for (l=0; l<L; l++)
    {
        for (k=0; k<n+1; k++)
        {
            B_re=B_im=0;
            for (m=0; m<M; m++)
            {
                amk_re=A[2*(m*N+k)+0];amk_im=A[2*(m*N+k)+1];
                amn_re=A[2*(m*N+n)+0];amn_im=A[2*(m*N+n)+1];
                B_re+=amk_re*amn_re+amk_im*amn_im;
                B_im+=amk_re*amn_im-amk_im*amn_re;
            }
            Z[0]=(int32_t)B_re;Z[1]=(int32_t)B_im;
            Z+=2;
        }
        A+=SA;
    }
}
#else
{
    #define VECLEN ((int)(sizeof(xb_vecNx16)/(2*sizeof(int16_t))))
    int l, k, m;
    int delta= 0;
    const xb_vecNx16 *restrict pA;
    const xb_vecNx16 *restrict pRepA;
    const xb_vecNx16 *restrict pA_tmp;
    const xb_vecNx16 *restrict pA2 = (const xb_vecNx16*)A;
          xb_vecNx16 *restrict pZ = (xb_vecNx16*)Z;

    valign A_align, Z_align;
    xb_vecNx16 A0, RepA;
    xb_vecNx16 Zl, Zh;
    xb_vecNx40 acc0;
    Z_align= BBE_ZALIGN();
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        delta= XT_XOR(delta, delta);
        __Pragma("loop_count min=1")
        for (k= n+1; k > 0; k-= VECLEN)
        {
            pA= (const xb_vecNx16*)XT_ADDX4(delta,(uintptr_t)pA2);
            delta= XT_ADDI_N(delta, VECLEN);
            pRepA= (const xb_vecNx16*)XT_ADDX4(n,(uintptr_t)pA2);
            pA_tmp= pA;
            A_align= BBE_LA_PP(pA_tmp);
            BBE_LAVNX16_XP(A0,A_align,pA_tmp,k*4);
            BBE_LPNX16_XP(RepA,pRepA,N*4);
            RepA= BBE_REPNX16C(RepA,0);
            pA= (const xb_vecNx16*)XT_ADDX4(N,(uintptr_t)pA);
            // mul
            acc0 = BBE_MULNX16J(RepA,A0);
            __Pragma("loop_count min=3")
            for (m=1; m < M; m++)
            {
                pA_tmp= pA;
                A_align= BBE_LA_PP(pA_tmp);
                BBE_LAVNX16_XP(A0,A_align,pA_tmp,k*4);
                BBE_LPNX16_XP(RepA,pRepA,N*4);
                RepA= BBE_REPNX16C(RepA,0);
                pA= (const xb_vecNx16*)XT_ADDX4(N,(uintptr_t)pA);
                BBE_MULANX16J(acc0,RepA,A0);
            }
            // store res
            Zl=BBE_MOVSVWL(acc0);
            Zh=BBE_MOVSVWH(acc0);
            BBE_SAVNX16_XP(Zl,Z_align,pZ,k*2*4);
            BBE_SAVNX16_XP(Zh,Z_align,pZ,k*2*4 - 2*BBE_SIMD_WIDTH);
        }
        pA2= (const xb_vecNx16 *)XT_ADDX2(SA,(uintptr_t)pA2);
    }
    BBE_SAPOS_FP(Z_align,pZ);
    #undef VECLEN
}
#endif

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A 
where adj(...) denotes the conjugate transpose of a matrix

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. 
Matrix R is stored in special format: only upper-diagonal elements are 
stored and they are written column by column. So, total number of elements
in one matrix R is the sum of arithmetic progression 1,2...N == ((N+1)*N)/2

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxpn() and cholbkwnxpn(), respectively.

Matrix sizes SA,SR,SD are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements. 
SA=size(M*N)
SR=size(((N+1)*N)/2)
SD=size(N)
Scratch size in bytes is defined by scratch allocation functions

Input:
 M, N         Dimensional parameters
 L            Number of matrices
 A[L][SA]     Sequence of L complex matrices A
Output:
 R[L][SR]     Sequence of L upper triangular complex matrices R
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Number of matrices L must be positive
3. M and N must be positive multiples of 4
4. Number of columns for input matrices A must not exceed the number
   of rows: N<=M.
---------------------------------------------------------------------------*/

void cholmxnn ( void    * pScr,
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
                int M,
                int N,
                int L )
{
    int16_t *         restrict R=(int16_t *      )_R;
    int16_t *         restrict D=(int16_t *      )_D;
    const int16_t * restrict   A=(const int16_t *)_A;
    static const tCholnIteration it[]=
    {
        {cholnConv0_7  ,cholnFwdrec0 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv0_7  ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv8_15 ,cholnFwdrec8 ,cholnDiagUpd8},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv8_15 ,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv16_23,cholnFwdrec16,cholnDiagUpd16},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv16_23,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv24_31,cholnFwdrec24,cholnDiagUpd24},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32},
        {cholnConv24_31,cholnFwdrec32,cholnDiagUpd32}
    };
    int32_t* Z=(int32_t*)pScr; // L columns of A'*A
    int16_t* y; // pointer to the new column in R
    int n;
    int SA=2*getSpace(M*N);
    int SR=2*getSpace((N*(N+1))>>1);
    int SD=2*getSpace(N);

    if (M<=0 || N<=0 || L<=0) return;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    y=R;

    for (n=0; n<XT_MIN(N,32); n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        it[n].conv(Z,A,n,M,N,L,SA);
        // make forward recursion to update n new column elements
        it[n].fwdrec(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        it[n].diagUpd(y,D+2*n,Z,n,L,SR,SD);
    }

    for (; n<N; n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        conv(Z,A,n,M,N,L,SA);
        // make forward recursion to update n new column elements
        cholnFwdrec(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        cholnDiagUpd(y,D+2*n,Z,n,L,SR,SD);
    }
} /* cholmxnn() */

#endif

#define CHOLN_SCRATCH(M,N,L) (2*(N)*(L)*sizeof(int32_t))

size_t cholmxnn_getScratchSize  (int M,int N, int L)
{
    M=M<0?0:M;
    N=N<0?0:N;
    L=L<0?0:L;
    return CHOLN_SCRATCH(M,N,L);
}

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
  NatureDSP_Baseband library. QR-based matrix decomposition and inversion for streaming order
    cqr_build_rMxNs
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
#include "qr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, qr_build_r8x8s, (void* pScr,int16_t * restrict V, int16_t * restrict R, int L))
#else

/*-----------------------------------------------------------------------
[c]qr_build_rMxNs

QR decomposition of MxN complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Fixed point representation of output matrices R is the same as for
input matrices, but Householder vectors V are always Q14.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]                  matrices A (L matrices of size MxN)
Output:
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
R[M*N][L]                  upper triangular matrices (L matrices
                           of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
void  qr_build_r8x8s ( void* pScr, int16_t * restrict V, int16_t * restrict R, int L)
{
    int16_t * pV = V;
    int16_t* pR;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH == 0);
    (void)pScr;

    pR = R;
    qrHouseholder8(pR, pV, 8, 8, L); qrUpdateR8(pR, pV, 0, 8, 8, L); pR += (8 + 1)*L; pV += 8 * BBE_SIMD_WIDTH;
    qrHouseholder7(pR, pV, 8, 8, L); qrUpdateR7(pR, pV, 1, 8, 8, L); pR += (8 + 1)*L; pV += 7 * BBE_SIMD_WIDTH;
    qrHouseholder6(pR, pV, 8, 8, L); qrUpdateR6(pR, pV, 2, 8, 8, L); pR += (8 + 1)*L; pV += 6 * BBE_SIMD_WIDTH;
    qrHouseholder5(pR, pV, 8, 8, L); qrUpdateR5(pR, pV, 3, 8, 8, L); pR += (8 + 1)*L; pV += 5 * BBE_SIMD_WIDTH;
    qrHouseholder4(pR, pV, 8, 8, L); qrUpdateR4(pR, pV, 4, 8, 8, L); pR += (8 + 1)*L; pV += 4 * BBE_SIMD_WIDTH;
    qrHouseholder3(pR, pV, 8, 8, L); qrUpdateR3(pR, pV, 5, 8, 8, L); pR += (8 + 1)*L; pV += 3 * BBE_SIMD_WIDTH;
    qrHouseholder2(pR, pV, 8, 8, L); qrUpdateR2(pR, pV, 6, 8, 8, L); pR += (8 + 1)*L; pV += 2 * BBE_SIMD_WIDTH;
} /* qr_build_r8x8s() */
#endif

size_t qr_build_r8x8s_getScratchSize (int M, int N,int L)
{
    (void)M;
    (void)N;
    (void)L;
    return 0;
} /* qr_build_r8x8s_getScratchSize() */

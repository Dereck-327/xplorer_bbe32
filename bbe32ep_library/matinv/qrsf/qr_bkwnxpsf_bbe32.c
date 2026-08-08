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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, real data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#if HAVE_VFPU
/*
    backward recursion: P==1
*/
static void qrsfBkwnx1(
                  float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            int N,int P,int L)
{
#if 0
    int m, k, _PL = P*L;
    int l;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));

    for (k = N - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t* pRt = R + (k*N + (k + 1))*L + l;
            // calculate y(m,:)-R(m,:)*X, 1xP
            float32_t r_re, x_re, B_re;
            B_re = (x[l + k*_PL + 0]);
            for (m = 0; m < N - k - 1; m++)
            {
                x_re = x[l + (k + 1 + m)*_PL];
                r_re = pRt[m*L];
                B_re -= (x_re*r_re);
            }
            x[l + k*_PL] = B_re*D[l + k*L];
        }
    }
#endif // 0

    int m, k, l;

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pX0;
    const xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pR0;
          xb_vecN_2xf32 * restrict pXw;

    xb_vecN_2xf32 Acc, X0, R0, D0;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));

    pD = (const xb_vecN_2xf32 *)(D + L*(N - 1));
    pX0 = (const xb_vecN_2xf32 *)(x + P*L*(N - 1));
    pXw = (xb_vecN_2xf32 *)(x + P*L*(N - 1));
    for (k = N - 1; k >= 0; k--)
    {
        pR0 = (const xb_vecN_2xf32 *)(R + L*(k*N + (k + 1)));
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            pX = pX0;
            pR = pR0;
            BBE_LVN_2XF32_XP(Acc, pX, P*L * sizeof(float32_t));

            for (m = 0; m < (N - k - 1); m++)
            {
                BBE_LVN_2XF32_XP(X0, pX, P*L * sizeof(float32_t));
                BBE_LVN_2XF32_XP(R0, pR, L * sizeof(float32_t));
                BBE_MULSN_2XF32(Acc, X0, R0);
            }

            BBE_LVN_2XF32_IP(D0, pD, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
            pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pX0);
            pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pR0);
        }
        pD = (const xb_vecN_2xf32 *)XT_ADDX4(-2 * L, (uintptr_t)pD);
        pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(-(P + 1)*L, (uintptr_t)pX0);
        pXw = (xb_vecN_2xf32 *)XT_ADDX4(-(P + 1)*L, (uintptr_t)pXw);
    }
}

/*-------------------------------------------------------------------------
[c]qr_bkwNxPsf

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
D[N*L]     reciprocals of main diagonal written in a special format
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/
void qr_bkwnxpsf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, const int N, const int P, int L)
{
    int p;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);

    for (p=0; p<P; p++)
    {
        qrsfBkwnx1(B, R, D,N,P,L);
        B+=L;
    }
}

// scratch memory needed for cbkwxxxxxs functions
size_t qr_bkwnxpsf_getScratchSize(int N, int P, int L) 
{
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    L=XT_MAX(L,0);
    return 0; 
}
#else
DISCARD_FUN(void, qr_bkwnxpsf, (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, const int N, const int P, int L))
size_t qr_bkwnxpsf_getScratchSize(int N, int P, int L) 
{
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    L=XT_MAX(L,0);
    return 0; 
}
#endif

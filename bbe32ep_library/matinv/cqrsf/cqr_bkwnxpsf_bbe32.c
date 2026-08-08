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
  QR decomposition, floating point, complex data, stream format
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
static void cqrsfBkwnx1(
                  float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            int N,int P,int L)
{
#if 0
    int m, k, _2L = 2 * L, _2PL = 2 * P*L;
    int l;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));

    for (k = N - 1; k >= 0; k--)
    {
        for (l = 0; l < L; l++)
        {
            const float32_t* pRt = R + (k*N + (k + 1))*_2L + 2 * l;
            // calculate y(m,:)-R(m,:)*X, 1xP
            float32_t r_re, r_im;
            float32_t x_re, x_im;
            float32_t B_re, B_im;
            B_re = (x[l * 2 + k*_2PL + 0]);
            B_im = (x[l * 2 + k*_2PL + 1]);

            for (m = 0; m < N - k - 1; m++)
            {
                x_re = x[l * 2 + (k + 1 + m)*_2PL + 0];
                x_im = x[l * 2 + (k + 1 + m)*_2PL + 1];
                r_re = pRt[m*_2L + 0];
                r_im = pRt[m*_2L + 1];
                B_re -= (x_re*r_re) - (x_im*r_im);
                B_im -= (x_re*r_im) + (x_im*r_re);
            }
            x[l * 2 + k*_2PL + 0] = B_re*D[l * 2 + k*_2L + 0];
            x[l * 2 + k*_2PL + 1] = B_im*D[l * 2 + k*_2L + 1];
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

    xb_vecN_2xf32 Acc, Acc1, Acc2, Acc3, X0, R0, D0;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));

    pD = (const xb_vecN_2xf32 *)(D + 2 * L*(N - 1));
    pX0 = (const xb_vecN_2xf32 *)(x + 2 * P*L*(N - 1));
    pXw = (xb_vecN_2xf32 *)(x + 2 * P*L*(N - 1));
    for (k = N - 1; k >= 0; k--)
    {
        pR0 = (const xb_vecN_2xf32 *)(R + 2 * L*(k*N + (k + 1)));
        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 4))
        {
            pX = pX0;
            pR = pR0;
            BBE_LVN_2XF32_XP(Acc, pX, 2 * P*L * sizeof(float32_t));
            Acc1 = BBE_ZERON_2XF32();
            Acc2 = BBE_ZERON_2XF32();
            Acc3 = BBE_ZERON_2XF32();

            for (m = 0; m < (N - k - 1) >> 1; m++)
            {
                BBE_LVN_2XF32_XP(X0, pX, 2 * P*L * sizeof(float32_t));
                BBE_LVN_2XF32_XP(R0, pR, 2 * L * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
                BBE_MULMASN_2XF32(Acc1, X0, R0, 2, 11);

                BBE_LVN_2XF32_XP(X0, pX, 2 * P*L * sizeof(float32_t));
                BBE_LVN_2XF32_XP(R0, pR, 2 * L * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc2, X0, R0, 3, 4);
                BBE_MULMASN_2XF32(Acc3, X0, R0, 2, 11);
            }
            if ((N - k - 1) & 1)
            {
                BBE_LVN_2XF32_XP(X0, pX, 2 * P*L * sizeof(float32_t));
                BBE_LVN_2XF32_XP(R0, pR, 2 * L * sizeof(float32_t));
                BBE_MULMASN_2XF32(Acc, X0, R0, 3, 4);
                BBE_MULMASN_2XF32(Acc1, X0, R0, 2, 11);
            }
            Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
            Acc = BBE_ADDN_2XF32(Acc, Acc1);
            Acc = BBE_ADDN_2XF32(Acc, Acc2);

            BBE_LVN_2XF32_IP(D0, pD, 2 * BBE_SIMD_WIDTH);
            Acc = BBE_MULN_2XF32(Acc, D0);
            BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
            pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pX0);
            pR0 = (const xb_vecN_2xf32 *)XT_ADDX4(8, (uintptr_t)pR0);
        }
        pD = (const xb_vecN_2xf32 *)XT_ADDX4(-4 * L, (uintptr_t)pD);
        pX0 = (const xb_vecN_2xf32 *)XT_ADDX4(-2 * (P + 1)*L, (uintptr_t)pX0);
        pXw = (xb_vecN_2xf32 *)XT_ADDX4(-2 * (P + 1)*L, (uintptr_t)pXw);
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
void cqr_bkwnxpsf  (void* pScr, complex_float* restrict _B, const complex_float* restrict _R, const complex_float* restrict _D, const int N, const int P, int L)
{
    float32_t *B=(float32_t *)_B;
    float32_t *R=(float32_t *)_R;
    float32_t *D=(float32_t *)_D;
    int p;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);

    for (p=0; p<P; p++)
    {
        cqrsfBkwnx1(B, R, D,N,P,L);
        B+=2*L;
    }
}
// scratch memory needed for cbkwxxxxxs functions
size_t cqr_bkwnxpsf_getScratchSize(int N, int P, int L) 
{
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    L=XT_MAX(L,0);
    return 0; 
}
#else
DISCARD_FUN(void, cqr_bkwnxpsf, (void* pScr, complex_float* restrict _B, const complex_float* restrict _R, const complex_float* restrict _D, const int N, const int P, int L))
// scratch memory needed for cbkwxxxxxs functions
size_t cqr_bkwnxpsf_getScratchSize(int N, int P, int L) 
{
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    L=XT_MAX(L,0);
    return 0; 
}
#endif

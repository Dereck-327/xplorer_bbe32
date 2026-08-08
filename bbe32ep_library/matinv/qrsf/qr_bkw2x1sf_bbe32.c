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
static void qrsfBkw2x1(
                  float32_t* restrict x, 
            const float32_t* restrict R,
            const float32_t* restrict D,
            int N,int P,int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pX;
    const xb_vecN_2xf32 * restrict pD;
    const xb_vecN_2xf32 * restrict pR;
          xb_vecN_2xf32 * restrict pXw;

    xb_vecN_2xf32 Acc, X0, R0, D0;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT(N == 2 && P == 1);

    pD = (const xb_vecN_2xf32 *)(D);
    pX = (const xb_vecN_2xf32 *)(x);
    pXw = (xb_vecN_2xf32 *)(x);
    pR = (const xb_vecN_2xf32 *)(R + L);

    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        X0 = BBE_LVN_2XF32_X(pX, L * sizeof(float32_t));
        D0 = BBE_LVN_2XF32_X(pD, L * sizeof(float32_t));
        X0 = BBE_MULN_2XF32(X0, D0);
        BBE_SVN_2XF32_X(X0, pXw, L * sizeof(float32_t));

        BBE_LVN_2XF32_IP(Acc, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(R0, pR, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(Acc, X0, R0);

        BBE_LVN_2XF32_IP(D0, pD, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULN_2XF32(Acc, D0);
        BBE_SVN_2XF32_IP(Acc, pXw, 2 * BBE_SIMD_WIDTH);
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
void qr_bkw2x1sf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D,int L)
{
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    qrsfBkw2x1(B, R, D,2,1,L);
}
size_t qr_bkw2x1sf_getScratchSize(int N, int P, int L) 
{ 
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(N==2 && P==1);
    L=XT_MAX(L,0);
    return 0; 
}
#else
DISCARD_FUN(void, qr_bkw2x1sf, (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D,int L))
size_t qr_bkw2x1sf_getScratchSize(int N, int P, int L) 
{ 
    (void)N,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(N==2 && P==1);
    L=XT_MAX(L,0);
    return 0; 
}
#endif

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

    LU decomposition for real matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "lunf_common.h"

#if HAVE_VFPU

/*--------------------------------------------------
update rows below i-th in matrix NxN (8x8, 16x16 only)

Input:
norm[L]   - normalization value for diagonal Aji, j=i+1...N-1
N           matrix size (8 or 16)
L           number of matrices
Input/output:
A[(N-i)*N]  matrix - pointer to the i-th row
Returns 
none:
--------------------------------------------------*/
// specialized version for N-i==2
void lunf_update2(float32_t* A, const float32_t* restrict norm, int i, int N, int L)
{
    int l;
    int SA=N*N;
    vboolN_2 bk,bi;
    xtfloat* restrict pNorm;
    xb_vecN_2xf32 * restrict pAw;
    xb_vecN_2xf32 * restrict pAj;
    xb_vecN_2xf32 * restrict pAi;
    xb_vecN_2xf32 Aj,Ai,Aji,NORM;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N-i==2);
    pNorm=(xtfloat*)norm;
    bi=BBE_LTRN_2I(7) &~ BBE_LTRN_2I(6);    /* 1 in i-th position */
    bk=~BBE_LTRN_2I(7);                    /* 1 for k>i          */
    pAw=pAj=(xb_vecN_2xf32*)&A[N + (i&8)];
    pAi=(xb_vecN_2xf32 *)&A[(i&8)];
    lastinc=(SA*sizeof(float32_t) - (N*sizeof(float32_t) * (2-1)));
    for (l=0; l<L; l++)
    {
        BBE_LVN_2XF32_XP(Ai,pAi,SA*sizeof(float32_t));
        BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
        NORM=BBE_REPN_2XF32(NORM,0);
        BBE_LVN_2XF32_XP(Aj,pAj,lastinc+N*sizeof(float32_t));
        BBE_MULN_2XF32T(Aj,Aj,NORM,bi);
        Aji=BBE_REPN_2XF32(Aj,6);
        BBE_MULSN_2XF32T(Aj,Aji,Ai,bk);
        BBE_SVN_2XF32_XP(Aj,pAw,lastinc+N*sizeof(float32_t));
    }
}
#endif

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
// specialized version N-i==12
void lunf_update12(float32_t* A, const float32_t* restrict norm, int i, int N, int L)
{
    int l;
    int SA=N*N;
    vboolN_2 bk,bi;
    xtfloat* restrict pNorm;
    xb_vecN_2xf32 * restrict pAw;
    xb_vecN_2xf32 * restrict pAj;
    xb_vecN_2xf32 * restrict pAi;
    xb_vecN_2xf32 Aj0,Aj1,Ai0,Ai1,Aji,NORM;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16);
    NASSERT(N-i==12);
    pNorm=(xtfloat*)norm;
    bi=BBE_LTRN_2I(5) &~ BBE_LTRN_2I(4);    /* 1 in i-th position */
    bk=~BBE_LTRN_2I(5);                    /* 1 for k>i          */
    pAw=pAj=(xb_vecN_2xf32*)&A[N ];
    pAi=(xb_vecN_2xf32 *)&A[0];
    lastinc=(SA*sizeof(float32_t) - (N*sizeof(float32_t) * (12-1)));
    for (l=0; l<L; l++)
    {
        int j;
        Ai1=BBE_LVN_2XF32_I (pAi,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(Ai0,pAi,SA*sizeof(float32_t));
        BBE_LSN_2XF32_IP(NORM,pNorm,sizeof(float32_t));
        NORM=BBE_REPN_2XF32(NORM,0);
        for (j=0; j<11; j++)
        {
            Aj1=BBE_LVN_2XF32_I (pAj,2*BBE_SIMD_WIDTH);
            BBE_LVN_2XF32_XP(Aj0,pAj,N*sizeof(float32_t));
            BBE_MULN_2XF32T(Aj0,Aj0,NORM,bi);
            Aji=BBE_REPN_2XF32(Aj0,4);
            BBE_MULSN_2XF32T(Aj0,Aji,Ai0,bk);
            BBE_MULSN_2XF32 (Aj1,Aji,Ai1);
            BBE_SVN_2XF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_XP(Aj0,pAw,N*sizeof(float32_t));
        }
        pAw+=lastinc/sizeof(xb_vecN_2xf32);
        pAj+=lastinc/sizeof(xb_vecN_2xf32);
    }
}
#endif

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

    LU decomposition for complex matrices (block ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "clunf_common.h"

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
// specialized version N-i==16
void clunf_update16(complex_float* A, const complex_float* restrict norm, int i, int N, int L)
{
    int l;
    int SA=N*N;
    vboolN_4 bk,bi;
    xtcomplexfloat* restrict pNorm;
    xb_vecN_4xcf32 * restrict pAw;
    xb_vecN_4xcf32 * restrict pAj;
    xb_vecN_4xcf32 * restrict pAi;
    xb_vecN_4xcf32 Aji,Aj0,Aj1,Aj2,Aj3,Ai0,Ai1,Ai2,Ai3,NORM;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16 || N==8);
    NASSERT(N-i==16);
    pNorm=(xtcomplexfloat*)norm;
    bi=BBE_LTRN_4I(1) ;    /* 1 in i-th position */
    bk=~BBE_LTRN_4I(1);    /* 1 for k>i          */
    pAw=pAj=(xb_vecN_4xcf32*)&A[N + (i&~3) ];
    pAi=(xb_vecN_4xcf32 *)&A[(i&~3)];
    lastinc=(SA*sizeof(complex_float) - (N*sizeof(complex_float) * (16-1)));
    for (l=0; l<L; l++)
    {
        int j;
        Ai1=BBE_LVN_4XCF32_I (pAi,1*2*BBE_SIMD_WIDTH);
        Ai2=BBE_LVN_4XCF32_I (pAi,2*2*BBE_SIMD_WIDTH);
        Ai3=BBE_LVN_4XCF32_I (pAi,3*2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Ai0,pAi,SA*sizeof(complex_float));
        BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
        NORM=BBE_REPN_4XCF32(NORM,0);
        for (j=0; j<16-1; j++)
        {
            Aj1=BBE_LVN_4XCF32_I (pAj,1*2*BBE_SIMD_WIDTH);
            Aj2=BBE_LVN_4XCF32_I (pAj,2*2*BBE_SIMD_WIDTH);
            Aj3=BBE_LVN_4XCF32_I (pAj,3*2*BBE_SIMD_WIDTH);
            BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
            BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
            Aji=BBE_REPN_4XCF32(Aj0,0);
            BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
            BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
            BBE_MULSN_4XCF32 (Aj2,Aji,Ai2);
            BBE_MULSN_4XCF32 (Aj3,Aji,Ai3);
            BBE_SVN_4XCF32_I (Aj1,pAw,1*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I (Aj2,pAw,2*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_I (Aj3,pAw,3*2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));
        }
        pAw+=lastinc/sizeof(xb_vecN_4xcf32);
        pAj+=lastinc/sizeof(xb_vecN_4xcf32);
    }
}
#endif

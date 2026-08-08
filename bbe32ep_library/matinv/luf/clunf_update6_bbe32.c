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
// specialized version N-i==6
void clunf_update6(complex_float* A, const complex_float* restrict norm, int i, int N, int L)
#if 1
{
    int l;
    int SA=N*N;
    vboolN_4 bk,bi;
    xtcomplexfloat* restrict pNorm;
    xb_vecN_4xcf32 * restrict pAw;
    xb_vecN_4xcf32 * restrict pAj;
    xb_vecN_4xcf32 * restrict pAi;
    xb_vecN_4xcf32 Aj0,Aj1,Ai0,Ai1,Aji,NORM;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16 || N==8);
    NASSERT(N-i==6);
    pNorm=(xtcomplexfloat*)norm;
    bi=BBE_LTRN_4I(3) &~ BBE_LTRN_4I(2);    /* 1 in i-th position */
    bk=~BBE_LTRN_4I(3);                    /* 1 for k>i          */
    pAw=pAj=(xb_vecN_4xcf32*)&A[N+(i&~3)];
    pAi=(xb_vecN_4xcf32 *)&A[(i&~3)];
    lastinc=(SA*sizeof(complex_float) - (N*sizeof(complex_float) * (6-1)));
    for (l=0; l<L; l++)
    {
        int j;
        Ai1=BBE_LVN_4XCF32_I (pAi,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Ai0,pAi,SA*sizeof(complex_float));
        BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
        NORM=BBE_REPN_4XCF32(NORM,0);
        for (j=0; j<6-1; j++)
        {
            Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
            BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
            BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
            Aji=BBE_REPN_4XCF32(Aj0,2);
            BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
            BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
            BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
            BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));
        }
        pAw+=lastinc/sizeof(xb_vecN_4xcf32);
        pAj+=lastinc/sizeof(xb_vecN_4xcf32);
    }
}
#else
{
    int l;
    int SA=N*N;
    vboolN_4 bk,bi;
    xtcomplexfloat* restrict pNorm;
    xb_vecN_4xcf32 * restrict pAw;
    xb_vecN_4xcf32 * restrict pAj;
    xb_vecN_4xcf32 * restrict pAi;
    xb_vecN_4xcf32 Aj0,Aj1,Ai0,Ai1,Aji,NORM;
    int lastinc;

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(N==16 || N==8);
    NASSERT(N-i==6);
    pNorm=(xtcomplexfloat*)norm;
    bi=BBE_LTRN_4I(3) &~ BBE_LTRN_4I(2);    /* 1 in i-th position */
    bk=~BBE_LTRN_4I(3);                    /* 1 for k>i          */
    pAw=pAj=(xb_vecN_4xcf32*)&A[N + (i&~3)];
    pAi=(xb_vecN_4xcf32 *)&A[(i&~3)];
    lastinc=(SA*sizeof(complex_float) - (N*sizeof(complex_float) * (6-1)));
    for (l=0; l<L; l++)
    {
        Ai1=BBE_LVN_4XCF32_I (pAi,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Ai0,pAi,SA*sizeof(complex_float));
        BBE_LSN_4XCF32_IP(NORM,pNorm,sizeof(complex_float));
        NORM=BBE_REPN_4XCF32(NORM,0);

        Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
        BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
        Aji=BBE_REPN_4XCF32(Aj0,2);
        BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
        BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
        BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));

        Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
        BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
        Aji=BBE_REPN_4XCF32(Aj0,2);
        BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
        BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
        BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));

        Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
        BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
        Aji=BBE_REPN_4XCF32(Aj0,2);
        BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
        BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
        BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));

        Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Aj0,pAj,N*sizeof(complex_float));
        BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
        Aji=BBE_REPN_4XCF32(Aj0,2);
        BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
        BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
        BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(Aj0,pAw,N*sizeof(complex_float));

        Aj1=BBE_LVN_4XCF32_I (pAj,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_XP(Aj0,pAj,lastinc+N*sizeof(complex_float));
        BBE_MULN_4XCF32T(Aj0,Aj0,NORM,bi);
        Aji=BBE_REPN_4XCF32(Aj0,2);
        BBE_MULSN_4XCF32T(Aj0,Aji,Ai0,bk);
        BBE_MULSN_4XCF32 (Aj1,Aji,Ai1);
        BBE_SVN_4XCF32_I (Aj1,pAw,2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(Aj0,pAw,lastinc+N*sizeof(complex_float));
    }
}
#endif
#endif

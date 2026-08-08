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
    BBE16EP code for Cholesky decomposition, block format
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

#if HAVE_CHOLN

#if 0
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

static void conv(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA)
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
#endif
/*---------------------------------------------------
   compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   n            number of column
   Output:
   Z[L][n+1][2] results
---------------------------------------------------*/
// n=16...23
void cholnConv16_23(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA)
{
#if 1
    int l, j;
    const xb_vecNx16 *restrict pA0;
    const xb_vecNx16 *restrict pA1;
          xb_vecNx16 *restrict pZ = (      xb_vecNx16*)Z;

    xb_vecNx16 A00, A01, A02, A10, A11, A12;
    xb_vecNx16 repA0, repA1;
    xb_vecNx16 Zl, Zh;
    xb_vecNx40 acc0, acc1, acc2;
    valign alignA, align;
    const vselN sel= BBE_MOVVSELNX16(BBE_MOVVA16C(((2*(n-16)+1) << 16) + 2*(n-16)),0);
    NASSERT(n>=16 && n<=23);
    NASSERT(N%4==0 && M%4==0 && SA==2*(N*M));
    (void)SA;

    align = BBE_ZALIGN();
    pA0 = (const xb_vecNx16*)A;
    __Pragma("loop_count min=1");
    for (l=0; l<L; l++)
    {
        acc0= acc1= acc2= BBE_ZERONX40();
        __Pragma("loop_count min=2, factor=2");
        for (j=0; j<(M>>1); j++)
        {
            pA1 = (const xb_vecNx16*)((int16_t *)(pA0)+N*2);
            BBE_LVNX16_IP(A00, pA0,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(A01, pA0,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(A02, pA0,N*2*sizeof(int16_t)*2-4*BBE_SIMD_WIDTH);
            
            alignA = BBE_LANX16_PP(pA1);
            BBE_LANX16_IP (A10, alignA, pA1);
            BBE_LANX16_IP (A11, alignA, pA1);
            BBE_LAVNX16_XP(A12, alignA, pA1, N*sizeof(int16_t)*2-4*BBE_SIMD_WIDTH);

            repA0 = BBE_SHFLNX16(A02,sel);
            repA1 = BBE_SHFLNX16(A12,sel);
            BBE_MULANX16J(acc0,repA0,A00);
            BBE_MULANX16J(acc1,repA0,A01);
            BBE_MULANX16J(acc2,repA0,A02);
            BBE_MULANX16J(acc0,repA1,A10);
            BBE_MULANX16J(acc1,repA1,A11);
            BBE_MULANX16J(acc2,repA1,A12);
        }
        // store res
        Zl=BBE_MOVSVWL(acc0);
        Zh=BBE_MOVSVWH(acc0);
        BBE_SANX16_IP(Zl,align,pZ);
        BBE_SANX16_IP(Zh,align,pZ);

        Zl=BBE_MOVSVWL(acc1);
        Zh=BBE_MOVSVWH(acc1);
        BBE_SANX16_IP(Zl,align,pZ);
        BBE_SANX16_IP(Zh,align,pZ);

        Zl=BBE_MOVSVWL(acc2);
        Zh=BBE_MOVSVWH(acc2);
        BBE_SAVNX16_XP(Zl,align,pZ,(n+1)*2*4 - 8*BBE_SIMD_WIDTH);
        BBE_SAVNX16_XP(Zh,align,pZ,(n+1)*2*4 - 10*BBE_SIMD_WIDTH);
    }
    BBE_SAPOS_FP(align,pZ);
#else
    conv(Z, A, n, M, N, L, SA);
#endif
}
#endif

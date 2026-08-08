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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    Cholesky backward recursion, block format, real floating point data
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
    */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "rcholnf_common.h"

#if (HAVE_VFPU)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-1)) m=LOG2_BBE_SIMD_WIDTH-1;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*----------------------------------------------------------------------------------
   reversing R matrices for easier readings by rows (diagonal elements are omitted):
   original R    transformed R
   0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
     2 4 7 b
       5 8 c
         9 d
           e

   Input:
   R[L][SR]        L input matrices
   Rt[L*N*(N-1)]   stream of L trasposed matrices
----------------------------------------------------------------------------------*/
void rcholnfTransformR(float32_t* Rt,const float32_t* R,int N,int L)
{
#if 0
    int SR = getSpace((N*(N + 1)) >> 1);
    int n, m, l;
    const float32_t* pR;
    for (l = 0; l < L; l++)
    {
        for (n = 0; n < N; n++)
        {
            pR = R + l*SR + (((N - n)*(N - n + 3) - 2) >> 1);
            for (m = 0; m < n; m++)
            {
                Rt[0] = pR[0];
                pR += 1 * (N - n + m + 1);
                Rt += 1;
            }
        }
    }
#endif // 0

    const int SR = getSpace((N*(N + 1)) >> 1);
    int n, m, l;
    const xtfloat * restrict pR;
          xtfloat * restrict pRw = (      xtfloat *)Rt;
    const xtfloat * restrict pR2 = (const xtfloat *)R;
    int delta, delta2, dpR;
    const int delta3 = (N + 1) * sizeof(float32_t);
    xb_vecN_2xf32 vTmp;

    __Pragma("loop_count min=1");
    for (l = 0; l < L; l++)
    {
        delta2 = delta3;
        dpR = N;
        __Pragma("loop_count min=4");
        for (n = 0; n < N; n++)
        {
            pR = (const xtfloat *)XT_ADDX2((dpR)*(dpR + 3) - 2, (uintptr_t)pR2);
            delta = delta2;
            for (m = 0; m < n; m++)
            {
                BBE_LSN_2XF32_XP(vTmp, pR, delta);
                BBE_SSN_2XF32_XP(vTmp, pRw, sizeof(float32_t));
                delta = XT_ADDX4(1, delta);
            }
            delta2 = XT_ADDX4(-1, delta2);
            dpR = XT_ADDI(dpR, -1);
        }
        pR2 = (const xtfloat *)XT_ADDX4(SR, (uintptr_t)pR2);
    }
}
#endif

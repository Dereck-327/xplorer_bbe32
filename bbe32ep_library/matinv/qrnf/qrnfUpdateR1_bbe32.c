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
  QR decomposition, floating point, real data, block format
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "qrnf_common.h"
#if HAVE_VFPU

/*-------------------------------------------------------
    partial update of R matrix
    Fi[L][SV] diagonal rotation matrix (only 0-th element filled)
    v[L][SV]  Housholder vector (M elements filled)
    SV,SD     strides for V/Fi and D
    M,N       matrix size
    L         number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    special case: N=M==1
-------------------------------------------------------*/
void qrnfUpdateR1(      float32_t* restrict R,
                  const float32_t* restrict v,
                  int SA, int L)
{
#if 0
    int l;
    float32_t A_re;
    // compute v'B first
    for (l = 0; l < L; l++)
    {
        float32_t v_re, b_re;
        v_re = v[l];
        b_re = R[l*SA];
        A_re = b_re*(1.f - 2.f*v_re*v_re);
        R[l*SA] = A_re;
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pV  = (const xb_vecN_2xf32 *)v;
    const xtfloat       * restrict pV_;
    const xtfloat       * restrict pR  = (const xtfloat       *)R;
          xtfloat       * restrict pRw = (      xtfloat       *)R;

    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;
    xb_vecN_2xf32 Acc, V0, R0, Z0;
    valign vV;

    vV = BBE_LAN_2XF32_PP(pV);
    // compute v'B first
    for (l = 0; l < L / (BBE_SIMD_WIDTH / 2); l++)
    {
        t7 = BBE_LSN_2XF32_X(pR, 7 * SA * sizeof(float32_t));
        t6 = BBE_LSN_2XF32_X(pR, 6 * SA * sizeof(float32_t));
        t5 = BBE_LSN_2XF32_X(pR, 5 * SA * sizeof(float32_t));
        t4 = BBE_LSN_2XF32_X(pR, 4 * SA * sizeof(float32_t));
        t3 = BBE_LSN_2XF32_X(pR, 3 * SA * sizeof(float32_t));
        t2 = BBE_LSN_2XF32_X(pR, 2 * SA * sizeof(float32_t));
        t1 = BBE_LSN_2XF32_X(pR,     SA * sizeof(float32_t));
        BBE_LSN_2XF32_XP(t0, pR, 8 * SA * sizeof(float32_t));

        t0 = BBE_SELN_2XF32I(t1, t0, BBE_SELI_INTERLEAVE_2_LO);
        t2 = BBE_SELN_2XF32I(t3, t2, BBE_SELI_INTERLEAVE_2_LO);
        t4 = BBE_SELN_2XF32I(t5, t4, BBE_SELI_INTERLEAVE_2_LO);
        t6 = BBE_SELN_2XF32I(t7, t6, BBE_SELI_INTERLEAVE_2_LO);
        t0 = BBE_SELN_2XF32I(t2, t0, BBE_SELI_INTERLEAVE_4_LO);
        t4 = BBE_SELN_2XF32I(t6, t4, BBE_SELI_INTERLEAVE_4_LO);
        R0 = BBE_SELN_2XF32I(t4, t0, BBE_SELI_EXTRACT_LO_HALVES);

        BBE_LAN_2XF32_IP(V0, vV, pV);

        Acc = BBE_MULN_2XF32(V0, V0);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);
        BBE_MULSN_2XF32(R0, Z0, R0);

        t0 = BBE_REPN_2XF32(R0, 0);
        t1 = BBE_REPN_2XF32(R0, 1);
        t2 = BBE_REPN_2XF32(R0, 2);
        t3 = BBE_REPN_2XF32(R0, 3);
        t4 = BBE_REPN_2XF32(R0, 4);
        t5 = BBE_REPN_2XF32(R0, 5);
        t6 = BBE_REPN_2XF32(R0, 6);
        t7 = BBE_REPN_2XF32(R0, 7);

        BBE_SSN_2XF32_X(t7, pRw, 7 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t6, pRw, 6 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t5, pRw, 5 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t4, pRw, 4 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t3, pRw, 3 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t2, pRw, 2 * SA * sizeof(float32_t));
        BBE_SSN_2XF32_X(t1, pRw,     SA * sizeof(float32_t));
        BBE_SSN_2XF32_XP(t0, pRw, 8 * SA * sizeof(float32_t));
    }

    pV_ = (const xtfloat *)pV;
    __Pragma("loop_count max=7");
    for (l *= (BBE_SIMD_WIDTH / 2); l < L; l++)
    {
        BBE_LSN_2XF32_XP(R0, pR, SA * sizeof(float32_t));
        BBE_LSN_2XF32_XP(V0, pV_, sizeof(float32_t));

        Acc = BBE_MULN_2XF32(V0, V0);
        Z0 = BBE_ADDN_2XF32(R0, R0);
        BBE_MULSN_2XF32(R0, Z0, Acc);

        BBE_SSN_2XF32_XP(R0, pRw, SA * sizeof(float32_t));
    }
}
#endif

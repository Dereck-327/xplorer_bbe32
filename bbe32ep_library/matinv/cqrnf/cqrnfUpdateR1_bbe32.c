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
  QR decomposition, floating point, complex data, block format
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "cqrnf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------
    partial update of R matrix
    Fi[L][SV][2] diagonal rotation matrix (only 0-th element filled)
    v[L][SV][2]  Housholder vector (M elements filled)
    SV,SD        strides for V/Fi and D
    M,N          matrix size
    L            number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    special case: N=M==1
-------------------------------------------------------*/
void cqrnfUpdateR1(       float32_t* restrict R,
                    const float32_t* restrict v,
                    int SA, int L)
{
#if 0
    int l;
    float32_t A_re, A_im;
    // compute v'B first
    for (l = 0; l < L; l++)
    {
        float32_t v_re, v_im, b_re, b_im, z_re, z_im;
        v_re = v[2 * l + 0];  v_im = v[2 * l + 1];
        b_re = R[l*SA + 0]; b_im = R[l*SA + 1];
        z_re = (v_re*b_re) + (v_im*b_im);
        z_im = (v_re*b_im) - (v_im*b_re);
        z_re *= 2.f;
        z_im *= 2.f;
        A_re = b_re - ((z_re*v_re) - (z_im*v_im));
        A_im = b_im - ((z_re*v_im) + (z_im*v_re));
        R[l*SA + 0] = A_re;
        R[l*SA + 1] = A_im;
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pV  = (const xb_vecN_2xf32 *)v;
    const long long     * restrict pV_;
    const long long     * restrict pR  = (const long long     *)R;
          long long     * restrict pRw = (      long long     *)R;

    xb_vecN_4x64 t0, t1, t2, t3;
    xb_vecN_2xf32 Acc, /*Acc1, Acc2, Acc3,*/ V0, R0, Z0;
    valign vV;

    vV = BBE_LAN_2XF32_PP(pV);
    // compute v'B first
    for (l = 0; l < L / (BBE_SIMD_WIDTH / 4); l++)
    {
        t3 = BBE_LSN_4X64_X(pR, 3 * SA * sizeof(float32_t));
        t2 = BBE_LSN_4X64_X(pR, 2 * SA * sizeof(float32_t));
        t1 = BBE_LSN_4X64_X(pR,     SA * sizeof(float32_t));
        BBE_LSN_4X64_XP(t0, pR, 4 * SA * sizeof(float32_t));

        t0 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t1), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_INTERLEAVE_4_LO));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t3), BBE_MOVNX16_FROMN_4X64(t2), BBE_SELI_INTERLEAVE_4_LO));
        R0 = BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16I(BBE_MOVNX16_FROMN_4X64(t2), BBE_MOVNX16_FROMN_4X64(t0), BBE_SELI_EXTRACT_LO_HALVES));

        BBE_LAN_2XF32_IP(V0, vV, pV);

        Acc = BBE_MULMN_2XF32(V0, R0, 0, 4);
        BBE_MULMASN_2XF32(Acc, V0, R0, 2, 11);

        Z0 = BBE_ADDN_2XF32(Acc, Acc);
        //Z0 = BBE_MULN_2XF32(Acc, 2.f);

        BBE_MULMASN_2XF32(R0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(R0, Z0, V0, 2, 11);

        t0 = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(R0));
        t1 = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(R0, BBE_SHFLI_REP_1X4)));
        t2 = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(R0, BBE_SHFLI_REP_2X4)));
        t3 = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BBE_SHFLN_2XF32I(R0, BBE_SHFLI_REP_3X4)));

        BBE_SSN_4X64_X(t3, pRw, 3 * SA * sizeof(float32_t));
        BBE_SSN_4X64_X(t2, pRw, 2 * SA * sizeof(float32_t));
        BBE_SSN_4X64_X(t1, pRw,     SA * sizeof(float32_t));
        BBE_SSN_4X64_XP(t0, pRw, 4 * SA * sizeof(float32_t));
    }

    pV_ = (const long long *)pV;
    __Pragma("loop_count max=3");
    for (l *= (BBE_SIMD_WIDTH / 4); l < L; l++)
    {
        BBE_LSN_4X64_XP(t0, pR, SA * sizeof(float32_t));
        R0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(t0));

        BBE_LSN_4X64_XP(t0, pV_, 2 * sizeof(float32_t));
        V0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4X64(t0));

        Acc = BBE_MULMN_2XF32(V0, R0, 0, 4);
        BBE_MULMASN_2XF32(Acc, V0, R0, 2, 11);

        Z0 = BBE_ADDN_2XF32(Acc, Acc);
        //Z0 = BBE_MULN_2XF32(Acc, 2.f);

        BBE_MULMASN_2XF32(R0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(R0, Z0, V0, 2, 11);

        t0 = BBE_MOVN_4X64_FROMNX16(BBE_MOVNX16_FROMN_2XF32(R0));

        BBE_SSN_4X64_XP(t0, pRw, SA * sizeof(float32_t));
    }
}
#endif

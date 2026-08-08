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
#include "common.h"
#include "qrsf_common.h"

#if HAVE_VFPU

/*---------------------------------------------------------
    update matrix B[L][MxP] by housholder vectors V[L][SV]

    Input:
    M,N,P,L     dimensions
    V[L][SV]    Householder vectors
    Input/output:
    B[L][MxP]   B matrices MxP
  
---------------------------------------------------------*/
//P==1, M=3
void qrsfUpdateB3(float32_t* B,const float32_t* V,int L)
{
    int l;

    const xb_vecN_2xf32 * restrict pBr;
          xb_vecN_2xf32 * restrict pB;
    const xb_vecN_2xf32 * restrict pV;

    xb_vecN_2xf32 Acc, V0, V1, V2, B0, B1, B2, Z0;

    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 2) == 0));

    pBr = (const xb_vecN_2xf32 *)B;
    pB = (xb_vecN_2xf32 *)B;
    pV = (const xb_vecN_2xf32 *)V;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
    {
        V2 = BBE_LVN_2XF32_X(pV, 2 * L * sizeof(float32_t));
        B2 = BBE_LVN_2XF32_X(pBr, 2 * L * sizeof(float32_t));
        Acc = BBE_MULN_2XF32(V2, B2);
        V1 = BBE_LVN_2XF32_X(pV, L * sizeof(float32_t));
        B1 = BBE_LVN_2XF32_X(pBr, L * sizeof(float32_t));
        BBE_MULAN_2XF32(Acc, V1, B1);
        BBE_LVN_2XF32_IP(V0, pV, 2 * BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(B0, pBr, 0);
        BBE_MULAN_2XF32(Acc, V0, B0);

        Z0 = BBE_ADDN_2XF32(Acc, Acc);
        //Z0 = BBE_MULN_2XF32(Acc, 2.f);
            
        B2 = BBE_LVN_2XF32_X(pBr, 2 * L * sizeof(float32_t));
        BBE_MULSN_2XF32(B2, Z0, V2);
        BBE_SVN_2XF32_X(B2, pB, 2 * L * sizeof(float32_t));
        B1 = BBE_LVN_2XF32_X(pBr, L * sizeof(float32_t));
        BBE_MULSN_2XF32(B1, Z0, V1);
        BBE_SVN_2XF32_X(B1, pB, L * sizeof(float32_t));
        BBE_LVN_2XF32_IP(B0, pBr, 2 * BBE_SIMD_WIDTH);
        BBE_MULSN_2XF32(B0, Z0, V0);
        BBE_SVN_2XF32_IP(B0, pB, 2 * BBE_SIMD_WIDTH);
    }
}
#endif

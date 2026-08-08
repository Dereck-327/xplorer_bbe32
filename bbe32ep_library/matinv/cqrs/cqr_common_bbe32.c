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
  NatureDSP_Baseband library. QR-based matrix decomposition and inversion for streaming order
    complex qr decompostion/inversion helper routine
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, cqrComputeLastRot,(int16_t* pFi, const int16_t* R, int M, int L))
#else

/*
    compute last rotation for square matrices
*/
void cqrComputeLastRot(int16_t* pFi, const int16_t* R, int M, int L)
{
    const xb_vecNx16* restrict pR;
    xb_vecNx16* restrict _pFi;
    int l;
    pR = (const xb_vecNx16*)(R + (M*M - 1)*L * 2);
    _pFi = (xb_vecNx16*)XT_ADD(2 * BBE_SIMD_WIDTH*(M - 1), (uintptr_t)pFi);
    for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        xb_vecNx16  x0, _fi, mant, tmp0, tmp1;
        vsaN exp, c_vec;
        xb_vecNx40 acc0, temp;
        BBE_LVNX16_IP(x0, pR, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MAGINX16C(x0, x0);
        acc0 = BBE_ADDNX40(acc0, acc0);
        c_vec = BBE_NSAENX40(acc0);
        acc0 = BBE_SLLNX40(acc0, c_vec);
        BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
        BBE_MULUUSNX16(acc0, tmp1, tmp0);
        acc0 = BBE_SRAINX40(acc0, 24);
        tmp0 = BBE_PACKLNX40(acc0);
        mant = BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
        exp = BBE_SUBSR1SAVSN(18 + 1, c_vec);
        temp = BBE_MULUSRNX16(mant, x0, exp);
        _fi = BBE_PACKVNX40(temp, exp);
        BBE_SVNX16_XP(_fi, _pFi, 2 * BBE_SIMD_WIDTH*(SIZE_OF_V(M, M) + SIZE_OF_FI(M)));
    }
}
#endif

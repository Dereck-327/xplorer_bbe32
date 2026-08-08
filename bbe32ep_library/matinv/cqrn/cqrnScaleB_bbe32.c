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
    BBE32 code for QR decomposition (calc_qb part), block format
    IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqrn_common.h"

#if HAVE_CQRN

/*-------------------------------------------------------------------------
  Scale B[SB] by 1 bit left with rounding
  Input/output:
  B[SB]        input/output matrix
---------------------------------------------------------------------------*/
void cqrnScaleB(int16_t* B,int SB)
{
    const xb_vecNx16 *restrict pB_ld;
          xb_vecNx16 *restrict pB_st;
    xb_vecNx16 b0, c0;
    int m;

    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT(SB%8==0);

    pB_st = (xb_vecNx16 *)B;
    pB_ld = pB_st;
    c0 = BBE_ZERONX16();
    /* scale down B by 1 bit right with rounding */
    __Pragma("loop_count min=1")
    for (m=0; m<(SB>>LOG2_BBE_SIMD_WIDTH); m++)
    {
        BBE_LVNX16_IP(b0, pB_ld, sizeof(int16_t)*BBE_SIMD_WIDTH);
        b0 = BBE_ADDSR1RNX16(b0, c0);
        BBE_SVNX16_IP(b0, pB_st, sizeof(int16_t)*BBE_SIMD_WIDTH);
    }
    if (SB&(BBE_SIMD_WIDTH/2))
    {
        valign al_ld, al_st;
        al_ld = BBE_LANX16_PP(pB_ld);
        al_st = BBE_ZALIGN();
        BBE_LAVNX16_XP(b0, al_ld, pB_ld, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
        b0 = BBE_ADDSR1RNX16(b0, c0);
        BBE_SAVNX16_XP(b0, al_st, pB_st, sizeof(int16_t)*(BBE_SIMD_WIDTH/2));
        BBE_SANX16POS_FP(al_st, pB_st);
    }
}

#endif

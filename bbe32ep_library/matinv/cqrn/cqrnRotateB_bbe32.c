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

#define sz_i16 sizeof(int16_t)

/*------------------------------------------------
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    introduces additional shift right !
    Input:
    Fi[L][SV]
    Input/output:
    B[L][SB]
    Temporary:
    pScr    - size in bytes N*2*BBE_SIMD_WIDTH
------------------------------------------------*/
void cqrnRotateB(int16_t* B,const int16_t* Fi,int N,int P,int SB,int L)
{
    const xb_vecNx16 * restrict pFi;
    const xb_vecNx16 * restrict pB_ld0;
    const xb_vecNx16 * restrict pB_ld1;
          xb_vecNx16 * restrict pB_st0;
          xb_vecNx16 * restrict pB_st1;
    xb_vecNx16 B0, B1, F0, F1;
    xb_vecNx40 ACC0, ACC1;
    valign al_ld0, al_ld1, al_st0, al_st1;
    vsaN rnd;
    int NbPtail;
    int l, n, p;

    NASSERT(L>0);
    NASSERT(P>1);
    NASSERT(N%4==0 && N>0);

    rnd = BBE_MOVVSA32(14);
    NbPtail = 2*sz_i16*(P & (BBE_SIMD_WIDTH/2-1));
    al_st0 = al_st1 = BBE_ZALIGN();

    __Pragma("loop_count min=1")
    for (l = 0; l < L; l++)
    {
        pFi = (const xb_vecNx16 *)(Fi + 2*l);
        pB_st0 = (   xb_vecNx16 *)(B + SB*l);
        
        /* Process matrix B by 2 rows per iteration */
        __Pragma("loop_count min=2")
        for (n = 0; n < (N >> 1); n++)
        {
            pB_st1 = (xb_vecNx16 *)((int16_t *)pB_st0 + 2*P);
            pB_ld0 = pB_st0;
            pB_ld1 = pB_st1;
            al_ld0 = BBE_LANX16_PP(pB_ld0);
            al_ld1 = BBE_LANX16_PP(pB_ld1);
            /* load rotation values for 2 rows */
            BBE_LPNX16_XP(F0, pFi, sz_i16*2*L);
            BBE_LPNX16_XP(F1, pFi, sz_i16*2*L);
            F0 = BBE_REPNX16C(F0, 0);
            F1 = BBE_REPNX16C(F1, 0);
            
            /* Innermost loop: process by 2*(SIMD_WIDTH/2) per iteration */
            for (p = 0; p < (P >> (LOG2_BBE_SIMD_WIDTH-1)); p++)
            {
                BBE_LANX16_IP(B0, al_ld0, pB_ld0);
                BBE_LANX16_IP(B1, al_ld1, pB_ld1);

                ACC0 = BBE_MULRNX16J(B0, F0, rnd);
                ACC1 = BBE_MULRNX16J(B1, F1, rnd);

                B0 = BBE_PACKVNX40(ACC0, rnd);
                B1 = BBE_PACKVNX40(ACC1, rnd);

                BBE_SANX16_IP(B0, al_st0, pB_st0);
                BBE_SANX16_IP(B1, al_st1, pB_st1);
            }
            /* process last P%(SIMD_WIDTH/2) values of 2 rows */
            BBE_LAVNX16_XP(B0, al_ld0, pB_ld0, NbPtail);
            BBE_LAVNX16_XP(B1, al_ld1, pB_ld1, NbPtail);

            ACC0 = BBE_MULRNX16J(B0, F0, rnd);
            ACC1 = BBE_MULRNX16J(B1, F1, rnd);

            B0 = BBE_PACKVNX40(ACC0, rnd);
            B1 = BBE_PACKVNX40(ACC1, rnd);
            /* save updated values */
            BBE_SAVNX16_XP(B0, al_st0, pB_st0, NbPtail);
            BBE_SAVNX16_XP(B1, al_st1, pB_st1, NbPtail);
            BBE_SANX16POS_FP(al_st0, pB_st0);
            BBE_SANX16POS_FP(al_st1, pB_st1);
            /* jump to the next 2 rows */
            pB_st0 = pB_st1;
        }
    }
}
#endif

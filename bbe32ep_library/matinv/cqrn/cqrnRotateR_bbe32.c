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
    BBE32 code for QR decomposition 
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
/*-------------------------------------------------------
    rotate R[L][SA] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV] diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    extraShift extra shift right
    Input/output:
    R[L][SA]  sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper trinagle matrix 
              are zeroed!
-------------------------------------------------------*/
void cqrnRotateR(int16_t* restrict R,const int16_t* restrict Fi,int N,int SA,int L, int extraShift)
{
    const xb_vecNx16 * restrict pFi;
    const xb_vecNx16 * restrict pR_ld0;
    const xb_vecNx16 * restrict pR_ld1;
    const xb_vecNx16 * restrict pR_ld2;
    const xb_vecNx16 * restrict pR_ld3;
          xb_vecNx16 * restrict pR_st0;
          xb_vecNx16 * restrict pR_st1;
          xb_vecNx16 * restrict pR_st2;
          xb_vecNx16 * restrict pR_st3;
    xb_vecNx16 R0, R1, R2, R3, F0, F1, F2, F3;
    xb_vecNx40 ACC0, ACC1, ACC2, ACC3;
    vsaN rnd;

    int l, m, n;

    rnd = BBE_MOVVSA32(15-extraShift);

    if (N & 4)
    {
        valign al_R_ld0, al_R_ld1, al_R_st0, al_R_st1;
        al_R_st0 = al_R_st1 = BBE_ZALIGN();
        __Pragma("loop_count min=1")
        for (l = 0; l < L; l++)
        {
            pFi = (const xb_vecNx16 *)(Fi + 2*l);

            pR_ld0 = (const xb_vecNx16 *)(R + SA*l);
            pR_ld1 = (const xb_vecNx16 *)((int16_t *)pR_ld0 + 2*N);

            /* Process matrix R by 2 rows per iteration */
            __Pragma("loop_count min=1")
            for (m = 0; m < (N >> 1); m++)
            {
                /* load rotation vectors for 2 rows */
                BBE_LPNX16_XP(F0, pFi, sz_i16*2*L);
                BBE_LPNX16_XP(F1, pFi, sz_i16*2*L);
                F0 = BBE_REPNX16C(F0, 0);
                F1 = BBE_REPNX16C(F1, 0);

                al_R_ld0 = BBE_LANX16_PP(pR_ld0);
                al_R_ld1 = BBE_LANX16_PP(pR_ld1);
                pR_st0 = (xb_vecNx16 *)pR_ld0;
                pR_st1 = (xb_vecNx16 *)pR_ld1;

                /* Innermost loop: process 4 rows by SIMD_WIDTH/2 values per iteration */
                for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH-1)); n++)
                {
                    /* load and rotate values */
                    BBE_LANX16_IP(R0, al_R_ld0, pR_ld0);
                    BBE_LANX16_IP(R1, al_R_ld1, pR_ld1);

                    ACC0 = BBE_MULRNX16J(R0, F0, rnd);
                    ACC1 = BBE_MULRNX16J(R1, F1, rnd);

                    R0 = BBE_PACKVNX40(ACC0, rnd);
                    R1 = BBE_PACKVNX40(ACC1, rnd);
                    /* save updated values */
                    BBE_SANX16_IP(R0, al_R_st0, pR_st0);
                    BBE_SANX16_IP(R1, al_R_st1, pR_st1);
                }

                /* Process last 4 values */
                /* load and rotate values */
                BBE_LAVNX16_XP(R0, al_R_ld0, pR_ld0, sz_i16*4*2);
                BBE_LAVNX16_XP(R1, al_R_ld1, pR_ld1, sz_i16*4*2);
                ACC0 = BBE_MULRNX16J(R0, F0, rnd);
                ACC1 = BBE_MULRNX16J(R1, F1, rnd);
                R0 = BBE_PACKVNX40(ACC0, rnd);
                R1 = BBE_PACKVNX40(ACC1, rnd);
                /* save updated values */
                BBE_SAVNX16_XP(R0, al_R_st0, pR_st0, sz_i16*4*2);
                BBE_SAVNX16_XP(R1, al_R_st1, pR_st1, sz_i16*4*2);
                BBE_SANX16POS_FP(al_R_st0, pR_st0);
                BBE_SANX16POS_FP(al_R_st1, pR_st1);

                /* jump to the next 2 rows of the matrix R */
                pR_ld0 = pR_ld1;
                pR_ld1 = (const xb_vecNx16 *)((int16_t *)pR_ld0 + 2*N);
            }
        }
    }
    else
    {
        __Pragma("loop_count min=1")
        for (l = 0; l < L; l++)
        {
            pFi = (const xb_vecNx16 *)(Fi + 2*l);

            pR_st0 = (xb_vecNx16 *)(R + SA*l);
            pR_st1 = (xb_vecNx16 *)((int16_t *)pR_st0 + 2*N);
            pR_st2 = (xb_vecNx16 *)((int16_t *)pR_st1 + 2*N);
            pR_st3 = (xb_vecNx16 *)((int16_t *)pR_st2 + 2*N);
            pR_ld0 = pR_st0;
            pR_ld1 = pR_st1;
            pR_ld2 = pR_st2;
            pR_ld3 = pR_st3;

            /* Process matrix R by 4 rows per iteration */
            __Pragma("loop_count min=1")
            for (m = 0; m < (N >> 2); m++)
            {
                /* load rotation vectors for 4 rows */
                BBE_LPNX16_XP(F0, pFi, sz_i16*2*L);
                BBE_LPNX16_XP(F1, pFi, sz_i16*2*L);
                BBE_LPNX16_XP(F2, pFi, sz_i16*2*L);
                BBE_LPNX16_XP(F3, pFi, sz_i16*2*L);
                F0 = BBE_REPNX16C(F0, 0);
                F1 = BBE_REPNX16C(F1, 0);
                F2 = BBE_REPNX16C(F2, 0);
                F3 = BBE_REPNX16C(F3, 0);

                /* Innermost loop: process 4 rows by SIMD_WIDTH/2 values per iteration */
                __Pragma("loop_count min=1")
                for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH-1)); n++)
                {
                    /* load and rotate values */
                    BBE_LVNX16_IP(R0, pR_ld0, sz_i16*BBE_SIMD_WIDTH);
                    BBE_LVNX16_IP(R1, pR_ld1, sz_i16*BBE_SIMD_WIDTH);
                    BBE_LVNX16_IP(R2, pR_ld2, sz_i16*BBE_SIMD_WIDTH);
                    BBE_LVNX16_IP(R3, pR_ld3, sz_i16*BBE_SIMD_WIDTH);

                    ACC0 = BBE_MULRNX16J(R0, F0, rnd);
                    ACC1 = BBE_MULRNX16J(R1, F1, rnd);
                    ACC2 = BBE_MULRNX16J(R2, F2, rnd);
                    ACC3 = BBE_MULRNX16J(R3, F3, rnd);

                    R0 = BBE_PACKVNX40(ACC0, rnd);
                    R1 = BBE_PACKVNX40(ACC1, rnd);
                    R2 = BBE_PACKVNX40(ACC2, rnd);
                    R3 = BBE_PACKVNX40(ACC3, rnd);
                    /* save updated values */
                    BBE_SVNX16_IP(R0, pR_st0, sz_i16*BBE_SIMD_WIDTH);
                    BBE_SVNX16_IP(R1, pR_st1, sz_i16*BBE_SIMD_WIDTH);
                    BBE_SVNX16_IP(R2, pR_st2, sz_i16*BBE_SIMD_WIDTH);
                    BBE_SVNX16_IP(R3, pR_st3, sz_i16*BBE_SIMD_WIDTH);
                }
                /* jump to the next 4 rows of the matrix R */
                pR_st0 = pR_st3;
                pR_st1 = (xb_vecNx16 *)((int16_t *)pR_st0 + 2*N);
                pR_st2 = (xb_vecNx16 *)((int16_t *)pR_st1 + 2*N);
                pR_st3 = (xb_vecNx16 *)((int16_t *)pR_st2 + 2*N);

                pR_ld0 = pR_st0;
                pR_ld1 = pR_st1;
                pR_ld2 = pR_st2;
                pR_ld3 = pR_st3;
            }
        }
    }
}

#endif

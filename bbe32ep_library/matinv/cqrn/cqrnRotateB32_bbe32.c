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
//    N==32
void cqrnRotateB32(void* pScr, int16_t* B,const int16_t* Fi,int L)
{
    #define N 32
    const xb_vecNx16 * restrict pf0;
    const xb_vecNx16 * restrict pf1;
    const xb_vecNx16 * restrict pf2;
    const xb_vecNx16 * restrict pf3;
    const xb_vecNx16 * restrict pf4;
    const xb_vecNx16 * restrict pf5;
    const xb_vecNx16 * restrict pf6;
    const xb_vecNx16 * restrict pf7;
    const xb_vecNx16 * restrict pscr_ld;
          xb_vecNx16 * restrict pscr_st;
    const xb_vecNx16 * restrict pB_ld;
          xb_vecNx16 * restrict pB_st;
    xb_vecNx16 B0, B1, B2, B3;
    xb_vecNx16 F0, F1, F2, F3, F4, F5, F6, F7;
    xb_vecNx40 ACC0, ACC1, ACC2, ACC3;
    valign al_f0, al_f1, al_f2, al_f3;
    vsaN rnd14;

    int m, l, L_, Lmod;
    
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0);

    rnd14 = BBE_MOVVSA32(14);
    Lmod = L&(BBE_SIMD_WIDTH/2-1);
    L_ = L-Lmod;
    
    /* Rotate all matrices B except the last L%8 */
    if (L_)
    {
        for (m=0; m<N; m+=8)
        {
            pB_st = (xb_vecNx16 *)(B+2*m);
            pB_ld = pB_st;
            pf0 = (const xb_vecNx16 *)(Fi+m*2*L);
            pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L);
            pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L);
            pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L);
            pf4 = (const xb_vecNx16 *)((int16_t *)pf3+2*L);
            pf5 = (const xb_vecNx16 *)((int16_t *)pf4+2*L);
            pf6 = (const xb_vecNx16 *)((int16_t *)pf5+2*L);
            pf7 = (const xb_vecNx16 *)((int16_t *)pf6+2*L);
            
            /* Innermost loop: process by 8 rows x 8 matrices */
            __Pragma("loop_count min=1")
            for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
            {
                /* load and transpose vectors Fi */
                al_f1 = BBE_LANX16_PP(pf1);
                al_f2 = BBE_LANX16_PP(pf2);
                al_f3 = BBE_LANX16_PP(pf3);
                BBE_LVNX16_IP(F0,        pf0, sz_i16*BBE_SIMD_WIDTH);
                BBE_LANX16_IP(F1, al_f1, pf1);
                BBE_LANX16_IP(F2, al_f2, pf2);
                BBE_LANX16_IP(F3, al_f3, pf3);
                al_f0 = BBE_LANX16_PP(pf4);
                al_f1 = BBE_LANX16_PP(pf5);
                al_f2 = BBE_LANX16_PP(pf6);
                al_f3 = BBE_LANX16_PP(pf7);
                BBE_LANX16_IP(F4, al_f0, pf4);
                BBE_LANX16_IP(F5, al_f1, pf5);
                BBE_LANX16_IP(F6, al_f2, pf6);
                BBE_LANX16_IP(F7, al_f3, pf7);
            
                BBE_DSELNX16I(F1, F0, F1, F0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F3, F2, F3, F2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F5, F4, F5, F4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F7, F6, F7, F6, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F2, F0, F2, F0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F3, F1, F3, F1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F6, F4, F6, F4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F7, F5, F7, F5, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F4, F0, F4, F0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F5, F1, F5, F1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F6, F2, F6, F2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(F7, F3, F7, F3, BBE_DSELI_DEINTERLEAVE_2);
                /* load and rotate matrices B */
                BBE_LVNX16_IP(B0, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B1, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B2, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B3, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
                ACC1 = BBE_MULRNX16J(B1, F1, rnd14);
                ACC2 = BBE_MULRNX16J(B2, F2, rnd14);
                ACC3 = BBE_MULRNX16J(B3, F3, rnd14);

                B0 = BBE_PACKVNX40(ACC0, rnd14);
                B1 = BBE_PACKVNX40(ACC1, rnd14);
                B2 = BBE_PACKVNX40(ACC2, rnd14);
                B3 = BBE_PACKVNX40(ACC3, rnd14);
                /* save updated values */
                BBE_SVNX16_IP(B0, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B1, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B2, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B3, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                
                BBE_LVNX16_IP(B0, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B1, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B2, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B3, pB_ld, (N/8)*sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16J(B0, F4, rnd14);
                ACC1 = BBE_MULRNX16J(B1, F5, rnd14);
                ACC2 = BBE_MULRNX16J(B2, F6, rnd14);
                ACC3 = BBE_MULRNX16J(B3, F7, rnd14);

                B0 = BBE_PACKVNX40(ACC0, rnd14);
                B1 = BBE_PACKVNX40(ACC1, rnd14);
                B2 = BBE_PACKVNX40(ACC2, rnd14);
                B3 = BBE_PACKVNX40(ACC3, rnd14);
                /* save updated values */
                BBE_SVNX16_IP(B0, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B1, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B2, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(B3, pB_st, (N/8)*sz_i16*BBE_SIMD_WIDTH);
            }
        }
    }

    if (Lmod)
    {
        /* transpose last vectors Fi and save it to the scratch */
        pscr_st = (xb_vecNx16 *)((int16_t *)pScr);
        
        pf0 = (const xb_vecNx16 *)(           Fi +2*L_);
        pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L );
        pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L );
        pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L );

        for (m=0; m<(N/8); m++)
        {
            al_f0 = BBE_LANX16_PP(pf0);
            al_f1 = BBE_LANX16_PP(pf1);
            al_f2 = BBE_LANX16_PP(pf2);
            al_f3 = BBE_LANX16_PP(pf3);
            BBE_LAVNX16_XP(F0, al_f0, pf0, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F1, al_f1, pf1, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F2, al_f2, pf2, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F3, al_f3, pf3, Lmod*sz_i16*2);
            pf0 = (const xb_vecNx16 *)((int16_t *)pf3+2*L_);
            pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L );
            pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L );
            pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L );
        
            al_f0 = BBE_LANX16_PP(pf0);
            al_f1 = BBE_LANX16_PP(pf1);
            al_f2 = BBE_LANX16_PP(pf2);
            al_f3 = BBE_LANX16_PP(pf3);
            BBE_LAVNX16_XP(F4, al_f0, pf0, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F5, al_f1, pf1, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F6, al_f2, pf2, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F7, al_f3, pf3, Lmod*sz_i16*2);
            pf0 = (const xb_vecNx16 *)((int16_t *)pf3+2*L_);
            pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L );
            pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L );
            pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L );

            BBE_DSELNX16I(F1, F0, F1, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F3, F2, F3, F2, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F5, F4, F5, F4, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F7, F6, F7, F6, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F2, F0, F2, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F3, F1, F3, F1, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F6, F4, F6, F4, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F7, F5, F7, F5, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F4, F0, F4, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F5, F1, F5, F1, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F6, F2, F6, F2, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F7, F3, F7, F3, BBE_DSELI_DEINTERLEAVE_2);

            BBE_SVNX16_IP(F0, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F1, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F2, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F3, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F4, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F5, pscr_st, 4*sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(F6, pscr_st, (1-24)*(int)sz_i16*BBE_SIMD_WIDTH);
        }
        
        /* Rotate last L%8 matrices B */
        pB_st = (xb_vecNx16 *)(B+2*N*L_);
        pB_ld = pB_st;
        pscr_ld = (const xb_vecNx16 *)pScr;

        __Pragma("no_unroll")
        for (l=0; l<Lmod; l++)
        {
            BBE_LVNX16_IP(F0, pscr_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(F1, pscr_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(F2, pscr_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(F3, pscr_ld, sz_i16*BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(B0, pB_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(B1, pB_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(B2, pB_ld, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(B3, pB_ld, sz_i16*BBE_SIMD_WIDTH);

            ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F1, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F2, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F3, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);

            BBE_SVNX16_IP(B0, pB_st, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(B1, pB_st, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(B2, pB_st, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(B3, pB_st, sz_i16*BBE_SIMD_WIDTH);
        }
    }
    #undef N
}
#endif

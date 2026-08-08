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

#define sz_i16 ((int)sizeof(int16_t))

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
void cqrnRotateB1(void* pScr, int16_t* B,const int16_t* Fi,int N,int SB,int L)
{
    const int16_t    * restrict pf_ref;
    const xb_vecNx16 * restrict pf0;
    const xb_vecNx16 * restrict pf1;
    const xb_vecNx16 * restrict pf2;
    const xb_vecNx16 * restrict pf3;
    const xb_vecNx16 * restrict pf4;
    const xb_vecNx16 * restrict pf5;
    const xb_vecNx16 * restrict pf6;
    const xb_vecNx16 * restrict pf7;
          xb_vecNx16 * restrict pscr0;
          xb_vecNx16 * restrict pscr1;
    const xb_vecNx16 * restrict pB_ld;
          xb_vecNx16 * restrict pB_st;

    xb_vecNx16 B0, B1, B2, B3;
    xb_vecNx16 F0, F1, F2, F3, F4, F5, F6, F7;
    xb_vecNx40 ACC0, ACC1, ACC2, ACC3;
    valign al_f0, al_f1, al_f2, al_f3;
    valign al_Bld, al_Bst;
    vsaN rnd14;

    int m, l, L_, Lmod;

    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0);
    NASSERT(N>0 && N%4==0);

    rnd14 = BBE_MOVVSA32(14);
    Lmod = L&(BBE_SIMD_WIDTH/2-1);
    L_ = L-Lmod;
    
    if (SB==8/*M==N==4*/)
    {
        pB_ld = (const xb_vecNx16 *)(B);
        pB_st = (      xb_vecNx16 *)(pB_ld);
        al_Bld = BBE_LANX16_PP(pB_ld);
        al_Bst = BBE_ZALIGN();
        pf0 = (const xb_vecNx16 *)(Fi);
        pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L );
        pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L );
        pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L );
        al_f0 = BBE_LANX16_PP(pf0);
        al_f1 = BBE_LANX16_PP(pf1);
        al_f2 = BBE_LANX16_PP(pf2);
        al_f3 = BBE_LANX16_PP(pf3);
        for (l = 0; l < (L+BBE_SIMD_WIDTH/2-1)/(BBE_SIMD_WIDTH/2); l++)
        {
            BBE_LAVNX16_XP(F0, al_f0, pf0, L*2*sz_i16-l*(BBE_SIMD_WIDTH/2)*sz_i16);
            BBE_LAVNX16_XP(F1, al_f1, pf1, L*2*sz_i16-l*(BBE_SIMD_WIDTH/2)*sz_i16);
            BBE_LAVNX16_XP(F2, al_f2, pf2, L*2*sz_i16-l*(BBE_SIMD_WIDTH/2)*sz_i16);
            BBE_LAVNX16_XP(F3, al_f3, pf3, L*2*sz_i16-l*(BBE_SIMD_WIDTH/2)*sz_i16);

            BBE_DSELNX16I(F1, F0, F1, F0, BBE_DSELI_INTERLEAVE_2);
            BBE_DSELNX16I(F3, F2, F3, F2, BBE_DSELI_INTERLEAVE_2);
            BBE_DSELNX16I(F2, F0, F2, F0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELNX16I(F3, F1, F3, F1, BBE_DSELI_INTERLEAVE_4);
            
            BBE_LAVNX16_XP(B0, al_Bld, pB_ld, (SB*L-0*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_LAVNX16_XP(B1, al_Bld, pB_ld, (SB*L-1*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_LAVNX16_XP(B2, al_Bld, pB_ld, (SB*L-2*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_LAVNX16_XP(B3, al_Bld, pB_ld, (SB*L-3*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);

            ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F2, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F1, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F3, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);
            
            BBE_SAVNX16_XP(B0, al_Bst, pB_st, (SB*L-0*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_SAVNX16_XP(B1, al_Bst, pB_st, (SB*L-1*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_SAVNX16_XP(B2, al_Bst, pB_st, (SB*L-2*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
            BBE_SAVNX16_XP(B3, al_Bst, pB_st, (SB*L-3*BBE_SIMD_WIDTH)*sz_i16-l*4*BBE_SIMD_WIDTH*sz_i16);
        }
        BBE_SANX16POS_FP(al_Bst, pB_st);
        return;
    }
    /* Rotate all matrices B except the last L%8 */
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        pf_ref = Fi+l*BBE_SIMD_WIDTH;
        pB_ld  = (const xb_vecNx16 *)(B+SB*l*(BBE_SIMD_WIDTH/2));
        pB_st  = (      xb_vecNx16 *)(pB_ld);
        al_Bst = BBE_ZALIGN();

        /* Innermost loop: process by 8 rows x 8 matrices */
        for (m=0; m<(N>>3); m++)
        {
            pf0 = (const xb_vecNx16 *)(pf_ref);
            pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L);
            pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L);
            pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L);
            pf4 = (const xb_vecNx16 *)((int16_t *)pf3+2*L);
            pf5 = (const xb_vecNx16 *)((int16_t *)pf4+2*L);
            pf6 = (const xb_vecNx16 *)((int16_t *)pf5+2*L);
            pf7 = (const xb_vecNx16 *)((int16_t *)pf6+2*L);
            pf_ref = (const int16_t *)pf7+2*L;
            /* load and transpose vectors Fi */
            al_f0 = BBE_LANX16_PP(pf0);
            al_f1 = BBE_LANX16_PP(pf1);
            al_f2 = BBE_LANX16_PP(pf2);
            al_f3 = BBE_LANX16_PP(pf3);
            BBE_LANX16_IP(F0, al_f0, pf0);
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
            /* load and rotate matrices B (first 4 of 8) */
            BBE_LVNX16_XP(B0, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B1, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B2, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B3, pB_ld, SB*sz_i16);

            ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F1, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F2, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F3, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);
            /* save updated values */
            BBE_SVNX16_XP(B0, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B1, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B2, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B3, pB_st, SB*sz_i16);
            
            /* load and rotate matrices B (next 4 of 8) */
            BBE_LVNX16_XP(B0, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B1, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B2, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B3, pB_ld, BBE_SIMD_WIDTH*sz_i16-7*SB*sz_i16);

            ACC0 = BBE_MULRNX16J(B0, F4, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F5, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F6, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F7, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);
            /* save updated values */
            BBE_SVNX16_XP(B0, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B1, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B2, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B3, pB_st, BBE_SIMD_WIDTH*sz_i16-7*SB*sz_i16);
        }
        if (N&4)
        {
            pf0 = (const xb_vecNx16 *)(pf_ref);
            pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L);
            pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L);
            pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L);
            /* load and transpose vectors Fi */
            al_f0 = BBE_LANX16_PP(pf0);
            al_f1 = BBE_LANX16_PP(pf1);
            al_f2 = BBE_LANX16_PP(pf2);
            al_f3 = BBE_LANX16_PP(pf3);
            BBE_LANX16_IP(F0, al_f0, pf0);
            BBE_LANX16_IP(F1, al_f1, pf1);
            BBE_LANX16_IP(F2, al_f2, pf2);
            BBE_LANX16_IP(F3, al_f3, pf3);
            F4 = F5 = F6 = F7 = BBE_ZERONX16();
        
            BBE_DSELNX16I(F1, F0, F1, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F3, F2, F3, F2, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F2, F0, F2, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F3, F1, F3, F1, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F4, F0, F4, F0, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F5, F1, F5, F1, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F6, F2, F6, F2, BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(F7, F3, F7, F3, BBE_DSELI_DEINTERLEAVE_2);
            /* load and rotate matrices B (first 4 of 8) */
            BBE_LVNX16_XP(B0, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B1, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B2, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B3, pB_ld, SB*sz_i16);

            ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F1, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F2, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F3, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);
            /* save updated values */
            BBE_SVNX16_XP(B0, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B1, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B2, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B3, pB_st, SB*sz_i16);
            
            /* load and rotate matrices B (next 4 of 8) */
            BBE_LVNX16_XP(B0, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B1, pB_ld, SB*sz_i16);
            BBE_LVNX16_XP(B2, pB_ld, SB*sz_i16);
            BBE_LVNX16_IP(B3, pB_ld, 0);

            ACC0 = BBE_MULRNX16J(B0, F4, rnd14);
            ACC1 = BBE_MULRNX16J(B1, F5, rnd14);
            ACC2 = BBE_MULRNX16J(B2, F6, rnd14);
            ACC3 = BBE_MULRNX16J(B3, F7, rnd14);

            B0 = BBE_PACKVNX40(ACC0, rnd14);
            B1 = BBE_PACKVNX40(ACC1, rnd14);
            B2 = BBE_PACKVNX40(ACC2, rnd14);
            B3 = BBE_PACKVNX40(ACC3, rnd14);
            /* save updated values */
            BBE_SVNX16_XP(B0, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B1, pB_st, SB*sz_i16);
            BBE_SVNX16_XP(B2, pB_st, SB*sz_i16);
            BBE_SVNX16_IP(B3, pB_st, 0);
        }
    }

    /* Process last L%8 matrices */
    if (Lmod)
    {
        /* transpose last vectors Fi and save it to the scratch */
        pscr0 = (xb_vecNx16 *)((int16_t *)pScr);
        
        pf0 = (const xb_vecNx16 *)(           Fi +2*L_);
        pf1 = (const xb_vecNx16 *)((int16_t *)pf0+2*L );
        pf2 = (const xb_vecNx16 *)((int16_t *)pf1+2*L );
        pf3 = (const xb_vecNx16 *)((int16_t *)pf2+2*L );
        
        //__Pragma("loop_count min=1")
        for (m=0; m<(N>>3); m++)
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

            BBE_SVNX16_IP(F0, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F1, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F2, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F3, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F4, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F5, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F6, pscr0, sz_i16*BBE_SIMD_WIDTH);
        }
        if (N&4)
        {
            al_f0 = BBE_LANX16_PP(pf0);
            al_f1 = BBE_LANX16_PP(pf1);
            al_f2 = BBE_LANX16_PP(pf2);
            al_f3 = BBE_LANX16_PP(pf3);
            BBE_LAVNX16_XP(F0, al_f0, pf0, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F1, al_f1, pf1, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F2, al_f2, pf2, Lmod*sz_i16*2);
            BBE_LAVNX16_XP(F3, al_f3, pf3, Lmod*sz_i16*2);

            BBE_DSELNX16I(F1, F0, F1, F0, BBE_DSELI_INTERLEAVE_2);
            BBE_DSELNX16I(F3, F2, F3, F2, BBE_DSELI_INTERLEAVE_2);
            BBE_DSELNX16I(F2, F0, F2, F0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELNX16I(F3, F1, F3, F1, BBE_DSELI_INTERLEAVE_4);

            BBE_SVNX16_IP(F0, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F2, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F1, pscr0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(F3, pscr0, sz_i16*BBE_SIMD_WIDTH);
        }

        /* Rotate last L%8 matrices B */
        for (l=0; l<Lmod; l++)
        {
            pscr0 = (xb_vecNx16 *)((int16_t*)pScr+l*BBE_SIMD_WIDTH);
            pscr1 = (xb_vecNx16 *)((int16_t*)pScr+(N>>3)*7*BBE_SIMD_WIDTH+l*4*2);
            pB_st = (xb_vecNx16 *)(B+SB*(l+L_));
            pB_ld = pB_st;
            
            for (m=0; m<(N>>3); m++)
            {
                BBE_LVNX16_IP(F0, pscr0, 7*sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(B0, pB_ld, sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
                B0 = BBE_PACKVNX40(ACC0, rnd14);

                BBE_SVNX16_IP(B0, pB_st, sz_i16*BBE_SIMD_WIDTH);
            }
            if (N&4)
            {
                al_f1 = BBE_LANX16_PP(pscr1);
                BBE_LAVNX16_XP(F0, al_f1, pscr1, sz_i16*4*2);
                BBE_LVNX16_IP(B0, pB_ld, sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16J(B0, F0, rnd14);
                B0 = BBE_PACKVNX40(ACC0, rnd14);

                BBE_SVNX16_IP(B0, pB_st, sz_i16*BBE_SIMD_WIDTH);
            }
        }
    }
}
#endif

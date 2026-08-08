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
    BBE32 code for QR decomposition (build_r part), block format
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
//    special case: M==1
void cqrnUpdateR1(int16_t* restrict R,
            const int16_t* restrict v,
                  int SA, int L)
{
    const xb_vecNx16 * restrict pV;
    const xb_vecNx16 * restrict pRld;
          xb_vecNx16 * restrict pRst;
    xb_vecNx16 V0, R0, R1, R2, R3, R4, R5, R6, R7, Z0, _1Q14;
    xb_vecNx40 ACC0;
    valign al_V;
    int l, Lmod;
    vsaN rnd14;
    
    rnd14 = BBE_MOVVSA32(14);
    _1Q14 = BBE_MOVVINT16(4);
    _1Q14 = BBE_SLLINX16(_1Q14, 12);
    Lmod = L & (BBE_SIMD_WIDTH/2-1);

    pV   = (const xb_vecNx16 *)v;
    pRld = (const xb_vecNx16 *)R;
    pRst = (      xb_vecNx16 *)R;
    al_V = BBE_LANX16_PP(pV);

    /* Update by SIMD_WIDTH/2 == 8 matrices per iteration */
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LANX16_IP(V0, al_V, pV);
        /* load 8 values and put it to the one register */
        BBE_LPNX16_XP(R0, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R1, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R2, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R3, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R4, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R5, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R6, pRld, SA*sz_i16);
        BBE_LPNX16_XP(R7, pRld, SA*sz_i16);
        R0 = BBE_SELNX16I(R1, R0, BBE_SELI_INTERLEAVE_2_LO);
        R2 = BBE_SELNX16I(R3, R2, BBE_SELI_INTERLEAVE_2_LO);
        R0 = BBE_SELNX16I(R2, R0, BBE_SELI_INTERLEAVE_4_LO);
        R4 = BBE_SELNX16I(R5, R4, BBE_SELI_INTERLEAVE_2_LO);
        R6 = BBE_SELNX16I(R7, R6, BBE_SELI_INTERLEAVE_2_LO);
        R4 = BBE_SELNX16I(R6, R4, BBE_SELI_INTERLEAVE_4_LO);
        R0 = BBE_SELNX16I(R4, R0, BBE_SELI_EXTRACT_LO_HALVES);

        /* R = R - v*v'*R */
        ACC0 = BBE_MULRNX16J(R0, V0, rnd14);
        Z0 = BBE_PACKVNX40(ACC0, rnd14);
        
        ACC0 = BBE_MULRNX16(R0, _1Q14, rnd14);
        BBE_MULSNX16C(ACC0, Z0, V0);
        R0 = BBE_PACKVNX40(ACC0, rnd14);

        /* save updated values */
        R1 = BBE_REPNX16C(R0, 1);
        R2 = BBE_REPNX16C(R0, 2);
        R3 = BBE_REPNX16C(R0, 3);
        R4 = BBE_REPNX16C(R0, 4);
        R5 = BBE_REPNX16C(R0, 5);
        R6 = BBE_REPNX16C(R0, 6);
        R7 = BBE_REPNX16C(R0, 7);
        BBE_SPNX16_XP(R0, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R1, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R2, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R3, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R4, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R5, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R6, pRst, SA*sz_i16);
        BBE_SPNX16_XP(R7, pRst, SA*sz_i16);
    }

    /* Update last L%8 matrices */
    __Pragma("no_unroll")
    for (l=0; l<Lmod; l++)
    {
        BBE_LPNX16_IP(V0, pV  , 2*sz_i16);
        BBE_LPNX16_XP(R0, pRld, SA*sz_i16);

        ACC0 = BBE_MULRNX16J(R0, V0, rnd14);
        Z0 = BBE_PACKVNX40(ACC0, rnd14);
        
        ACC0 = BBE_MULRNX16(R0, _1Q14, rnd14);
        BBE_MULSNX16C(ACC0, Z0, V0);
        R0 = BBE_PACKVNX40(ACC0, rnd14);
        
        BBE_SPNX16_XP(R0, pRst, SA*sz_i16);
    }
}
#endif

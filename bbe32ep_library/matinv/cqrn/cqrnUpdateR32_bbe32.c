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
    Temporary:
    Z[2*N*L]
-------------------------------------------------------*/
void cqrnUpdateR32(int16_t* restrict Z,
                   int16_t* restrict R,
             const int16_t* restrict v,
                   int SA, int M, int N, int N0, int L)
{
          int16_t    *          R_aligned;
    const xb_vecNx16 * restrict pV;
    const xb_vecNx16 * restrict pRld;
          xb_vecNx16 * restrict pRst;
    xb_vecNx16 V0, R0, R1, R2, R3, Z0, Z1, Z2, Z3, _1Q14;
    xb_vecNx40 ACC0, ACC1, ACC2, ACC3;
    int l, m, n;
    vsaN rnd14;
    
    NASSERT_ALIGN(Z,2*BBE_SIMD_WIDTH);
    NASSERT(N>16 && N<=32);
    NASSERT(M<=32);
    NASSERT(L>0);

    rnd14 = BBE_MOVVSA32(14);
    _1Q14 = BBE_MOVVINT16(4);
    _1Q14 = BBE_SLLINX16(_1Q14, 12);

    /* offset address of matrix R to make loads/stores aligned */
    R_aligned = (int16_t *)((uintptr_t)R & ~(sz_i16*BBE_SIMD_WIDTH-1));

    /* Update L matrices as */
    /* R = R - v*v'*R       */

    if (M==1)
    {
        pV = (const xb_vecNx16 *)(v);
        __Pragma("loop_count min=1");
        for (l=0; l<L; l++)
        {
            BBE_LPNX16_IP(V0, pV, sz_i16*2);
            V0 = BBE_REPNX16C(V0, 0);

            pRld = (const xb_vecNx16 *)(R_aligned+l*SA);
            pRst = (      xb_vecNx16 *)(pRld);
                
            __Pragma("loop_count min=1");
            for (n=0; n<((N+BBE_SIMD_WIDTH/2-1)>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
            {
                /* Z=v'*R */
                BBE_LVNX16_IP(R0, pRld, sz_i16*BBE_SIMD_WIDTH);
                ACC0 = BBE_MULRNX16J(R0, V0, rnd14);
                Z0 = BBE_PACKVNX40(ACC0, rnd14);
    
                /* R=R-v*Z */
                ACC0 = BBE_MULRNX16(R0, _1Q14, rnd14);
                BBE_MULSNX16C(ACC0, Z0, V0);
                R0 = BBE_PACKVNX40(ACC0, rnd14);

                BBE_SVNX16_IP(R0, pRst, sz_i16*BBE_SIMD_WIDTH);
            }
        }
        return;
    }

    /* case for M>1 */
    if (((N+7)& ~7) & 8)
    {
        __Pragma("loop_count min=1");
        for (l=0; l<L; l++)
        {
            /* Z=v'*R */
            pV   = (const xb_vecNx16 *)(v+2*M*l);
            pRld = (const xb_vecNx16 *)(R_aligned+l*SA);
            ACC0 = BBE_ZERONX40();
            ACC0 = ACC1 = ACC2 = BBE_RNDADJNX40(ACC0, rnd14);

            __Pragma("loop_count min=2");
            for (m=0; m<M; m++)
            {
                BBE_LPNX16_IP(V0, pV, sz_i16*2);
                V0 = BBE_REPNX16C(V0, 0);
                BBE_LVNX16_IP(R0, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R1, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(R2, pRld, sz_i16*2*N0-2*sz_i16*BBE_SIMD_WIDTH);

                BBE_MULANX16J(ACC0, R0, V0);
                BBE_MULANX16J(ACC1, R1, V0);
                BBE_MULANX16J(ACC2, R2, V0);
            }
            Z0 = BBE_PACKVNX40(ACC0, rnd14);
            Z1 = BBE_PACKVNX40(ACC1, rnd14);
            Z2 = BBE_PACKVNX40(ACC2, rnd14);

            /* R=R-v*Z */
            pV   = (const xb_vecNx16 *)(v+2*M*l);
            pRld = (const xb_vecNx16 *)(R_aligned+l*SA);
            pRst = (      xb_vecNx16 *)(pRld);

            __Pragma("loop_count min=2");
            for (m=0; m<M; m++)
            {
                BBE_LPNX16_IP(V0, pV, sz_i16*2);
                V0 = BBE_REPNX16C(V0, 0);
                BBE_LVNX16_IP(R0, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R1, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(R2, pRld, sz_i16*2*N0-2*sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16(R0, _1Q14, rnd14);
                ACC1 = BBE_MULRNX16(R1, _1Q14, rnd14);
                ACC2 = BBE_MULRNX16(R2, _1Q14, rnd14);
                BBE_MULSNX16C(ACC0, Z0, V0);
                BBE_MULSNX16C(ACC1, Z1, V0);
                BBE_MULSNX16C(ACC2, Z2, V0);
                R0 = BBE_PACKVNX40(ACC0, rnd14);
                R1 = BBE_PACKVNX40(ACC1, rnd14);
                R2 = BBE_PACKVNX40(ACC2, rnd14);

                BBE_SVNX16_IP(R0, pRst, sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(R1, pRst, sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_XP(R2, pRst, sz_i16*2*N0-2*sz_i16*BBE_SIMD_WIDTH);
            }
        }
    }
    else
    {
        __Pragma("loop_count min=1");
        for (l=0; l<L; l++)
        {
            /* Z=v'*R */
            pV   = (const xb_vecNx16 *)(v+2*M*l);
            pRld = (const xb_vecNx16 *)(R_aligned+l*SA);
            ACC0 = BBE_ZERONX40();
            ACC0 = ACC1 = ACC2 = ACC3 = BBE_RNDADJNX40(ACC0, rnd14);

            __Pragma("loop_count min=2");
            for (m=0; m<M; m++)
            {
                BBE_LPNX16_IP(V0, pV, sz_i16*2);
                V0 = BBE_REPNX16C(V0, 0);
                BBE_LVNX16_IP(R0, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R1, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R2, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(R3, pRld, sz_i16*2*N0-3*sz_i16*BBE_SIMD_WIDTH);

                BBE_MULANX16J(ACC0, R0, V0);
                BBE_MULANX16J(ACC1, R1, V0);
                BBE_MULANX16J(ACC2, R2, V0);
                BBE_MULANX16J(ACC3, R3, V0);
            }
            Z0 = BBE_PACKVNX40(ACC0, rnd14);
            Z1 = BBE_PACKVNX40(ACC1, rnd14);
            Z2 = BBE_PACKVNX40(ACC2, rnd14);
            Z3 = BBE_PACKVNX40(ACC3, rnd14);

            /* R=R-v*Z */
            pV   = (const xb_vecNx16 *)(v+2*M*l);
            pRld = (const xb_vecNx16 *)(R_aligned+l*SA);
            pRst = (      xb_vecNx16 *)(pRld);

            __Pragma("loop_count min=2");
            for (m=0; m<M; m++)
            {
                BBE_LPNX16_IP(V0, pV, sz_i16*2);
                V0 = BBE_REPNX16C(V0, 0);
                BBE_LVNX16_IP(R0, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R1, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(R2, pRld, sz_i16*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(R3, pRld, sz_i16*2*N0-3*sz_i16*BBE_SIMD_WIDTH);

                ACC0 = BBE_MULRNX16(R0, _1Q14, rnd14);
                ACC1 = BBE_MULRNX16(R1, _1Q14, rnd14);
                ACC2 = BBE_MULRNX16(R2, _1Q14, rnd14);
                ACC3 = BBE_MULRNX16(R3, _1Q14, rnd14);
                BBE_MULSNX16C(ACC0, Z0, V0);
                BBE_MULSNX16C(ACC1, Z1, V0);
                BBE_MULSNX16C(ACC2, Z2, V0);
                BBE_MULSNX16C(ACC3, Z3, V0);
                R0 = BBE_PACKVNX40(ACC0, rnd14);
                R1 = BBE_PACKVNX40(ACC1, rnd14);
                R2 = BBE_PACKVNX40(ACC2, rnd14);
                R3 = BBE_PACKVNX40(ACC3, rnd14);

                BBE_SVNX16_IP(R0, pRst, sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(R1, pRst, sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(R2, pRst, sz_i16*BBE_SIMD_WIDTH);
                BBE_SVNX16_XP(R3, pRst, sz_i16*2*N0-3*sz_i16*BBE_SIMD_WIDTH);
            }
        }
    }
}
#endif

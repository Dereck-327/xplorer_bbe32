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
void cqrnUpdateR16(int16_t* restrict Z,
                   int16_t* restrict R,
             const int16_t* restrict v,
                   int SA, int M,int N, int N0, int L)
{
          int16_t    *          R_aligned;
          xb_vecNx16 * restrict pZ;
    const xb_vecNx16 * restrict pV;
    const xb_vecNx16 * restrict pRld0;
          xb_vecNx16 * restrict pRst0;
    xb_vecNx16 V0, V1, V2, V3, V4, V5, V6, V7;
    xb_vecNx16 Z0, Z1;
    xb_vecNx16 R00, R10, R20, R30, R40, R50, R60, R70;
    xb_vecNx16 R01, R11, R21, R31, R41, R51, R61, R71;
    xb_vecNx40 ACC0, ACC1;
    xb_vecNx16 T0, _1Q14;
    vselN shfl;
    vsaN rnd14;
    valign al_pV;
    int l, m, offs;
    
    NASSERT_ALIGN(Z,2*BBE_SIMD_WIDTH);
    NASSERT(N>8 && N<=16);
    NASSERT(M<=16);
    NASSERT(L>0);

    rnd14 = BBE_MOVVSA32(14);
    _1Q14 = BBE_MOVVINT16(4);
    _1Q14 = BBE_SLLINX16(_1Q14, 12);

    /* offset address of matrix R to make loads/stores aligned */
    R_aligned = (int16_t *)((uintptr_t)R & ~(sz_i16*BBE_SIMD_WIDTH-1));

    /* Update L matrices as */
    /* R = R - v*v'*R       */

    if (M<=4) /* M=1..4 */
    {
        offs = (4-M)*2;
        T0 = BBE_SEQNX16();
        T0 = BBE_SUBNX16(T0, BBE_MOVVA16(offs));
        shfl = BBE_MOVVSELNX16(T0, 0);

        pRld0 = (const xb_vecNx16 *)(R_aligned-offs*N0);
        pZ    = (      xb_vecNx16 *)(Z);
        pV    = (const xb_vecNx16 *)(v);
        al_pV = BBE_LANX16_PP(pV);
        
        __Pragma("loop_count min=1")
        for (l=0; l<L; l++)
        {
            /* Z=v'*R */
            BBE_LAVNX16_XP(T0, al_pV, pV, sz_i16*2*M);
            T0 = BBE_SHFLNX16(T0, shfl);
            V0 = BBE_REPNX16C(T0, 0);
            V1 = BBE_REPNX16C(T0, 1);
            V2 = BBE_REPNX16C(T0, 2);
            V3 = BBE_REPNX16C(T0, 3);
            
            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R10, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R11, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R20, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R21, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R30, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R31, pRld0, sz_i16*SA-3*sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
        
            ACC0 = BBE_MULRNX16J(R00, V0, rnd14);
            BBE_MULANX16J (ACC0, R10, V1);
            BBE_MULANX16J (ACC0, R20, V2);
            BBE_MULANX16J (ACC0, R30, V3);
            ACC1 = BBE_MULRNX16J(R01, V0, rnd14);
            BBE_MULANX16J (ACC1, R11, V1);
            BBE_MULANX16J (ACC1, R21, V2);
            BBE_MULANX16J (ACC1, R31, V3);

            Z0 = BBE_PACKVNX40(ACC0, rnd14);
            Z1 = BBE_PACKVNX40(ACC1, rnd14);

            BBE_SVNX16_IP(Z0, pZ, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(Z1, pZ, sz_i16*BBE_SIMD_WIDTH);
        }

        pRld0 = (const xb_vecNx16 *)(R_aligned-offs*N0);
        pRst0 = (      xb_vecNx16 *)(pRld0);
        pZ    = (      xb_vecNx16 *)(Z);
        pV    = (const xb_vecNx16 *)(v);
        al_pV = BBE_LANX16_PP(pV);

        __Pragma("loop_count min=1")
        for (l=0; l<L; l++)
        {
            /* R=R-v*Z */
            BBE_LAVNX16_XP(T0, al_pV, pV, sz_i16*2*M);
            T0 = BBE_SHFLNX16(T0, shfl);
            V0 = BBE_REPNX16C(T0, 0);
            V1 = BBE_REPNX16C(T0, 1);
            V2 = BBE_REPNX16C(T0, 2);
            V3 = BBE_REPNX16C(T0, 3);

            BBE_LVNX16_IP(Z0, pZ, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(Z1, pZ, sz_i16*BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R00, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R01, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V0);
            BBE_MULSNX16C(ACC1, Z1, V0);
            R00 = BBE_PACKVNX40(ACC0, rnd14);
            R01 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R00, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R01, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R10, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R11, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R10, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R11, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V1);
            BBE_MULSNX16C(ACC1, Z1, V1);
            R10 = BBE_PACKVNX40(ACC0, rnd14);
            R11 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R10, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R11, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R20, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R21, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R20, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R21, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V2);
            BBE_MULSNX16C(ACC1, Z1, V2);
            R20 = BBE_PACKVNX40(ACC0, rnd14);
            R21 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R20, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R21, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R30, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R31, pRld0, sz_i16*SA-3*sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R30, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R31, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V3);
            BBE_MULSNX16C(ACC1, Z1, V3);
            R30 = BBE_PACKVNX40(ACC0, rnd14);
            R31 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R30, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R31, pRst0, sz_i16*SA-3*sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
        }
        return;
    }
    if (M<=8) /* M=5..8 */
    {
        offs = (8-M)*2;
        T0 = BBE_SEQNX16();
        T0 = BBE_SUBNX16(T0, BBE_MOVVA16(offs));
        shfl = BBE_MOVVSELNX16(T0, 0);

        pRld0 = (const xb_vecNx16 *)(R_aligned-offs*N0);
        pZ    = (      xb_vecNx16 *)(Z);
        pV    = (const xb_vecNx16 *)(v);
        al_pV = BBE_LANX16_PP(pV);
        
        __Pragma("loop_count min=1")
        for (l=0; l<L; l++)
        {
            /* Z=v'*R */
            BBE_LAVNX16_XP(T0, al_pV, pV, sz_i16*2*M);
            T0 = BBE_SHFLNX16(T0, shfl);
            V0 = BBE_REPNX16C(T0, 0);
            V1 = BBE_REPNX16C(T0, 1);
            V2 = BBE_REPNX16C(T0, 2);
            V3 = BBE_REPNX16C(T0, 3);
            V4 = BBE_REPNX16C(T0, 4);
            V5 = BBE_REPNX16C(T0, 5);
            V6 = BBE_REPNX16C(T0, 6);
            V7 = BBE_REPNX16C(T0, 7);
            
            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R10, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R11, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R20, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R21, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R30, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R31, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R40, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R41, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R50, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R51, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R60, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R61, pRld0, sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(R70, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R71, pRld0, sz_i16*SA-7*sz_i16*2*N0-sz_i16*BBE_SIMD_WIDTH);
        
            ACC0 = BBE_MULRNX16J(R00, V0, rnd14);
            BBE_MULANX16J (ACC0, R10, V1);
            BBE_MULANX16J (ACC0, R20, V2);
            BBE_MULANX16J (ACC0, R30, V3);
            BBE_MULANX16J (ACC0, R40, V4);
            BBE_MULANX16J (ACC0, R50, V5);
            BBE_MULANX16J (ACC0, R60, V6);
            BBE_MULANX16J (ACC0, R70, V7);
            ACC1 = BBE_MULRNX16J(R01, V0, rnd14);
            BBE_MULANX16J (ACC1, R11, V1);
            BBE_MULANX16J (ACC1, R21, V2);
            BBE_MULANX16J (ACC1, R31, V3);
            BBE_MULANX16J (ACC1, R41, V4);
            BBE_MULANX16J (ACC1, R51, V5);
            BBE_MULANX16J (ACC1, R61, V6);
            BBE_MULANX16J (ACC1, R71, V7);

            Z0 = BBE_PACKVNX40(ACC0, rnd14);
            Z1 = BBE_PACKVNX40(ACC1, rnd14);

            BBE_SVNX16_IP(Z0, pZ, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(Z1, pZ, sz_i16*BBE_SIMD_WIDTH);
        }

        pRld0 = (const xb_vecNx16 *)(R_aligned-offs*N0);
        pRst0 = (      xb_vecNx16 *)(pRld0);
        pZ    = (      xb_vecNx16 *)(Z);
        pV    = (const xb_vecNx16 *)(v);
        al_pV = BBE_LANX16_PP(pV);

        __Pragma("loop_count min=1")
        for (l=0; l<L; l++)
        {
            /* R=R-v*Z */
            BBE_LAVNX16_XP(T0, al_pV, pV, sz_i16*2*M);
            T0 = BBE_SHFLNX16(T0, shfl);
            V0 = BBE_REPNX16C(T0, 0);
            V1 = BBE_REPNX16C(T0, 1);
            V2 = BBE_REPNX16C(T0, 2);
            V3 = BBE_REPNX16C(T0, 3);
            V4 = BBE_REPNX16C(T0, 4);
            V5 = BBE_REPNX16C(T0, 5);
            V6 = BBE_REPNX16C(T0, 6);
            V7 = BBE_REPNX16C(T0, 7);

            BBE_LVNX16_IP(Z0, pZ, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(Z1, pZ, sz_i16*BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R00, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R01, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V0);
            BBE_MULSNX16C(ACC1, Z1, V0);
            R00 = BBE_PACKVNX40(ACC0, rnd14);
            R01 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R00, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R01, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R10, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R11, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R10, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R11, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V1);
            BBE_MULSNX16C(ACC1, Z1, V1);
            R10 = BBE_PACKVNX40(ACC0, rnd14);
            R11 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R10, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R11, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R20, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R21, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R20, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R21, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V2);
            BBE_MULSNX16C(ACC1, Z1, V2);
            R20 = BBE_PACKVNX40(ACC0, rnd14);
            R21 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R20, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R21, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R30, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R31, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R30, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R31, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V3);
            BBE_MULSNX16C(ACC1, Z1, V3);
            R30 = BBE_PACKVNX40(ACC0, rnd14);
            R31 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R30, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R31, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(R40, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R41, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R40, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R41, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V4);
            BBE_MULSNX16C(ACC1, Z1, V4);
            R40 = BBE_PACKVNX40(ACC0, rnd14);
            R41 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R40, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R41, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R50, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R51, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R50, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R51, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V5);
            BBE_MULSNX16C(ACC1, Z1, V5);
            R50 = BBE_PACKVNX40(ACC0, rnd14);
            R51 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R50, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R51, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R60, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R61, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R60, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R61, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V6);
            BBE_MULSNX16C(ACC1, Z1, V6);
            R60 = BBE_PACKVNX40(ACC0, rnd14);
            R61 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R60, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R61, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            
            BBE_LVNX16_IP(R70, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R71, pRld0, sz_i16*SA-7*sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
            ACC0 = BBE_MULRNX16(R70, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R71, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V7);
            BBE_MULSNX16C(ACC1, Z1, V7);
            R70 = BBE_PACKVNX40(ACC0, rnd14);
            R71 = BBE_PACKVNX40(ACC1, rnd14);
            BBE_SVNX16_IP(R70, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R71, pRst0, sz_i16*SA-7*sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
        }
        return;
    }

    /* M>8 */
    __Pragma("loop_count min=1");
    for (l=0; l<L; l++)
    {
        /* Z=v'*R */
        pRld0 = (const xb_vecNx16 *)(R_aligned+l*SA);
        pV    = (const xb_vecNx16 *)(v+2*M*l);
        ACC0 = BBE_ZERONX40();
        ACC0 = ACC1 = BBE_RNDADJNX40(ACC0, rnd14);
        
        __Pragma("loop_count min=4");
        for (m=0; m<M; m++)
        {
            BBE_LPNX16_IP(V0, pV, sz_i16*2);
            V0 = BBE_REPNX16C(V0, 0);
            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);

            BBE_MULANX16J(ACC0, R00, V0);
            BBE_MULANX16J(ACC1, R01, V0);
        }
        Z0 = BBE_PACKVNX40(ACC0, rnd14);
        Z1 = BBE_PACKVNX40(ACC1, rnd14);
    
        /* R=R-v*Z */
        pV    = (const xb_vecNx16 *)(v+2*M*l);
        pRld0 = (const xb_vecNx16 *)(R_aligned+l*SA);
        pRst0 = (      xb_vecNx16 *)(pRld0);

        __Pragma("loop_count min=4");
        for (m=0; m<M; m++)
        {
            BBE_LPNX16_IP(V0, pV, sz_i16*2);
            V0 = BBE_REPNX16C(V0, 0);
            BBE_LVNX16_IP(R00, pRld0, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(R01, pRld0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);

            ACC0 = BBE_MULRNX16(R00, _1Q14, rnd14);
            ACC1 = BBE_MULRNX16(R01, _1Q14, rnd14);
            BBE_MULSNX16C(ACC0, Z0, V0);
            BBE_MULSNX16C(ACC1, Z1, V0);
            R00 = BBE_PACKVNX40(ACC0, rnd14);
            R01 = BBE_PACKVNX40(ACC1, rnd14);

            BBE_SVNX16_IP(R00, pRst0, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(R01, pRst0, sz_i16*2*N0-1*sz_i16*BBE_SIMD_WIDTH);
        }
    }
}
#endif

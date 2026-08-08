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
    BBE16EP code for QR decomposition (build_r part), block format
    IntegrIT, 2006-2016
*/


#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cqrn_common.h"

#if HAVE_CQRN

#define sz_i16 sizeof(int16_t)

/*-------------------------------------------------------
    take column from sequence of block ordered matrices
    and put it to the linear array
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
/*
    M=1...4
*/
void cqrnTakeColumn4(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L)
{
    const xb_vecNx16 * restrict pA;
          xb_vecNx16 * restrict pX;
    xb_vecNx16 T0, T1, T2, T3;
    vselN offset;
    valign al_X;
    int l, NbM;
    
    NASSERT(M>=1 && M<=4);
    NASSERT(L>0);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);

    NbM = sz_i16*2*M;
    T0 = BBE_SEQNX16();
    T0 = BBE_ADDNX16(T0, 2*(4-M));
    offset = BBE_MOVVSELNX16(T0, 0);
    pA = (const xb_vecNx16 *)(A-2*(4-M)*N);
    pX = (      xb_vecNx16 *)(x);
    al_X = BBE_ZALIGN();
    
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        /* load 4 values and put it to the one register */
        BBE_LPNX16_XP(T0, pA, N*2*sz_i16);
        BBE_LPNX16_XP(T1, pA, N*2*sz_i16);
        BBE_LPNX16_XP(T2, pA, N*2*sz_i16);
        BBE_LPNX16_XP(T3, pA, SA*sz_i16-3*N*2*sz_i16);
        T0 = BBE_SELNX16I(T1, T0, BBE_SELI_INTERLEAVE_2_LO);
        T2 = BBE_SELNX16I(T3, T2, BBE_SELI_INTERLEAVE_2_LO);
        T0 = BBE_SELNX16I(T2, T0, BBE_SELI_INTERLEAVE_4_LO);

        /* save values */
        T0 = BBE_SELNX16(T0, T0, offset);
        BBE_SAVNX16_XP(T0, al_X, pX, NbM);
    }
    BBE_SANX16POS_FP(al_X, pX);
}

#endif

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
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)

    special case: M==1
-------------------------------------------------------*/
void cqrnHousholder1      (int16_t* restrict v,
                           int16_t* restrict Fi,
                           int16_t *restrict D,
                           const int16_t* restrict x, 
                           int SD, int L)
{
    const xb_vecNx16 * restrict pX;
          xb_vecNx16 * restrict pV;
          xb_vecNx16 * restrict pF;
          xb_vecNx16 * restrict pD;
    xb_vecNx16 X0, V0, D0, INVD0, F0;
    xb_vecNx16 U0, U1, U2, U3, U4, U5, U6, U7;
    xb_vecNx16 drsqrt, d2rsqrt, zero, _1re;
    vsaN D0_exp, U0_exp, norm_exp;
    xb_vecNx40 ACC0, ACC1;
    xb_vecNx40 _0x40000000; 
    vboolN mask;
    vsaN rnd, rnd16;
    valign alX, alF, alV;
    int ___0x40000000=0x40000000;
    int L_, Ltail, l;

    Ltail = L & (BBE_SIMD_WIDTH/2-1);
    L_ = L-Ltail;
    rnd16 = BBE_MOVVSA32(16);
    zero = BBE_ZERONX16();
    _1re = BBE_MOVVA16C(0x00007FFF);
    pV = (      xb_vecNx16 *)(v);
    pF = (      xb_vecNx16 *)(Fi);
    pD = (      xb_vecNx16 *)(D);
    pX = (const xb_vecNx16 *)(x);
    alX = BBE_LANX16_PP(pX);
    alF = BBE_ZALIGN();

    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LANX16_IP(X0, alX, pX);

        // compute d=1/sqrt(x'x)
        ACC0 = BBE_MULNX16J(X0, X0);
        ACC1 = BBE_SLSINX40(ACC0,5);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40 (ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        D0 = BBE_PACKLNX40(ACC0);
        D0 = BBE_SHFLNX16I(D0, BBE_SHFLI_DUPLICATE_1_EVEN);
        mask = BBE_NEQNX16(D0, zero);
        D0_exp = BBE_SUBSR1SAVSN(16, norm_exp );
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);
        U0 = BBE_MOVVVS(D0_exp);
        U0 = BBE_SELNX16I(U0, D0, BBE_SELI_INTERLEAVE_1_EVEN);
        U1 = BBE_REPNX16C(U0, 1);
        U2 = BBE_REPNX16C(U0, 2);
        U3 = BBE_REPNX16C(U0, 3);
        U4 = BBE_REPNX16C(U0, 4);
        U5 = BBE_REPNX16C(U0, 5);
        U6 = BBE_REPNX16C(U0, 6);
        U7 = BBE_REPNX16C(U0, 7);
        BBE_SPNX16_XP(U0, pD, SD*sz_i16);
        BBE_SPNX16_XP(U1, pD, SD*sz_i16);
        BBE_SPNX16_XP(U2, pD, SD*sz_i16);
        BBE_SPNX16_XP(U3, pD, SD*sz_i16);
        BBE_SPNX16_XP(U4, pD, SD*sz_i16);
        BBE_SPNX16_XP(U5, pD, SD*sz_i16);
        BBE_SPNX16_XP(U6, pD, SD*sz_i16);
        BBE_SPNX16_XP(U7, pD, SD*sz_i16);

        // compute fi=x0/sqrt(x0'x0)
        rnd = BBE_ADDSAVSN(1, D0_exp);
        ACC0 = BBE_MULUSNX16(D0, X0);
        ACC0 = BBE_NEGNX40(ACC0);
        ACC0 = BBE_RNDADJNX40(ACC0, rnd);
        F0 = BBE_PACKVNX40(ACC0, rnd);
        F0 = BBE_MOVNX16T(F0, _1re, mask );
        BBE_SANX16_IP(F0, alF, pF);
    }
    BBE_SANX16POS_FP(alF, pF);
    __Pragma("no_unroll");
    for (l=0; l<Ltail; l++)
    {
        BBE_LPNX16_IP(X0, pX, 2*sz_i16);

        // compute d=1/sqrt(x'x)
        ACC0 = BBE_MULNX16J(X0, X0);
        ACC1 = BBE_SLSINX40(ACC0,5);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40 (ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        D0 = BBE_PACKLNX40(ACC0);
        D0 = BBE_SHFLNX16I(D0, BBE_SHFLI_DUPLICATE_1_EVEN);
        mask = BBE_NEQNX16(D0, zero);
        D0_exp = BBE_SUBSR1SAVSN(16, norm_exp );
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);
        U0 = BBE_MOVVVS(D0_exp);
        U0 = BBE_SELNX16I(U0, D0, BBE_SELI_INTERLEAVE_1_LO);
        BBE_SPNX16_XP(U0, pD, SD*sz_i16);

        // compute fi=x0/sqrt(x0'x0)
        rnd = BBE_ADDSAVSN(1, D0_exp);
        ACC0 = BBE_MULUSNX16(D0, X0);
        ACC0 = BBE_NEGNX40(ACC0);
        ACC0 = BBE_RNDADJNX40(ACC0, rnd);
        F0 = BBE_PACKVNX40(ACC0, rnd);
        F0 = BBE_MOVNX16T(F0, _1re, mask );
        BBE_SPNX16_IP(F0, pF, 2*sz_i16);
    }



    pX   = (const xb_vecNx16 *)(x);
    pF   = (      xb_vecNx16 *)(Fi);
    pD   = (      xb_vecNx16 *)(D);
    alX = BBE_LANX16_PP(pX);
    alF = BBE_LANX16_PP(pF);
    alV = BBE_ZALIGN();
    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LANX16_IP(X0, alX, pX);
        BBE_LANX16_IP(F0, alF, pF);
        BBE_LPNX16_XP(U0, pD, SD*sz_i16);
        BBE_LPNX16_XP(U1, pD, SD*sz_i16);
        BBE_LPNX16_XP(U2, pD, SD*sz_i16);
        BBE_LPNX16_XP(U3, pD, SD*sz_i16);
        BBE_LPNX16_XP(U4, pD, SD*sz_i16);
        BBE_LPNX16_XP(U5, pD, SD*sz_i16);
        BBE_LPNX16_XP(U6, pD, SD*sz_i16);
        BBE_LPNX16_XP(U7, pD, SD*sz_i16);
        U0 = BBE_SELNX16I(U1, U0, BBE_SELI_INTERLEAVE_2_LO);
        U2 = BBE_SELNX16I(U3, U2, BBE_SELI_INTERLEAVE_2_LO);
        U0 = BBE_SELNX16I(U2, U0, BBE_SELI_INTERLEAVE_4_LO);
        U4 = BBE_SELNX16I(U5, U4, BBE_SELI_INTERLEAVE_2_LO);
        U6 = BBE_SELNX16I(U7, U6, BBE_SELI_INTERLEAVE_2_LO);
        U4 = BBE_SELNX16I(U6, U4, BBE_SELI_INTERLEAVE_4_LO);
        U0 = BBE_SELNX16I(U4, U0, BBE_SELI_EXTRACT_LO_HALVES);
        D0 = U0;
        D0_exp = BBE_MOVVSV(U0, 0);
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_ODD);

        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        // find 1/d
        _0x40000000 = BBE_MOVWA32(___0x40000000);
        BBE_DIVNX32U_4STEP0_0(_0x40000000, D0);
        BBE_DIVNX16U_4STEP_0(D0);
        BBE_DIVNX16U_4STEP_0(D0);
		INVD0 = BBE_DIVNX16U_4STEPN_0(D0);
		INVD0 = BBE_SHFLNX16I(INVD0, BBE_SHFLI_DUPLICATE_1_EVEN);
        // find v=x+fi/d
        rnd = BBE_ADDSAVSN(-14, D0_exp);
        ACC0 = BBE_MULUSRNX16(INVD0, F0, rnd16);
        ACC1 = BBE_UNPKNVNX16(X0, rnd);
        ACC0 = BBE_SUBNX40(ACC1, ACC0);
        V0 = BBE_PACKVNX40(ACC0, rnd16);
        
        // normalization
        // u = 1/sqrt(v'*v)
        ACC0 = BBE_MULNX16J(V0, V0);
        ACC1 = BBE_SLSINX40(ACC0,6);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40(ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SUBSR1SAVSN(17, norm_exp );
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);
        // v = v/sqrt(v'*v)
        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);
        BBE_SANX16_IP(V0, alV, pV);
    }
    BBE_SANX16POS_FP(alV, pV);
    __Pragma("no_unroll");
    for (l=0; l<Ltail; l++)
    {
        BBE_LPNX16_IP(X0, pX, 2*sz_i16);
        BBE_LPNX16_IP(F0, pF, 2*sz_i16);
        BBE_LPNX16_XP(D0, pD, SD*sz_i16);
        D0_exp = BBE_MOVVSV(D0, 0);
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_ODD);

        // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
        // find 1/d
        _0x40000000 = BBE_MOVWA32(___0x40000000);
        BBE_DIVNX32U_4STEP0_0(_0x40000000, D0);
        BBE_DIVNX16U_4STEP_0(D0);
        BBE_DIVNX16U_4STEP_0(D0);
		INVD0 = BBE_DIVNX16U_4STEPN_0(D0);
		INVD0 = BBE_SHFLNX16I(INVD0, BBE_SHFLI_DUPLICATE_1_EVEN);
        // find v=x+fi/d
        rnd = BBE_ADDSAVSN(-14, D0_exp);
        ACC0 = BBE_MULUSRNX16(INVD0, F0, rnd16);
        ACC1 = BBE_UNPKNVNX16(X0, rnd);
        ACC0 = BBE_SUBNX40(ACC1, ACC0);
        V0 = BBE_PACKVNX40(ACC0, rnd16);
        
        // normalization
        // u = 1/sqrt(v'*v)
        ACC0 = BBE_MULNX16J(V0, V0);
        ACC1 = BBE_SLSINX40(ACC0,6);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40(ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SUBSR1SAVSN(17, norm_exp );
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);
        // v = v/sqrt(v'*v)
        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);
        BBE_SPNX16_IP(V0, pV, 2*sz_i16);
    }
}
#endif

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
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)
    temporary:
    pScr[]    scratch, defined by cqrnHousholder_getScratchSz()
-------------------------------------------------------*/
void cqrnHousholder8(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L)
{
    const xb_vecNx16 * restrict pX;
          xb_vecNx16 * restrict pVScr;
          xb_vecNx16 * restrict pV;
          xb_vecNx16 * restrict pU;
          xb_vecNx16 * restrict pF;
          xb_vecNx16 * restrict pD;
          xb_vecNx16 * restrict pX0;  // !!!
          xb_vecNx16 * restrict pXabs;// !!!
          xb_vecNx16 * restrict pdiag;// !!!
    xb_vecNx16 X0, X1;
    xb_vecNx16 V0, D0, INVD0, F0, T0, T1;
    xb_vecNx16 U0, U1, U2, U3, U4, U5, U6, U7;
    xb_vecNx16 drsqrt, d2rsqrt, zero, _1re;
    xb_vecNx16 _0x4000;
    vsaN D0_exp, U0_exp, U1_exp, U2_exp, U3_exp, U4_exp, U5_exp, U6_exp, U7_exp, norm_exp;
    xb_vecNx40 ACC0, ACC1;
    xb_vecNx40 _0x40000000;
    vboolN mask;
    vsaN rnd, rnd16;
    valign al_pX, al_pV, al_pF;
    int ___0x40000000=0x40000000;
    int16_t  *x0;     /* [SL]   */
    int16_t  *diag_v; /* [SL]   */
    int16_t  *xabs_u; /* [SL]*2 */
    int L_, Lmod, l, NbM;

    NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x   , 2*BBE_SIMD_WIDTH);
    NASSERT(M<=8);

    NbM = M*2*sz_i16;
    L_ = L & ~(BBE_SIMD_WIDTH/2-1);
    Lmod = L-L_;

    rnd16 = BBE_MOVVSA32(16);
    zero = BBE_ZERONX16();
    _0x4000= BBE_MOVVA16(0x4000);
    _1re = BBE_MOVVA16C(0x00007FFF);

    /* allocate temporary arrays in scratch */
    {
        int SL;
        SL = ((L+BBE_SIMD_WIDTH/2-1) & ~(BBE_SIMD_WIDTH/2-1))<<1;
        x0      = (int16_t*)pScr;
        diag_v  = x0+SL;
        xabs_u  = diag_v+SL;
    }
    
    /* copy first element of x to x0 */
    pX0   = (      xb_vecNx16 *)(x0);
    pXabs = (      xb_vecNx16 *)(xabs_u);
    pX    = (const xb_vecNx16 *)(x);
    al_pX = BBE_LANX16_PP(pX);
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        BBE_SPNX16_IP(X0, pX0, sz_i16*2);

        ACC0 = BBE_MULNX16J(X0, X0);
        ACC0 = BBE_MOVNX40_FROMC40(BBE_RADDNX40C(ACC0));
        U0 = BBE_MOVSVWL(ACC0);
        BBE_SPNX16_IP(U0, pXabs, 2*sz_i16);
    }

    pXabs = (xb_vecNx16 *)(xabs_u);
    pD    = (xb_vecNx16 *)(D);
    pdiag = (xb_vecNx16 *)(diag_v);
    // compute d=1/sqrt(x'x)
    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(X0, pXabs, sz_i16*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(X1, X0, X0, X0, BBE_DSELI_INTERLEAVE_2);

        ACC0 = BBE_MOVSWV(X1, X0);
        ACC1 = BBE_SLSINX40(ACC0,5);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40 (ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        D0 = BBE_PACKLNX40(ACC0);
        D0_exp = BBE_SUBSR1SAVSN(16, norm_exp );
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
        BBE_SVNX16_IP(U0, pdiag, sz_i16*BBE_SIMD_WIDTH);
    }
    for (l=0; l<Lmod; l++)
    {
        BBE_LPNX16_IP(X0, pXabs, 2*sz_i16);

        ACC0 = BBE_MOVSWVL(X0);
        ACC1 = BBE_SLSINX40(ACC0,5);
        norm_exp = BBE_NSAENX40(ACC1);
        ACC1 = BBE_SLLNX40 (ACC1, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        D0 = BBE_PACKLNX40(ACC0);
        D0_exp = BBE_SUBSR1SAVSN(16, norm_exp );
        U0 = BBE_MOVVVS(D0_exp);
        U0 = BBE_SELNX16I(U0, D0, BBE_SELI_INTERLEAVE_1_LO);
        BBE_SPNX16_XP(U0, pD, SD*sz_i16);
        BBE_SPNX16_IP(U0, pdiag, 2*sz_i16);
    }

    // compute fi=x0/sqrt(x0'x0)
    pF  = (xb_vecNx16 *)(Fi);
    al_pF = BBE_ZALIGN();
    pX0 = (xb_vecNx16 *)(x0);
    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(X0, pX0, sz_i16*BBE_SIMD_WIDTH);

        // compute d=1/sqrt(x0'x0)
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
        D0_exp = BBE_SUBSR1SAVSN(17, norm_exp );
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        // compute fi=x0/sqrt(x0'x0)
        rnd = D0_exp;
        ACC0 = BBE_MULUSNX16(D0, X0);
        ACC0 = BBE_NEGNX40(ACC0);
        ACC0 = BBE_RNDADJNX40(ACC0, rnd);
        F0 = BBE_PACKVNX40(ACC0, rnd);
        F0 = BBE_MOVNX16T(F0, _1re, mask );
        BBE_SANX16_IP(F0, al_pF, pF);
    }
    BBE_SANX16POS_FP(al_pF, pF);
    for (l=0; l<Lmod; l++)
    {
        BBE_LPNX16_IP(X0, pX0, 2*sz_i16);

        // compute d=1/sqrt(x0'x0)
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
        D0_exp = BBE_SUBSR1SAVSN(17, norm_exp );
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        // compute fi=x0/sqrt(x0'x0)
        rnd = D0_exp;
        ACC0 = BBE_MULUSNX16(D0, X0);
        ACC0 = BBE_NEGNX40(ACC0);
        ACC0 = BBE_RNDADJNX40(ACC0, rnd);
        F0 = BBE_PACKVNX40(ACC0, rnd);
        F0 = BBE_MOVNX16T(F0, _1re, mask );
        BBE_SPNX16_IP(F0, pF, 2*sz_i16);
    }

    // compute reciprocal of main diagonal elements
    // update v=x+fi/d
    pdiag = (xb_vecNx16 *)(diag_v);
    pVScr = (xb_vecNx16 *)(diag_v);
    pX0   = (xb_vecNx16 *)(x0);
    pXabs = (xb_vecNx16 *)(xabs_u);
    pU    = (xb_vecNx16 *)(xabs_u);
    pF    = (xb_vecNx16 *)(Fi);
    al_pF = BBE_LANX16_PP(pF);
    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(D0, pdiag, sz_i16*BBE_SIMD_WIDTH);
        D0_exp = BBE_MOVVSV(D0, 0);
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_ODD);
        // find 1/d
        _0x40000000 = BBE_MOVWA32(___0x40000000);
        BBE_DIVNX32U_4STEP0_0(_0x40000000, D0);
        BBE_DIVNX16U_4STEP_0(D0);
        BBE_DIVNX16U_4STEP_0(D0);
		INVD0 = BBE_DIVNX16U_4STEPN_0(D0);
		INVD0 = BBE_SHFLNX16I(INVD0, BBE_SHFLI_DUPLICATE_1_EVEN);
        // find v=x+fi/d
        BBE_LVNX16_IP(X0, pX0, sz_i16*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(F0, al_pF, pF);
        D0_exp = BBE_SUBSAVSN(15, D0_exp);
        ACC1 = BBE_MULNX16(X0, _0x4000);
        ACC1 = BBE_SLSNX40(ACC1, D0_exp);
        ACC1 = BBE_RNDADJNX40(ACC1, rnd16);
        ACC0 = BBE_MULUSNX16(INVD0, F0);
        ACC0 = BBE_SUBNX40(ACC1, ACC0);
        V0 = BBE_PACKVNX40(ACC0, rnd16);
        
        // u = 1/sqrt(v'*v)
        BBE_LVNX16_IP(U0, pXabs, sz_i16*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(U1, U0, U0, U0, BBE_DSELI_INTERLEAVE_2);
        ACC0 = BBE_MOVSWV(U1, U0);

        BBE_MULSNX16J(ACC0, X0, X0);
        D0_exp = BBE_ADDSAVSN(-2, D0_exp);
        norm_exp = BBE_ADDSVSN(D0_exp, D0_exp);
        ACC0 = BBE_SLSNX40(ACC0, norm_exp);
        BBE_MULANX16J(ACC0, V0, V0);

        norm_exp = BBE_NSAENX40(ACC0);
        ACC1 = BBE_SLLNX40(ACC0, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, d2rsqrt, drsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0_exp = BBE_SUBSR1SAVSN(20, norm_exp );
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        // normalization coefficients
        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);
        U0_exp = BBE_SUBSVSN(U0_exp, D0_exp);
        BBE_SVNX16_IP(V0, pVScr, sz_i16*BBE_SIMD_WIDTH);
        T0 = U0;
        T1 = BBE_MOVVVS(U0_exp);
        T0 = BBE_SELNX16I(T1, T0, BBE_SELI_INTERLEAVE_1_EVEN);
        BBE_SVNX16_IP(T0, pU, sz_i16*BBE_SIMD_WIDTH);
    }

    // renormalize Housholder vectors to sum(abs(v).^2)==1
    pVScr = (      xb_vecNx16 *)(diag_v);
    pU    = (      xb_vecNx16 *)(xabs_u);
    pX    = (const xb_vecNx16 *)(x);
    pV    = (      xb_vecNx16 *)(v);
    al_pX = BBE_LANX16_PP(pX);
    al_pV = BBE_ZALIGN();
    T0 = BBE_SEQNX16();
    T0 = BBE_SHFLNX16I(T0, BBE_SHFLI_DOUBLE_1_LO);
    mask = BBE_NEQNX16(T0, zero);
    for (l=0; l<(L_>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(V0, pVScr, sz_i16*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(T0, pU   , sz_i16*BBE_SIMD_WIDTH);

        U0 = BBE_REPNX16(T0, 1);
        U1 = BBE_REPNX16(T0, 3);
        U2 = BBE_REPNX16(T0, 5);
        U3 = BBE_REPNX16(T0, 7);
        U4 = BBE_REPNX16(T0, 9);
        U5 = BBE_REPNX16(T0, 11);
        U6 = BBE_REPNX16(T0, 13);
        U7 = BBE_REPNX16(T0, 15);
        U0_exp = BBE_MOVVSV(U0, 0);
        U1_exp = BBE_MOVVSV(U1, 0);
        U2_exp = BBE_MOVVSV(U2, 0);
        U3_exp = BBE_MOVVSV(U3, 0);
        U4_exp = BBE_MOVVSV(U4, 0);
        U5_exp = BBE_MOVVSV(U5, 0);
        U6_exp = BBE_MOVVSV(U6, 0);
        U7_exp = BBE_MOVVSV(U7, 0);
        U0 = BBE_REPNX16(T0, 0);
        U1 = BBE_REPNX16(T0, 2);
        U2 = BBE_REPNX16(T0, 4);
        U3 = BBE_REPNX16(T0, 6);
        U4 = BBE_REPNX16(T0, 8);
        U5 = BBE_REPNX16(T0, 10);
        U6 = BBE_REPNX16(T0, 12);
        U7 = BBE_REPNX16(T0, 14);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U0, X0, U0_exp);
        X0 = BBE_PACKVNX40(ACC0, U0_exp);
        T1 = V0;
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U1, X0, U1_exp);
        X0 = BBE_PACKVNX40(ACC0, U1_exp);
        T1 = BBE_REPNX16C(V0, 1);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U2, X0, U2_exp);
        X0 = BBE_PACKVNX40(ACC0, U2_exp);
        T1 = BBE_REPNX16C(V0, 2);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U3, X0, U3_exp);
        X0 = BBE_PACKVNX40(ACC0, U3_exp);
        T1 = BBE_REPNX16C(V0, 3);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U4, X0, U4_exp);
        X0 = BBE_PACKVNX40(ACC0, U4_exp);
        T1 = BBE_REPNX16C(V0, 4);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U5, X0, U5_exp);
        X0 = BBE_PACKVNX40(ACC0, U5_exp);
        T1 = BBE_REPNX16C(V0, 5);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U6, X0, U6_exp);
        X0 = BBE_PACKVNX40(ACC0, U6_exp);
        T1 = BBE_REPNX16C(V0, 6);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);

        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U7, X0, U7_exp);
        X0 = BBE_PACKVNX40(ACC0, U7_exp);
        T1 = BBE_REPNX16C(V0, 7);
        X0 = BBE_MOVNX16T(X0, T1, mask);
        BBE_SAVNX16_XP(X0, al_pV, pV, NbM);
    }
    for (l=0; l<Lmod; l++)
    {
        BBE_LPNX16_IP(D0, pdiag, 2*sz_i16);
        D0_exp = BBE_MOVVSV(D0, 0);
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_ODD);
        BBE_LPNX16_IP(F0, pF, 2*sz_i16);
        BBE_LPNX16_IP(X0, pX0, 2*sz_i16);
        // find 1/d
        _0x40000000 = BBE_MOVWA32(___0x40000000);
        BBE_DIVNX32U_4STEP0_0(_0x40000000, D0);
        BBE_DIVNX16U_4STEP_0(D0);
        BBE_DIVNX16U_4STEP_0(D0);
		INVD0 = BBE_DIVNX16U_4STEPN_0(D0);
		INVD0 = BBE_SHFLNX16I(INVD0, BBE_SHFLI_DUPLICATE_1_EVEN);
        // find v=x+fi/d
        D0_exp = BBE_SUBSAVSN(15, D0_exp);
        ACC1 = BBE_MULNX16(X0, _0x4000);
        ACC1 = BBE_SLSNX40(ACC1, D0_exp);
        ACC1 = BBE_RNDADJNX40(ACC1, rnd16);
        ACC0 = BBE_MULUSNX16(INVD0, F0);
        ACC0 = BBE_SUBNX40(ACC1, ACC0);
        V0 = BBE_PACKVNX40(ACC0, rnd16);
        
        // u = 1/sqrt(v'*v)
        BBE_LPNX16_IP(U0, pXabs, 2*sz_i16);
        ACC0 = BBE_MOVSWVL(U0);
        BBE_MULSNX16J(ACC0, X0, X0);
        D0_exp = BBE_ADDSAVSN(-2, D0_exp);
        norm_exp = BBE_ADDSVSN(D0_exp, D0_exp);
        ACC0 = BBE_SLSNX40(ACC0, norm_exp);
        BBE_MULANX16J(ACC0, V0, V0);

        norm_exp = BBE_NSAENX40(ACC0);
        ACC1 = BBE_SLLNX40(ACC0, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, d2rsqrt, drsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0_exp = BBE_SUBSR1SAVSN(20, norm_exp );
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        // normalization
        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);

        U0_exp = BBE_SUBSVSN(U0_exp, D0_exp);
        U0 = BBE_REPNX16C(U0, 0);
        T0 = BBE_MOVVVS(U0_exp);
        T0 = BBE_REPNX16C(T0, 0);
        U0_exp = BBE_MOVVSV(T0, 0);
        BBE_LAVNX16_XP(X0, al_pX, pX, NbM);
        ACC0 = BBE_MULUSRNX16(U0, X0, U0_exp);
        U0 = BBE_PACKVNX40(ACC0, U0_exp);
        U0 = BBE_MOVNX16T(U0, V0, mask);
        BBE_SAVNX16_XP(U0, al_pV, pV, NbM);
    }
    BBE_SAVNX16POS_FP(al_pV, pV);

}
#endif

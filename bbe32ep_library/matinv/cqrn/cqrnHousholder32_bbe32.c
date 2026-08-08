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
void cqrnHousholder32(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L)
{
    const xb_vecNx16 *  pX;
          xb_vecNx16 *  pV;
          xb_vecNx16 *  pF;
          xb_vecNx16 *  pD;
          xb_vecNx16 * pX0;  // !!!
          xb_vecNx16 * pXabs;// !!!
          xb_vecNx16 *  pdiag;// !!!
    xb_vecNx16 X0, X1, X2, X3, Xdiag;
    xb_vecNx16 V0, D0, INVD0, F0, T0;
    xb_vecNx16 U0, U1, U2, U3, U4, U5, U6, U7;
    xb_vecNx16 drsqrt, d2rsqrt, zero, _1re;
    xb_vecNx16 _0x4000;
    vsaN D0_exp, U0_exp, norm_exp;
    xb_vecNx40 ACC0, ACC1, ACC2, ACC3;
    xb_vecNx40 _0x40000000;
    vboolN mask;
    vsaN rnd, rnd16;
    valign al_pX, al_pV, al_pF;
    int ___0x40000000=0x40000000;
    int16_t  *x0;     /* [SL]   */
    int16_t  *diag;   /* [SL]   */
    int16_t  *xabs;   /* [SL]*2 */
    int L_, Lmod, l, NbM0, NbM1;

    NASSERT_ALIGN(pScr, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x   , 2*BBE_SIMD_WIDTH);
    NASSERT(M>16 && M<=32);

    {
        int sz=(int)cqrnHousholder_getScratchSz(M, L);
        xb_vecNx16 * pWR = (xb_vecNx16 *)pScr;
        xb_vecNx16 veczero = BBE_ZERONX16();
        for (l = 0; l < sz; l += 2 * BBE_SIMD_WIDTH)
        {
            BBE_SVNX16_IP(veczero, pWR, 2 * BBE_SIMD_WIDTH);
        }
    }
    NbM0 = M*2*sz_i16-2*sz_i16*BBE_SIMD_WIDTH;
    NbM1 = M*2*sz_i16-3*sz_i16*BBE_SIMD_WIDTH;
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
        diag    = x0+SL;
        xabs    = diag+SL;
    }
    
    /* copy first element of x to x0 */
    pX0   = (      xb_vecNx16 *)(x0);
    pXabs = (      xb_vecNx16 *)(xabs);
    pX    = (const xb_vecNx16 *)(x);
    al_pX = BBE_LANX16_PP(pX);
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        BBE_LANX16_IP (X0, al_pX, pX);
        BBE_LANX16_IP (X1, al_pX, pX);
        BBE_LAVNX16_XP(X2, al_pX, pX, NbM0);
        BBE_LAVNX16_XP(X3, al_pX, pX, NbM1);
        BBE_SPNX16_IP(X0, pX0, sz_i16*2);

        ACC0 = BBE_MULNX16J(X0, X0);
        BBE_MULANX16J(ACC0, X1, X1);
        BBE_MULANX16J(ACC0, X2, X2);
        BBE_MULANX16J(ACC0, X3, X3);
        ACC0 = BBE_MOVNX40_FROMC40(BBE_RADDNX40C(ACC0));
        U0 = BBE_MOVSVWL(ACC0);
        BBE_SPNX16_IP(U0, pXabs, 2*sz_i16);
    }

    pXabs = (xb_vecNx16 *)(xabs);
    pD    = (xb_vecNx16 *)(D);
    pdiag = (xb_vecNx16 *)(diag);
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
    __Pragma("no_unroll")
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
    // update v=x+fi/d and renormalize to sum(abs(v).^2)==1
    pdiag = (xb_vecNx16 *)(diag);
    pX0   = (xb_vecNx16 *)(x0);
    pXabs = (xb_vecNx16 *)(xabs);
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
        X0 = BBE_LVNX16_I(pX0, 0);
        BBE_LANX16_IP(F0, al_pF, pF);
        D0_exp = BBE_SUBSAVSN(15, D0_exp);
        ACC1 = BBE_MULNX16(X0, _0x4000);
        ACC1 = BBE_SLSNX40(ACC1, D0_exp);
        ACC1 = BBE_RNDADJNX40(ACC1, rnd16);
        ACC0 = BBE_MULUSNX16(INVD0, F0);
        ACC0 = BBE_SUBNX40(ACC1, ACC0);
        V0 = BBE_PACKVNX40(ACC0, rnd16);
        
        // u = 1/sqrt(v'*v)
        U0 = BBE_LVNX16_I(pXabs, 0);
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
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0_exp = BBE_SUBSR1SAVSN(20, norm_exp );
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);
        BBE_SVNX16_IP(V0, pX0, sz_i16*BBE_SIMD_WIDTH);

        U0_exp = BBE_SUBSVSN(U0_exp, D0_exp);
        U1 = BBE_MOVVVS(U0_exp);
        T0 = BBE_SELNX16I(U1, U0, BBE_SELI_INTERLEAVE_1_EVEN);
        BBE_SVNX16_IP(T0, pXabs, sz_i16*BBE_SIMD_WIDTH);
    }
    __Pragma("no_unroll")
    for (l=0; l<Lmod; l++)
    {
        BBE_LPNX16_IP(D0, pdiag, 2*sz_i16);
        D0_exp = BBE_MOVVSV(D0, 0);
        D0_exp = BBE_SHFLVSNI(D0_exp, BBE_VSA_SHFLI_DUPLICATE_1_ODD);
        BBE_LPNX16_IP(F0, pF, 2*sz_i16);
        X0 = BBE_LPNX16_I(pX0, 0);
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
        U0 = BBE_LPNX16_I(pXabs, 0);
        ACC0 = BBE_MOVSWVL(U0);
        BBE_MULSNX16J(ACC0, X0, X0);
        D0_exp = BBE_ADDSAVSN(-2, D0_exp);
        norm_exp = BBE_ADDSVSN(D0_exp, D0_exp);
        ACC0 = BBE_SLSNX40(ACC0, norm_exp);
        BBE_MULANX16J(ACC0, V0, V0);

        norm_exp = BBE_NSAENX40(ACC0);
        ACC1 = BBE_SLLNX40(ACC0, norm_exp);
        BBE_RSQRTLUNX40_0(ACC0, drsqrt, d2rsqrt, ACC1 );
        BBE_MULUUSNX16   (ACC0, drsqrt, d2rsqrt);
        ACC0 = BBE_SRAINX40(ACC0, 23);
        U0 = BBE_PACKLNX40(ACC0);
        U0_exp = BBE_SUBSR1SAVSN(20, norm_exp );
        U0 = BBE_SHFLNX16I(U0, BBE_SHFLI_DUPLICATE_1_EVEN);
        U0_exp = BBE_SHFLVSNI(U0_exp, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);

        ACC0 = BBE_MULUSRNX16(U0, V0, U0_exp);
        V0 = BBE_PACKVNX40(ACC0, U0_exp);
        BBE_SPNX16_IP(V0, pX0, 2*sz_i16);

        U0_exp = BBE_SUBSVSN(U0_exp, D0_exp);
        U1 = BBE_MOVVVS(U0_exp);
        T0 = BBE_SELNX16I(U1, U0, BBE_SELI_INTERLEAVE_1_EVEN);
        BBE_SPNX16_IP(T0, pXabs, 2*sz_i16);
    }

    // normalization
    pX0   = (xb_vecNx16 *)(x0);
    pXabs = (xb_vecNx16 *)(xabs);
    pX = (const xb_vecNx16 *)(x);
    al_pX = BBE_LANX16_PP(pX);
    pV = (xb_vecNx16 *)(v);
    al_pV = BBE_ZALIGN();
    __Pragma("loop_count min=1")
    for (l=0; l<L; l++)
    {
        BBE_LPNX16_IP(V0, pX0  , 2*sz_i16);
        BBE_LPNX16_IP(T0, pXabs, 2*sz_i16);

        D0 = BBE_REPNX16(T0, 0);
        T0 = BBE_REPNX16(T0, 1);
        D0_exp = BBE_MOVVSV(T0, 0);
        BBE_LAVNX16_XP(Xdiag, al_pX, pX, 2*sz_i16);
        BBE_LAVNX16_XP(X0, al_pX, pX, BBE_SIMD_WIDTH*sz_i16-2*sz_i16);
        BBE_LANX16_IP (X1, al_pX, pX);
        BBE_LAVNX16_XP(X2, al_pX, pX, NbM0);
        BBE_LAVNX16_XP(X3, al_pX, pX, NbM1);

        ACC0 = BBE_MULUSRNX16(D0, X0, D0_exp);
        ACC1 = BBE_MULUSRNX16(D0, X1, D0_exp);
        ACC2 = BBE_MULUSRNX16(D0, X2, D0_exp);
        ACC3 = BBE_MULUSRNX16(D0, X3, D0_exp);
        U0 = BBE_PACKVNX40(ACC0, D0_exp);
        U1 = BBE_PACKVNX40(ACC1, D0_exp);
        U2 = BBE_PACKVNX40(ACC2, D0_exp);
        U3 = BBE_PACKVNX40(ACC3, D0_exp);
        U0 = BBE_SELNX16I(U0, V0, BBE_SELI_PACK_2);
        BBE_SANX16_IP(U0, al_pV, pV);
        BBE_SANX16_IP(U1, al_pV, pV);
        BBE_SAVNX16_XP(U2, al_pV, pV, NbM0);
        BBE_SAVNX16_XP(U3, al_pV, pV, NbM1);
    }
    BBE_SAVNX16POS_FP(al_pV, pV);

}
#endif

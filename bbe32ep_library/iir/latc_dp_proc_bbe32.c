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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
    NatureDSP_Baseband library. IIR part
    Lattice complex block IIR w/ real coefficients
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Baseband_iir.h"
/* Common utility and macros declarations. */
#include "common.h"

#define sz_i16 sizeof(int16_t)

/* Lattice complex block IIR processing function, Low Noise fixed-point implementation. */

void latc_dp_proc2( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, cf0, cf1, r0, t0, t1;
    xb_vecNx16 dl0_l, dl0_h, dl1_l, dl1_h;
    int32_t cf01;
    xb_vecNx40 ACCR, ACCR_l, ACCR_h, ACCD;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    cf01 = (((int32_t)coef[0])<<16);
    cf01 |= (uint16_t)coef[1];
    cf0 = BBE_LSNX16_I(coef, 0*sz_i16);
    cf1 = BBE_LSNX16_I(coef, 1*sz_i16);
    cf0 = BBE_REPNX16(cf0, 0);
    cf1 = BBE_REPNX16(cf1, 0);
    // Load delay elements
    dl1_l = BBE_LPNX16_I(d, 0*sz_i16);
    dl1_h = BBE_LPNX16_I(d, 2*sz_i16);
    dl0_l = BBE_LPNX16_I(d, 4*sz_i16);
    dl0_h = BBE_LPNX16_I(d, 6*sz_i16);
    BBE_DSELNX16I(t0, t1, dl0_h, dl0_l, BBE_DSELI_INTERLEAVE_1);
    ACCD = BBE_MOVSWV(t0, t1);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q44 <- Q29*Q15 + Q29*Q15
        ACCR_l = BBE_MULUSRNX16(dl0_l, cf0, rnd15);
        BBE_MULUSANX16 (ACCR_l, dl1_l, cf1);
        ACCR_h = BBE_MULNX16PR(dl0_h, dl1_h, cf01);
        // Q29 <- Q44 - 15 w/ rounding
        ACCR_l = BBE_SRAINX40(ACCR_l, 15);
        ACCR_h = BBE_SLLINX40(ACCR_h, 1 );

        ACCR = BBE_SUBNX40(ACCR, ACCR_h);
        ACCR = BBE_SUBNX40(ACCR, ACCR_l);

        // Update the delay line
        // Q14 <- Q29 - 15 w/ rounding
        r0 = BBE_PACKQNX40(ACCR);
        // Q29 <- Q29 + Q14*Q15
        BBE_MULANX16(ACCD, r0, cf0);
        dl1_l = BBE_PACKLNX40(ACCD);
        dl1_h = BBE_PACKVNX40(ACCD, rnd16);
        ACCD = ACCR;
        dl0_l = BBE_PACKLNX40(ACCR);
        dl0_h = BBE_PACKVNX40(ACCR, rnd16);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SPNX16_I(dl1_l, d, 0*sz_i16);
    BBE_SPNX16_I(dl1_h, d, 2*sz_i16);
    BBE_SPNX16_I(dl0_l, d, 4*sz_i16);
    BBE_SPNX16_I(dl0_h, d, 6*sz_i16);

} // latc_dp_proc2()

void latc_dp_proc4( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r012, t0, t1;
    xb_vecNx16 dl0123_l, dl0_l, dl1_l, dl2_l, dl3_l;
    xb_vecNx16 dl0123_h, dl0_h, dl1_h, dl2_h, dl3_h;
    xb_vecNx16 cf0123, cf0, cf1, cf2, cf3;
    int32_t    cf23;
    xb_vecNx40 ACCR, ACCDL, ACC_l, ACC_h;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf23 = (((int32_t)coef[2])<<16);
    cf23 |= (uint16_t)coef[3];
    cf0123 = BBE_LVNX16_I(CF, 0);
    cf0 = BBE_REPNX16(cf0123, 0);
    cf1 = BBE_REPNX16(cf0123, 1);
    cf2 = BBE_REPNX16(cf0123, 2);
    cf3 = BBE_REPNX16(cf0123, 3);
    cf0123 = BBE_SELNX16I(cf0123, cf0123, BBE_SELI_INTERLEAVE_1_LO);
    cf0 = BBE_SELNX16I(zero, cf0, BBE_SELI_ROTATE_LEFT_4);
    cf0 = BBE_SELNX16I(cf0, cf0, BBE_SELI_ROTATE_RIGHT_2);
    cf1 = BBE_SELNX16I(zero, cf1, BBE_SELI_ROTATE_LEFT_6);
    cf1 = BBE_SELNX16I(cf1, cf1, BBE_SELI_ROTATE_RIGHT_2);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    dl0123_l = BBE_LVNX16_I(DL, 0*sz_i16*BBE_SIMD_WIDTH);
    dl0123_h = BBE_LVNX16_I(DL, 1*sz_i16*BBE_SIMD_WIDTH);
    BBE_DSELNX16I(t0, t1, dl0123_h, dl0123_l, BBE_DSELI_INTERLEAVE_1);
    ACCDL = BBE_MOVSWV(t0, t1);
    dl0_l = BBE_REPNX16C(dl0123_l, 0);
    dl1_l = BBE_REPNX16C(dl0123_l, 1);
    dl2_l = BBE_REPNX16C(dl0123_l, 2);
    dl3_l = BBE_REPNX16C(dl0123_l, 3);
    dl0_h = BBE_REPNX16C(dl0123_h, 0);
    dl1_h = BBE_REPNX16C(dl0123_h, 1);
    dl2_h = BBE_REPNX16C(dl0123_h, 2);
    dl3_h = BBE_REPNX16C(dl0123_h, 3);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        Xin = BBE_REPNX16C(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q44 <- Q44 + Q29*Q15
        ACC_l = BBE_MULUSRNX16(dl3_l, cf3, rnd15);
        BBE_MULUSANX16(ACC_l,  dl2_l, cf2);
        BBE_MULUSANX16(ACC_l,  dl1_l, cf1);
        BBE_MULUSANX16(ACC_l,  dl0_l, cf0);
        ACC_h = BBE_MULNX16PR(dl2_h, dl3_h, cf23);
        BBE_MULANX16  (ACC_h, dl1_h, cf1);
        BBE_MULANX16  (ACC_h, dl0_h, cf0);
        // Q29 <- Q44 - 15 w/ rounding
        ACC_l = BBE_SRAINX40(ACC_l, 15);
        ACC_h = BBE_SLLINX40(ACC_h, 1 );

        ACCR = BBE_SUBNX40(ACCR, ACC_l);
        ACCR = BBE_SUBNX40(ACCR, ACC_h);
        // Q14 <- Q29 - 15 w/ rounding
        r012 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        BBE_MULANX16(ACCDL, r012, cf0123);
        dl0123_l = BBE_PACKLNX40(ACCDL);
        dl0123_h = BBE_PACKVNX40(ACCDL, rnd16);
        dl3_l = BBE_REPNX16C(dl0123_l, 2);
        dl2_l = BBE_REPNX16C(dl0123_l, 1);
        dl1_l = BBE_REPNX16C(dl0123_l, 0);
        dl0_l = BBE_PACKLNX40(ACCR);
        dl3_h = BBE_REPNX16C(dl0123_h, 2);
        dl2_h = BBE_REPNX16C(dl0123_h, 1);
        dl1_h = BBE_REPNX16C(dl0123_h, 0);
        dl0_h = BBE_PACKVNX40(ACCR, rnd16);
        ACCDL = BBE_SELNX40I(ACCDL, ACCR, BBE_W_SELI_ROTATE_LEFT_2);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    dl0123_l = BBE_PACKLNX40(ACCDL);
    dl0123_h = BBE_PACKVNX40(ACCDL, rnd16);
    BBE_SVNX16_I(dl0123_l, DL, 0*sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(dl0123_h, DL, 1*sz_i16*BBE_SIMD_WIDTH);
} // latc_dp_proc4()

void latc_dp_proc6( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{

    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r01234;
    xb_vecNx16 dl012345_l, dl2345_l, dl45_l;
    xb_vecNx16 dl012345_h, dl2345_h, dl45_h;
    xb_vecNx16 cf012345, cf2345, cf45;
    xb_vecNx40 ACCR, ACCDL, ACC_l, ACC_h, ACC0, ACC1;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf012345 = BBE_LVNX16_I(CF, 0);
    cf012345 = BBE_SELNX16I(cf012345, cf012345, BBE_SELI_INTERLEAVE_1_LO);
    cf2345   = BBE_SELNX16I(zero, cf012345, BBE_SELI_ROTATE_RIGHT_4);
    cf45     = BBE_SELNX16I(zero, cf012345, BBE_SELI_ROTATE_RIGHT_8);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    dl012345_l = BBE_LVNX16_I(DL, 0*sz_i16*BBE_SIMD_WIDTH);
    dl012345_h = BBE_LVNX16_I(DL, 1*sz_i16*BBE_SIMD_WIDTH);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        Xin = BBE_REPNX16C(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);

        // Q44 <- Q44 + Q29*Q15
        dl2345_l = BBE_SELNX16I(dl012345_l, dl012345_l, BBE_SELI_ROTATE_RIGHT_4);
        dl2345_h = BBE_SELNX16I(dl012345_h, dl012345_h, BBE_SELI_ROTATE_RIGHT_4);
        dl45_l   = BBE_SELNX16I(dl012345_l, dl012345_l, BBE_SELI_ROTATE_RIGHT_8);
        dl45_h   = BBE_SELNX16I(dl012345_h, dl012345_h, BBE_SELI_ROTATE_RIGHT_8);

        ACC_l = BBE_MULUSRNX16(dl012345_l, cf012345 , rnd15);
        BBE_MULUSANX16 (ACC_l, dl2345_l  , cf2345  );
        BBE_MULUSANX16 (ACC_l, dl45_l    , cf45    );

        ACC_h = BBE_MULNX16(dl012345_h, cf012345);
        BBE_MULANX16(ACC_h, dl2345_h  , cf2345  );
        BBE_MULANX16(ACC_h, dl45_h    , cf45    );

        // Q29 <- Q44 - 15 w/ rounding
        ACC_l = BBE_SRAINX40(ACC_l , 15);
        ACC_h = BBE_SLLINX40(ACC_h , 1);
        ACC0 = BBE_ADDNX40(ACC_h , ACC_l );
        ACC1 = BBE_SELNX40I(ACC0, ACC0 , BBE_W_SELI_ROTATE_RIGHT_2);

        ACCR  = BBE_SUBNX40(ACCR , ACC0 );
        ACCR  = BBE_SUBNX40(ACCR , ACC1 );
        // Q14 <- Q29 - 15 w/ rounding
        r01234 = BBE_PACKQNX40(ACCR );

        // Update the delay line
        BBE_DSELNX16I(dl012345_h, dl012345_l, dl012345_h, dl012345_l, BBE_DSELI_INTERLEAVE_1);
        ACCDL = BBE_MOVSWV(dl012345_h, dl012345_l);
        BBE_MULANX16(ACCDL, r01234, cf012345);
        ACC0 = BBE_REPNX40C(ACCR, 0);
        ACCDL = BBE_SELNX40I(ACCDL, ACC0, BBE_W_SELI_ROTATE_LEFT_2);
        dl012345_l = BBE_PACKLNX40(ACCDL);
        dl012345_h = BBE_PACKVNX40(ACCDL, rnd16);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SVNX16_I(dl012345_l, DL, 0*sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(dl012345_h, DL, 1*sz_i16*BBE_SIMD_WIDTH);

} // latc_dp_proc6()

void latc_dp_proc8( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
#if 1
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r0123456;
    xb_vecNx16 dl01234567_l, dl234567_l, dl4567_l, dl67_l;
    xb_vecNx16 dl01234567_h, dl234567_h, dl4567_h, dl67_h;
    xb_vecNx16 cf01234567, cf234567, cf4567, cf67;
    xb_vecNx40 ACCR, ACCDL, ACC_l, ACC_h, ACC0, ACC1;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf01234567 = BBE_LVNX16_I(CF, 0);
    cf01234567 = BBE_SELNX16I(cf01234567, cf01234567, BBE_SELI_INTERLEAVE_1_LO);
    cf234567   = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_4);
    cf4567     = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_8);
    cf67       = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_12);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    dl01234567_l = BBE_LVNX16_I(DL, 0*sz_i16*BBE_SIMD_WIDTH);
    dl01234567_h = BBE_LVNX16_I(DL, 1*sz_i16*BBE_SIMD_WIDTH);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        Xin = BBE_REPNX16C(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q44 <- Q44 + Q29*Q15
        dl234567_l = BBE_SELNX16I(dl01234567_l, dl01234567_l, BBE_SELI_ROTATE_RIGHT_4);
        dl234567_h = BBE_SELNX16I(dl01234567_h, dl01234567_h, BBE_SELI_ROTATE_RIGHT_4);
        dl4567_l   = BBE_SELNX16I(dl01234567_l, dl01234567_l, BBE_SELI_ROTATE_RIGHT_8);
        dl4567_h   = BBE_SELNX16I(dl01234567_h, dl01234567_h, BBE_SELI_ROTATE_RIGHT_8);
        dl67_l     = BBE_SELNX16I(dl01234567_l, dl01234567_l, BBE_SELI_ROTATE_RIGHT_12);
        dl67_h     = BBE_SELNX16I(dl01234567_h, dl01234567_h, BBE_SELI_ROTATE_RIGHT_12);

        ACC_l = BBE_MULUSRNX16(dl01234567_l, cf01234567, rnd15);
        BBE_MULUSANX16 (ACC_l, dl234567_l  , cf234567  );
        BBE_MULUSANX16 (ACC_l, dl4567_l    , cf4567    );
        BBE_MULUSANX16 (ACC_l, dl67_l      , cf67      );

        ACC_h = BBE_MULNX16(dl01234567_h, cf01234567);
        BBE_MULANX16(ACC_h, dl234567_h  , cf234567  );
        BBE_MULANX16(ACC_h, dl4567_h    , cf4567    );
        BBE_MULANX16(ACC_h, dl67_h      , cf67      );
        // Q29 <- Q44 - 15 w/ rounding
        ACC_l = BBE_SRAINX40(ACC_l, 15);
        ACC_h = BBE_SLLINX40(ACC_h, 1);
        ACC0 = BBE_ADDNX40(ACC_h, ACC_l);
        ACC1 = BBE_SELNX40I(ACC0, ACC0, BBE_W_SELI_ROTATE_RIGHT_2);

        ACCR = BBE_SUBNX40(ACCR, ACC0);
        ACCR = BBE_SUBNX40(ACCR, ACC1);
        // Q14 <- Q29 - 15 w/ rounding
        r0123456 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        BBE_DSELNX16I(dl01234567_h, dl01234567_l, dl01234567_h, dl01234567_l, BBE_DSELI_INTERLEAVE_1);
        ACCDL = BBE_MOVSWV(dl01234567_h, dl01234567_l);
        BBE_MULANX16(ACCDL, r0123456, cf01234567);
        ACC0 = BBE_REPNX40C(ACCR, 0);
        ACCDL = BBE_SELNX40I(ACCDL, ACC0, BBE_W_SELI_ROTATE_LEFT_2);
        dl01234567_l = BBE_PACKLNX40(ACCDL);
        dl01234567_h = BBE_PACKVNX40(ACCDL, rnd16);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SVNX16_I(dl01234567_l, DL, 0*sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(dl01234567_h, DL, 1*sz_i16*BBE_SIMD_WIDTH);
#else
    xb_vecNx16 g, Ed1;

    xb_vecNx16 X, Y, Coef, Coef_, Zero_;

    xb_vecNx16 D_l0, D_h0, D_l1, D_h1;

    xb_vecNx40 T_, T_0, U_, D_, Zero; 

    const xb_vecNx16 *COEF = (const xb_vecNx16 *)coef;
          xb_vecNx16 *SECT = (      xb_vecNx16 *)d;

    vsaN rnd15 = BBE_MOVVSA32(15);
    vsaN rnd16 = BBE_MOVVSA32(16);

    int j;

    NASSERT( r && d && x && coef );
    NASSERT_ALIGN4(r);
    NASSERT_ALIGN4(x);

    Zero  = BBE_MOVQINT40(0);
    Zero_ = BBE_MOVQINT16(0);
    
    g = BBE_MOVVA16(gain);
    Ed1 = BBE_MOVVINT16(1);

    Coef = BBE_LVNX16_I(COEF, 0);
    Coef = BBE_SELNX16I(Coef, Coef, BBE_SELI_INTERLEAVE_1_LO);

    Coef_ = BBE_SELNX16I(Zero_, Coef, BBE_SELI_ROTATE_RIGHT_8);

    D_l1 = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
    D_h1 = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);

    for ( j=0; j<N; j++ )
    {
        BBE_LPNX16_IP(X, x, 4);
        X = BBE_REPNX16C(X, 0);

        T_ = BBE_MULUSNX16(g, X);

        // Q30 <- Q30 - ( Q30*Q15 - 15 w/ rounding )
        U_ = BBE_MULUSRNX16(D_l1, Coef, rnd15);
        U_ = BBE_SRAINX40(U_, 15);

        BBE_MULANX16(U_, D_h1, Coef);
        BBE_MULANX16(U_, D_h1, Coef);

        D_l0 = BBE_SELNX16I(Zero_, D_l1, BBE_SELI_ROTATE_RIGHT_8);
        D_h0 = BBE_SELNX16I(Zero_, D_h1, BBE_SELI_ROTATE_RIGHT_8);

        // Q30 <- Q30 - ( Q30*Q15 - 15 w/ rounding )
        T_0 = BBE_MULUSRNX16(D_l0, Coef_, rnd15);
        T_0 = BBE_SRAINX40(T_0, 15);

        BBE_MULANX16(T_0, D_h0, Coef_);
        BBE_MULANX16(T_0, D_h0, Coef_);

        U_ = BBE_ADDNX40(T_0, U_);
        T_ = BBE_SUBNX40(T_, U_);

        U_ = BBE_SELNX40I(Zero, U_, BBE_W_SELI_ROTATE_RIGHT_2);
        T_ = BBE_SUBNX40(T_, U_);

        U_ = BBE_SELNX40I(Zero, U_, BBE_W_SELI_ROTATE_RIGHT_2);
        T_ = BBE_SUBNX40(T_, U_);

        U_ = BBE_SELNX40I(Zero, U_, BBE_W_SELI_ROTATE_RIGHT_2);
        T_ = BBE_SUBNX40(T_, U_);

        Y = BBE_PACKQNX40(T_);

        BBE_SPNX16_IP(Y, r, 4);

        D_ = BBE_UNPKSNX16(D_h1);
        D_ = BBE_SLLINX40(D_, 16);
        BBE_MULUSANX16(D_, D_l1, Ed1);

        BBE_MULANX16(D_, Y, Coef);

        U_ = BBE_REPNX40C(T_, 0);
        D_ = BBE_SELNX40I(D_, U_, BBE_W_SELI_ROTATE_LEFT_2);

        D_l1 = BBE_PACKLNX40(D_); D_h1 = BBE_PACKVNX40(D_, rnd16);
    }

    BBE_SVNX16_I(D_l1, SECT, 0*2*BBE_SIMD_WIDTH);
    BBE_SVNX16_I(D_h1, SECT, 1*2*BBE_SIMD_WIDTH);
#endif
} // latc_dp_proc8()

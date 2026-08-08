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
    Lattice real block IIR
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

/* Lattice real block IIR processing function, double precision. */
void latr_dp_proc2( int16_t * restrict r,
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
    // Load delay elements
    dl1_l = BBE_LSNX16_I(d, 0*sz_i16);
    dl1_h = BBE_LSNX16_I(d, 1*sz_i16);
    dl0_l = BBE_LSNX16_I(d, 2*sz_i16);
    dl0_h = BBE_LSNX16_I(d, 3*sz_i16);
    BBE_DSELNX16I(t0, t1, dl0_h, dl0_l, BBE_DSELI_INTERLEAVE_1);
    ACCD = BBE_MOVSWV(t0, t1);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LSNX16_IP(Xin, x, sz_i16);
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
        BBE_SSNX16_IP(Rout, r, sz_i16);
    }
    // Save delay elements
    BBE_SSNX16_I(dl1_l, d, 0*sz_i16);
    BBE_SSNX16_I(dl1_h, d, 1*sz_i16);
    BBE_SSNX16_I(dl0_l, d, 2*sz_i16);
    BBE_SSNX16_I(dl0_h, d, 3*sz_i16);

} // latr_dp_proc2()

void latr_dp_proc4( int16_t * restrict r,
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
    cf0 = BBE_SELNX16I(zero, cf0, BBE_SELI_ROTATE_LEFT_2);
    cf0 = BBE_SELNX16I(cf0, cf0, BBE_SELI_ROTATE_RIGHT_1);
    cf1 = BBE_SELNX16I(zero, cf1, BBE_SELI_ROTATE_LEFT_3);
    cf1 = BBE_SELNX16I(cf1, cf1, BBE_SELI_ROTATE_RIGHT_1);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    t0 = BBE_LVNX16_I(DL, 0);
    dl0123_l = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_LO_HALVES);
    dl0123_h = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_HI_HALVES);
    BBE_DSELNX16I(t0, t1, dl0123_h, dl0123_l, BBE_DSELI_INTERLEAVE_1);
    ACCDL = BBE_MOVSWV(t0, t1);
    dl0_l = BBE_REPNX16(dl0123_l, 0);
    dl1_l = BBE_REPNX16(dl0123_l, 1);
    dl2_l = BBE_REPNX16(dl0123_l, 2);
    dl3_l = BBE_REPNX16(dl0123_l, 3);
    dl0_h = BBE_REPNX16(dl0123_h, 0);
    dl1_h = BBE_REPNX16(dl0123_h, 1);
    dl2_h = BBE_REPNX16(dl0123_h, 2);
    dl3_h = BBE_REPNX16(dl0123_h, 3);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LSNX16_IP(Xin, x, sz_i16);
        Xin = BBE_REPNX16(Xin, 0);
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
        dl3_l = BBE_REPNX16(dl0123_l, 2);
        dl2_l = BBE_REPNX16(dl0123_l, 1);
        dl1_l = BBE_REPNX16(dl0123_l, 0);
        dl0_l = BBE_PACKLNX40(ACCR);
        dl3_h = BBE_REPNX16(dl0123_h, 2);
        dl2_h = BBE_REPNX16(dl0123_h, 1);
        dl1_h = BBE_REPNX16(dl0123_h, 0);
        dl0_h = BBE_PACKVNX40(ACCR, rnd16);
        ACCDL = BBE_SELNX40I(ACCDL, ACCR, BBE_W_SELI_ROTATE_LEFT_1);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SSNX16_IP(Rout, r, sz_i16);
    }
    // Save delay elements
    dl0123_l = BBE_PACKLNX40(ACCDL);
    dl0123_h = BBE_PACKVNX40(ACCDL, rnd16);
    t0 = BBE_SELNX16I(dl0123_h, dl0123_l, BBE_SELI_EXTRACT_LO_HALVES);
    BBE_SVNX16_I(t0, DL, 0);
} // latr_dp_proc4()

void latr_dp_proc6( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r01234, t0, t1;
    xb_vecNx16 dl012345_l, dl2345_l, dl45_l;
    xb_vecNx16 dl012345_h, dl2345_h, dl45_h;
    xb_vecNx16 cf012345, cf2345, cf45;
    xb_vecNx40 ACCR, ACCDL, ACC_l, ACC_h, ACC0, ACC1, ZERO;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    ZERO = BBE_ZERONX40();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf012345 = BBE_LVNX16_I(CF, 0);
    cf2345 = BBE_SELNX16I(zero, cf012345, BBE_SELI_ROTATE_RIGHT_2);
    cf45   = BBE_SELNX16I(zero, cf012345, BBE_SELI_ROTATE_RIGHT_4);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    t0 = BBE_LVNX16_I(DL, 0);
    dl012345_l = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_LO_HALVES);
    dl012345_h = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_HI_HALVES);
    BBE_DSELNX16I(t0, t1, dl012345_h, dl012345_l, BBE_DSELI_INTERLEAVE_1);
    ACCDL = BBE_MOVSWV(t0, t1);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LSNX16_IP(Xin, x, sz_i16);
        Xin = BBE_REPNX16(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q44 <- Q44 + Q29*Q15
        dl2345_l = BBE_SELNX16I(zero, dl012345_l, BBE_SELI_ROTATE_RIGHT_2);
        dl2345_h = BBE_SELNX16I(zero, dl012345_h, BBE_SELI_ROTATE_RIGHT_2);
        dl45_l   = BBE_SELNX16I(zero, dl012345_l, BBE_SELI_ROTATE_RIGHT_4);
        dl45_h   = BBE_SELNX16I(zero, dl012345_h, BBE_SELI_ROTATE_RIGHT_4);
        ACC_l = BBE_MULUSRNX16(dl012345_l, cf012345, rnd15);
        BBE_MULUSANX16 (ACC_l, dl2345_l  , cf2345  );
        BBE_MULUSANX16 (ACC_l, dl45_l    , cf45    );
        ACC_h = BBE_MULNX16(dl012345_h, cf012345);
        BBE_MULANX16(ACC_h, dl2345_h  , cf2345  );
        BBE_MULANX16(ACC_h, dl45_h    , cf45    );
        // Q29 <- Q44 - 15 w/ rounding
        ACC_l = BBE_SRAINX40(ACC_l, 15);
        ACC_h = BBE_SLLINX40(ACC_h, 1 );
        ACC0 = BBE_ADDNX40(ACC_h, ACC_l);
        ACC1 = BBE_SELNX40I(ZERO, ACC0, BBE_W_SELI_ROTATE_RIGHT_1);

        ACCR = BBE_SUBNX40(ACCR, ACC0);
        ACCR = BBE_SUBNX40(ACCR, ACC1);
        // Q14 <- Q29 - 15 w/ rounding
        r01234 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        BBE_MULANX16(ACCDL, r01234, cf012345);
        ACC0 = BBE_REPNX40(ACCR, 0);
        ACCDL = BBE_SELNX40I(ACCDL, ACC0, BBE_W_SELI_ROTATE_LEFT_1);
        dl012345_l = BBE_PACKLNX40(ACCDL);
        dl012345_h = BBE_PACKVNX40(ACCDL, rnd16);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SSNX16_IP(Rout, r, sz_i16);
    }
    // Save delay elements
    dl012345_l = BBE_PACKLNX40(ACCDL);
    dl012345_h = BBE_PACKVNX40(ACCDL, rnd16);
    t0 = BBE_SELNX16I(dl012345_h, dl012345_l, BBE_SELI_EXTRACT_LO_HALVES);
    BBE_SVNX16_I(t0, DL, 0);
} // latr_dp_proc6()

void latr_dp_proc8( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r0123456, t0, t1;
    xb_vecNx16 dl01234567_l, dl234567_l, dl4567_l, dl67_l;
    xb_vecNx16 dl01234567_h, dl234567_h, dl4567_h, dl67_h;
    xb_vecNx16 cf01234567, cf234567, cf4567, cf67;
    xb_vecNx40 ACCR, ACCDL, ACC_l, ACC_h, ACC0, ACC1, ZERO;
    int n;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd15 = BBE_MOVVSA32(15);

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    ZERO = BBE_ZERONX40();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf01234567 = BBE_LVNX16_I(CF, 0);
    cf234567 = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_2);
    cf4567   = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_4);
    cf67     = BBE_SELNX16I(zero, cf01234567, BBE_SELI_ROTATE_RIGHT_6);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    t0 = BBE_LVNX16_I(DL, 0);
    dl01234567_l = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_LO_HALVES);
    dl01234567_h = BBE_SELNX16I(t0, t0, BBE_SELI_EXTRACT_HI_HALVES);
    BBE_DSELNX16I(t0, t1, dl01234567_h, dl01234567_l, BBE_DSELI_INTERLEAVE_1);
    ACCDL = BBE_MOVSWV(t0, t1);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LSNX16_IP(Xin, x, sz_i16);
        Xin = BBE_REPNX16(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q44 <- Q44 + Q29*Q15
        dl234567_l = BBE_SELNX16I(zero, dl01234567_l, BBE_SELI_ROTATE_RIGHT_2);
        dl234567_h = BBE_SELNX16I(zero, dl01234567_h, BBE_SELI_ROTATE_RIGHT_2);
        dl4567_l   = BBE_SELNX16I(zero, dl01234567_l, BBE_SELI_ROTATE_RIGHT_4);
        dl4567_h   = BBE_SELNX16I(zero, dl01234567_h, BBE_SELI_ROTATE_RIGHT_4);
        dl67_l     = BBE_SELNX16I(zero, dl01234567_l, BBE_SELI_ROTATE_RIGHT_6);
        dl67_h     = BBE_SELNX16I(zero, dl01234567_h, BBE_SELI_ROTATE_RIGHT_6);
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
        ACC_h = BBE_SLLINX40(ACC_h, 1 );
        ACC0 = BBE_ADDNX40(ACC_h, ACC_l);
        ACC1 = BBE_SELNX40I(ZERO, ACC0, BBE_W_SELI_ROTATE_RIGHT_1);

        ACCR = BBE_SUBNX40(ACCR, ACC0);
        ACCR = BBE_SUBNX40(ACCR, ACC1);
        // Q14 <- Q29 - 15 w/ rounding
        r0123456 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        BBE_MULANX16(ACCDL, r0123456, cf01234567);
        ACC0 = BBE_REPNX40(ACCR, 0);
        ACCDL = BBE_SELNX40I(ACCDL, ACC0, BBE_W_SELI_ROTATE_LEFT_1);
        dl01234567_l = BBE_PACKLNX40(ACCDL);
        dl01234567_h = BBE_PACKVNX40(ACCDL, rnd16);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SSNX16_IP(Rout, r, sz_i16);
    }
    // Save delay elements
    dl01234567_l = BBE_PACKLNX40(ACCDL);
    dl01234567_h = BBE_PACKVNX40(ACCDL, rnd16);
    t0 = BBE_SELNX16I(dl01234567_h, dl01234567_l, BBE_SELI_EXTRACT_LO_HALVES);
    BBE_SVNX16_I(t0, DL, 0);
} // latr_dp_proc8()

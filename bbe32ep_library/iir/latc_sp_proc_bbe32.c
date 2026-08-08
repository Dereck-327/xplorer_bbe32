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

/* Lattice complex block IIR processing function, Fast fixed-point implementation. */

void latc_sp_proc2( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, cf0, r0, dl0, dl1;
    int32_t cf01;
    xb_vecNx40 ACCR;
    int n;

    NASSERT( r && d && x && coef );
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    cf01 = (((int32_t)(-coef[0]))<<16);
    cf01 |= (uint16_t)(-coef[1]);
    cf0 = BBE_LSNX16_I(coef, 0);
    cf0 = BBE_REPNX16(cf0, 0);
    // Load delay elements
    dl1 = BBE_LPNX16_I(d, 0);
    dl0 = BBE_LPNX16_I(d, 2*sz_i16);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q29 <- Q29 + Q14*Q15 + Q14*Q15
        BBE_MULANX16PR(ACCR, dl0, dl1, cf01);
        // Q14 <- Q29 - 15 w/ rounding
        r0 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        dl1 = dl0 + BBE_MULNX16PACKQ(r0, cf0);
        dl0 = r0;

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SPNX16_I(dl1, d, 0);
    BBE_SPNX16_I(dl0, d, 2*sz_i16);

} // latc_sp_proc2()

void latc_sp_proc4( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, cf0, cf1, cf2, r0, r1, r2, dl0, dl1, dl2, dl3;
    int32_t cf23;
    xb_vecNx40 ACCR;
    int n;

    NASSERT( r && d && x && coef );
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    cf23 = (((int32_t)(-coef[2]))<<16);
    cf23 |= (uint16_t)(-coef[3]);
    cf0 = BBE_LSNX16_I(coef, 0*sz_i16);
    cf1 = BBE_LSNX16_I(coef, 1*sz_i16);
    cf2 = BBE_LSNX16_I(coef, 2*sz_i16);
    cf0 = BBE_REPNX16(cf0, 0);
    cf1 = BBE_REPNX16(cf1, 0);
    cf2 = BBE_REPNX16(cf2, 0);
    // Load delay elements
    dl3 = BBE_LPNX16_I(d, 0*sz_i16);
    dl2 = BBE_LPNX16_I(d, 2*sz_i16);
    dl1 = BBE_LPNX16_I(d, 4*sz_i16);
    dl0 = BBE_LPNX16_I(d, 6*sz_i16);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q29 <- Q29 + Q14*Q15 + Q14*Q15
        BBE_MULANX16PR(ACCR, dl2, dl3, cf23);
        // Q14 <- Q29 - 15 w/ rounding
        r2 = BBE_PACKQNX40(ACCR);
        BBE_MULSNX16(ACCR, dl1, cf1);
        r1 = BBE_PACKQNX40(ACCR);
        BBE_MULSNX16(ACCR, dl0, cf0);
        r0 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        dl3 = dl2 + BBE_MULNX16PACKQ(r2, cf2);
        dl2 = dl1 + BBE_MULNX16PACKQ(r1, cf1);
        dl1 = dl0 + BBE_MULNX16PACKQ(r0, cf0);
        dl0 = r0;

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SPNX16_I(dl3, d, 0*sz_i16);
    BBE_SPNX16_I(dl2, d, 2*sz_i16);
    BBE_SPNX16_I(dl1, d, 4*sz_i16);
    BBE_SPNX16_I(dl0, d, 6*sz_i16);
} // latc_sp_proc4()

void latc_sp_proc6( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r01234, dl012345;
    xb_vecNx16 dl0, dl1, dl2, dl3, dl4, dl5;
    xb_vecNx16 cf012345, cf0, cf1, cf2, cf3;
    int32_t    cf45;
    xb_vecNx40 ACCR;
    int n;

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf012345 = BBE_LVNX16_I(CF, 0);
    cf0 = BBE_REPNX16(cf012345, 0);
    cf1 = BBE_REPNX16(cf012345, 1);
    cf2 = BBE_REPNX16(cf012345, 2);
    cf3 = BBE_REPNX16(cf012345, 3);
    cf012345 = BBE_SELNX16I(cf012345, cf012345, BBE_SELI_INTERLEAVE_1_LO);
    cf0 = BBE_SELNX16I(zero, cf0, BBE_SELI_ROTATE_LEFT_2);
    cf1 = BBE_SELNX16I(zero, cf1, BBE_SELI_ROTATE_LEFT_4);
    cf2 = BBE_SELNX16I(zero, cf2, BBE_SELI_ROTATE_LEFT_6);
    cf3 = BBE_SELNX16I(zero, cf3, BBE_SELI_ROTATE_LEFT_8);
    cf45 = (((int32_t)(-coef[4]))<<16);
    cf45 |= (uint16_t)(-coef[5]);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    dl012345 = BBE_LVNX16_I(DL, 0);
    dl0 = BBE_REPNX16C(dl012345, 0);
    dl1 = BBE_REPNX16C(dl012345, 1);
    dl2 = BBE_REPNX16C(dl012345, 2);
    dl3 = BBE_REPNX16C(dl012345, 3);
    dl4 = BBE_REPNX16C(dl012345, 4);
    dl5 = BBE_REPNX16C(dl012345, 5);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        Xin = BBE_REPNX16C(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q29 <- Q29 + Q14*Q15 + Q14*Q15
        BBE_MULANX16PR(ACCR, dl4, dl5, cf45);
        BBE_MULSNX16  (ACCR, dl3,      cf3);
        BBE_MULSNX16  (ACCR, dl2,      cf2);
        BBE_MULSNX16  (ACCR, dl1,      cf1);
        BBE_MULSNX16  (ACCR, dl0,      cf0);
        // Q14 <- Q29 - 15 w/ rounding
        r01234 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        dl012345 = dl012345 + BBE_MULNX16PACKQ(r01234, cf012345);
        dl5 = BBE_REPNX16C(dl012345, 4);
        dl4 = BBE_REPNX16C(dl012345, 3);
        dl3 = BBE_REPNX16C(dl012345, 2);
        dl2 = BBE_REPNX16C(dl012345, 1);
        dl1 = BBE_REPNX16C(dl012345, 0);
        dl0 = r01234;
        dl012345 = BBE_SELNX16I(dl012345, r01234, BBE_SELI_PACK_2);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SVNX16_I(dl012345, DL, 0);
} // latc_sp_proc6()

void latc_sp_proc8( int16_t * restrict r,
                    int16_t * restrict d,
              const int16_t *          x,
              const int16_t *          coef,
                    int16_t            gain,
                    int N )
{
    const xb_vecNx16 * restrict CF;
          xb_vecNx16 * restrict DL;
    xb_vecNx16 Xin, Rout;
    xb_vecNx16 g, zero, r0123456, dl01234567;
    xb_vecNx16 dl0, dl1, dl2, dl3, dl4, dl5, dl6, dl7;
    xb_vecNx16 cf01234567, cf0, cf1, cf2, cf3, cf4, cf5;
    int32_t    cf67;
    xb_vecNx40 ACCR;
    int n;

    NASSERT( r && d && x && coef );

    zero = BBE_ZERONX16();
    // Load the input gain and reflection coefficients
    g = BBE_MOVVA16(gain);
    CF = (xb_vecNx16 *)coef;
    cf01234567 = BBE_LVNX16_I(CF, 0);
    cf0 = BBE_REPNX16(cf01234567, 0);
    cf1 = BBE_REPNX16(cf01234567, 1);
    cf2 = BBE_REPNX16(cf01234567, 2);
    cf3 = BBE_REPNX16(cf01234567, 3);
    cf4 = BBE_REPNX16(cf01234567, 4);
    cf5 = BBE_REPNX16(cf01234567, 5);
    cf01234567 = BBE_SELNX16I(cf01234567, cf01234567, BBE_SELI_INTERLEAVE_1_LO);
    cf0 = BBE_SELNX16I(zero, cf0, BBE_SELI_ROTATE_LEFT_2);
    cf1 = BBE_SELNX16I(zero, cf1, BBE_SELI_ROTATE_LEFT_4);
    cf2 = BBE_SELNX16I(zero, cf2, BBE_SELI_ROTATE_LEFT_6);
    cf3 = BBE_SELNX16I(zero, cf3, BBE_SELI_ROTATE_LEFT_8);
    cf4 = BBE_SELNX16I(zero, cf4, BBE_SELI_ROTATE_LEFT_10);
    cf5 = BBE_SELNX16I(zero, cf5, BBE_SELI_ROTATE_LEFT_12);
    cf67 = (((int32_t)(-coef[6]))<<16);
    cf67 |= (uint16_t)(-coef[7]);
    // Load delay elements
    DL = (xb_vecNx16 *)d;
    dl01234567 = BBE_LVNX16_I(DL, 0);
    dl0 = BBE_REPNX16C(dl01234567, 0);
    dl1 = BBE_REPNX16C(dl01234567, 1);
    dl2 = BBE_REPNX16C(dl01234567, 2);
    dl3 = BBE_REPNX16C(dl01234567, 3);
    dl4 = BBE_REPNX16C(dl01234567, 4);
    dl5 = BBE_REPNX16C(dl01234567, 5);
    dl6 = BBE_REPNX16C(dl01234567, 6);
    dl7 = BBE_REPNX16C(dl01234567, 7);

    __Pragma("loop_count min=1")
    for (n = 0; n < N; n++)
    {
        BBE_LPNX16_IP(Xin, x, 2*sz_i16);
        Xin = BBE_REPNX16C(Xin, 0);
        // Compute the sample
        // Q29 <- Q15*Q15 - 1
        ACCR = BBE_MULNX16(Xin, g);
        ACCR = BBE_SRAINX40(ACCR, 1);
        // Q29 <- Q29 + Q14*Q15 + Q14*Q15
        BBE_MULANX16PR(ACCR, dl6, dl7, cf67);
        BBE_MULSNX16  (ACCR, dl5,      cf5);
        BBE_MULSNX16  (ACCR, dl4,      cf4);
        BBE_MULSNX16  (ACCR, dl3,      cf3);
        BBE_MULSNX16  (ACCR, dl2,      cf2);
        BBE_MULSNX16  (ACCR, dl1,      cf1);
        BBE_MULSNX16  (ACCR, dl0,      cf0);
        // Q14 <- Q29 - 15 w/ rounding
        r0123456 = BBE_PACKQNX40(ACCR);

        // Update the delay line
        dl01234567 = dl01234567 + BBE_MULNX16PACKQ(r0123456, cf01234567);
        dl7 = BBE_REPNX16C(dl01234567, 6);
        dl6 = BBE_REPNX16C(dl01234567, 5);
        dl5 = BBE_REPNX16C(dl01234567, 4);
        dl4 = BBE_REPNX16C(dl01234567, 3);
        dl3 = BBE_REPNX16C(dl01234567, 2);
        dl2 = BBE_REPNX16C(dl01234567, 1);
        dl1 = BBE_REPNX16C(dl01234567, 0);
        dl0 = r0123456;
        dl01234567 = BBE_SELNX16I(dl01234567, r0123456, BBE_SELI_PACK_2);

        // Format and store the output sample
        // Q15 <- Q29 - 14 w/ rounding
        ACCR = BBE_SLLINX40(ACCR, 1);
        Rout = BBE_PACKQNX40(ACCR);
        BBE_SPNX16_IP(Rout, r, 2*sz_i16);
    }
    // Save delay elements
    BBE_SVNX16_I(dl01234567, DL, 0);
} // latc_sp_proc8()


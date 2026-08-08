/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Communications
    zcseq_gen2()
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

/*-------------------------------------------------------------------------
  zcseq_gen2()

  Operation : This function generates u'th root Zadoff-Chu(ZC) sequence
              as defined in section 5.7.2 of 3GPP 36.211 V8.4.0(2008-09)
              xu(n) = e^(-j*pi*u*n*(n+1)/Nzc)
              where 0 <= n <= Nzc-1

  Accuracy: 44 (1.3e-3)

  Output:
    r    Pointer to output complex ZC sequence.
         Real and imaginary outputs are interleaved.
         Output is in Q15 format.
  Input:
    u    Root no. to be used for ZC sequence generation
         Range should be (1 <= u <= Nzc-1)
    Nzc  Length of ZC sequence (839 or 139)

  Restrictions:
    r[]  Must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

#define inzc_139  15449523  // 2/139 in Q30
#define inzc_839  655251268  // 2/839 in Q38
static const int16_t ALIGN(32) sincos_table[8][BBE_SIMD_WIDTH] = {
    /*Cosine*/
    /* cos(pi*(2n+1)/16)*32768, n=8..15,0..7 */
    { 32138, 27246, 18205, 6393, -6393, -18205, -27246, -32138,
    -32138, -27246, -18205, -6393, 6393, 18205, 27246, 32138 },
    /* -(1/4)*pi^1*sin(pi*(2n+1)/16)*32768,   n=8..15,0..7 */
    { -5021, -14298, -21399, -25241, -25241, -21399, -14298, -5021,
    5021, 14298, 21399, 25241, 25241, 21399, 14298, 5021 },
    /* -(1/32)*pi^2*cos(pi*(2n+1)/16)*32768, n=8..15,0..7 */
    { -9912, -8403, -5615, -1972, 1972, 5615, 8403, 9912,
    9912, 8403, 5615, 1972, -1972, -5615, -8403, -9912 },
    /* (1/384)*pi^3*sin(pi*(2n+1)/16)*32768, n=8..15,0..7 */
    { 516, 1470, 2200, 2595, 2595, 2200, 1470, 516,
    -516, -1470, -2200, -2595, -2595, -2200, -1470, -516 },
    /* (1/4)*pi^1*cos(pi*(2n+1)/16)*32768,   n=8..15,0..7 */
    /*Sine*/
    { 6393, 18204, 27246, 32138, 32138, 27246, 18204, 6393,
    -6393, -18204, -27246, -32138, -32138, -27246, -18204, -6393 },
    /* (1/8)*pi^1*cos(pi*(2n+1)/32)*32768 */
    { 25241, 21399, 14298, 5021, -5021, -14298, -21399, -25241,
    -25241, -21399, -14298, -5021, 5021, 14298, 21399, 25241 },
    /* -(1/32)*pi^2*sin(pi*(2n+1)/16)*32768, n=8..15,0..7 */
    { -1972, -5615, -8403, -9912, -9912, -8403, -5615, -1972,
    1972, 5615, 8403, 9912, 9912, 8403, 5615, 1972 },
    /* -(1/384)*pi^3*cos(pi*(2n+1)/16)*32768, n=8..15,0..7 */
    { -2595, -2200, -1470, -516, 516, 1470, 2200, 2595,
    2595, 2200, 1470, 516, -516, -1470, -2200, -2595 }
};

void zcseq_gen2 ( int16_t * restrict r, const uint16_t u, uint16_t Nzc )
{
    static const int16_t ALIGN(32) arith_prog[16] =
    {
        0,
        0 + 1,
        0 + 1 + 2,
        0 + 1 + 2 + 3,
        0 + 1 + 2 + 3 + 4,
        0 + 1 + 2 + 3 + 4 + 5,
        0 + 1 + 2 + 3 + 4 + 5 + 6,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12 + 13,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12 + 13 + 14,
        0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12 + 13 + 14 + 15
    };
    xb_vecNx16    * restrict T = (xb_vecNx16 *)sincos_table;
    /* Q0.30 x Q10.0
    * max value = 2*(nzc-1)/nzc which is less than 2, hence fits in Q1.29
    */
    uint32_t inc_tmp;
    xb_vecNx16 * pr = (xb_vecNx16 *)(r);
    xb_vecNx40 P, PINC, PTMP;
    valign  r_align;
    int i;
    int inc, inc_17;
    vsaN shr15;

    xb_vecNx16 PH, CS, SN, x_srl, tmp;
    vselN sel;
    xb_vecNx16 TABLE0, TABLE1, TABLE2, TABLE3;
    xb_vecNx16 LO, HI;

    const xb_vecNx16 * restrict t0_ptr;
    const xb_vecNx16 * restrict t1_ptr;
    const xb_vecNx16 * restrict t2_ptr;
    const xb_vecNx16 * restrict t3_ptr;

    shr15 = BBE_MOVVSA32(15);
    inc_tmp = (Nzc == 839) ? 655251268UL : 3955077798UL;
    inc = ((int64_t)inc_tmp * u) >> 8; // multiply always fit in 40 bit range
    inc_17 = ((inc << 4) + inc + 1) >> 1; //((inc << 3) + inc + 1) >> 1 ;


    /* for NZC = 839, we have Q30 format
    while for NZC = 139 we have Q22 format
    */
    {
        xb_vecNx16 t, lo, hi;
        t = BBE_SEQNX16();
        PINC = BBE_MOVWA32(inc);
        lo = BBE_PACKLNX40(PINC);
        PINC = BBE_SRAINX40(PINC, 16);
        hi = BBE_PACKLNX40(PINC);
        PINC = BBE_MULNX16(hi, t);
        PINC = BBE_SLLINX40(PINC, 16);
        BBE_MULUSANX16(PINC, lo, t);

        t = BBE_LVNX16_I((const xb_vecNx16 *)arith_prog, 0);
        P = BBE_MOVWA32(inc);
        lo = BBE_PACKLNX40(P);
        P = BBE_SRAINX40(P, 16);
        hi = BBE_PACKLNX40(P);
        P = BBE_MULNX16(hi, t);
        P = BBE_SLLINX40(P, 16);
        BBE_MULUSANX16(P, lo, t);
    }

    /** for n > 15
    theta[n] = theta[n-16] + 16{(delta[n] + 17p/2}

    delta[n] = delta[n-16] + 16p
    p  = 2u/Nzc

    theta[n] = theta mod 2.0
    **/
    for (i = 0; i<(Nzc&(~15)); i += 16)
    {
        PTMP = BBE_RNDADJNX40(P, shr15);
        PTMP = BBE_SRLINX40(PTMP, 15);
        PH = BBE_PACKLNX40(PTMP);

        PTMP = BBE_MOVWA32(inc_17);
        XT_MOVEQZ(inc_17, inc_17, inc_17);
        PTMP = BBE_ADDNX40(PINC, PTMP);
        PTMP = BBE_SLLINX40(PTMP, 4);
        P = BBE_ADDNX40(P, PTMP);

        //COMPUTE new inc
        PTMP = BBE_MOVWA32(inc);
        //inc+=zero;
        XT_MOVEQZ(inc, inc, inc);
        PTMP = BBE_SLLINX40(PTMP, 4);
        PINC = BBE_ADDNX40(PINC, PTMP);

        /*
        exp(-j*pi*q*m.*(m+1)/Nrszc) =
        cos(phases[m]) - j*sin(phases[m]);
        */
        CS = BBE_ZERONX16();
        SN = BBE_ZERONX16();
        t0_ptr = (const xb_vecNx16 *)T + 0;
        t1_ptr = (const xb_vecNx16 *)T + 1;
        t2_ptr = (const xb_vecNx16 *)T + 2;
        t3_ptr = (const xb_vecNx16 *)T + 3;

        TABLE3 = BBE_LVNX16_I(t3_ptr, 0);
        TABLE2 = BBE_LVNX16_I(t2_ptr, 0);
        TABLE1 = BBE_LVNX16_I(t1_ptr, 0);
        TABLE0 = BBE_LVNX16_I(t0_ptr, 0);
        //PH = BBE_POLY8X20_SU(PH, 12, 2);
        // ---- POLY.SU --- //
        x_srl = BBE_SRLINX16(PH, 12);
        sel = xb_vecNx16_rtor_vselN(x_srl);
        //x_srl = x_srl - tmp1 ;
        PH = BBE_POLYNX16_OFF(PH, 12, 2);
        // ---- POLY.SU --- //

        // ---- POLY.STEP --- //
        CS = BBE_SHFLNX16(TABLE3, sel);
        tmp = BBE_SHFLNX16(TABLE2, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        tmp = BBE_SHFLNX16(TABLE1, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        tmp = BBE_SHFLNX16(TABLE0, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        // ---- POLY.STEP --- //
        t0_ptr = (const xb_vecNx16 *)T + 4;
        t1_ptr = (const xb_vecNx16 *)T + 5;
        t2_ptr = (const xb_vecNx16 *)T + 6;
        t3_ptr = (const xb_vecNx16 *)T + 7;

        TABLE3 = BBE_LVNX16_I(t3_ptr, 0);
        TABLE2 = BBE_LVNX16_I(t2_ptr, 0);
        TABLE1 = BBE_LVNX16_I(t1_ptr, 0);
        TABLE0 = BBE_LVNX16_I(t0_ptr, 0);

        // ---- POLY.STEP --- //
        SN = BBE_SHFLNX16(TABLE3, sel);
        tmp = BBE_SHFLNX16(TABLE2, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        tmp = BBE_SHFLNX16(TABLE1, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        tmp = BBE_SHFLNX16(TABLE0, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        // ---- POLY.STEP --- //
        SN = BBE_NEGNX16(SN);
        BBE_DSELNX16I(HI, LO, SN, CS, BBE_DSELI_INTERLEAVE_1);

        BBE_SVNX16_IP(LO, pr, (2 * BBE_SIMD_WIDTH));
        BBE_SVNX16_IP(HI, pr, (2 * BBE_SIMD_WIDTH));
    }
    r_align = BBE_ZALIGN();

    /* Since Nrszc is prime, there are no chances to have
    the job complete at the moment, so... aga-ain... */
    {
        PTMP = BBE_RNDADJNX40(P, shr15);
        PTMP = BBE_SRLINX40(PTMP, 15);
        PH = BBE_PACKLNX40(PTMP);

        /*
        exp(-j*pi*q*m.*(m+1)/Nrszc) =
        cos(phases[m]) - j*sin(phases[m]);
        */
        CS = BBE_ZERONX16();
        SN = BBE_ZERONX16();
        t0_ptr = (const xb_vecNx16 *)T + 0;
        t1_ptr = (const xb_vecNx16 *)T + 1;
        t2_ptr = (const xb_vecNx16 *)T + 2;
        t3_ptr = (const xb_vecNx16 *)T + 3;

        TABLE3 = BBE_LVNX16_I(t3_ptr, 0);
        TABLE2 = BBE_LVNX16_I(t2_ptr, 0);
        TABLE1 = BBE_LVNX16_I(t1_ptr, 0);
        TABLE0 = BBE_LVNX16_I(t0_ptr, 0);
        //PH = BBE_POLY8X20_SU(PH, 12, 2);
        // ---- POLY.SU --- //
        x_srl = BBE_SRLINX16(PH, 12);
        sel = xb_vecNx16_rtor_vselN(x_srl);
        //x_srl = x_srl - tmp1 ;
        PH = BBE_POLYNX16_OFF(PH, 12, 2);
        // ---- POLY.SU --- //

        // ---- POLY.STEP --- //
        CS = BBE_SHFLNX16(TABLE3, sel);
        tmp = BBE_SHFLNX16(TABLE2, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        tmp = BBE_SHFLNX16(TABLE1, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        tmp = BBE_SHFLNX16(TABLE0, sel); CS = BBE_MULNX16PACKQ(PH, CS); CS = BBE_ADDNX16(CS, tmp);
        // ---- POLY.STEP --- //
        t0_ptr = (const xb_vecNx16 *)T + 4;
        t1_ptr = (const xb_vecNx16 *)T + 5;
        t2_ptr = (const xb_vecNx16 *)T + 6;
        t3_ptr = (const xb_vecNx16 *)T + 7;

        TABLE3 = BBE_LVNX16_I(t3_ptr, 0);
        TABLE2 = BBE_LVNX16_I(t2_ptr, 0);
        TABLE1 = BBE_LVNX16_I(t1_ptr, 0);
        TABLE0 = BBE_LVNX16_I(t0_ptr, 0);

        // ---- POLY.STEP --- //
        SN = BBE_SHFLNX16(TABLE3, sel);
        tmp = BBE_SHFLNX16(TABLE2, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        tmp = BBE_SHFLNX16(TABLE1, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        tmp = BBE_SHFLNX16(TABLE0, sel); SN = BBE_MULNX16PACKQ(PH, SN); SN = BBE_ADDNX16(SN, tmp);
        // ---- POLY.STEP --- //
        SN = BBE_NEGNX16(SN);
        BBE_DSELNX16I(HI, LO, SN, CS, BBE_DSELI_INTERLEAVE_1);

        if (Nzc == 839)
        {
            BBE_SAVNX16_XP(LO, r_align, pr, 2 * 14);
        }
        else    //case nzc = 139
        {
            BBE_SANX16_IP(LO, r_align, pr);
            BBE_SAVNX16_XP(HI, r_align, pr, 2 * 6);
        }
        BBE_SAPOS_FP(r_align, pr);

    }
} /* zcseq_gen2() */

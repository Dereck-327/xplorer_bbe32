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
    Zadov-Chu sequence generator
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

#if !(HAVE_DIV && 1)
DISCARD_FUN(void, zcseq_gen, (int16_t * restrict r, int u, int v, int M))
#else
/*-------------------------------------------------------------------------
Zadov-Chu sequence generator

Function generates a complex exponential with some properties. Implements
the base algorithm from LTE standard see para 5.5.1.1, 3GPP TS 36.211
V8.8.0 (2009-09)

Accuracy: 9 (2.7e-4)

Input:
  u       Group number, 0..29
  v       Base sequence number within the group, 0..1
  M       Size of sequence (36..1320 in steps of 12)

Output:
  r[2*M]  Output complex Zadoff-Chu sequence, Q15

Return value:
  none

Restrictions:
  r[]     Must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/

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

static const int primes[] =
{
    31, 47, 59, 71, 83, 89, 107, 113,
    131, 139, 151, 167, 179, 191, 199, 211,
    227, 239, 251, 263, 271, 283, 293, 311,
    317, 331, 347, 359, 367, 383, 389, 401,
    419, 431, 443, 449, 467, 479, 491, 503,
    509, 523, 523, 547, 563, 571, 587, 599,
    607, 619, 631, 647, 659, 661, 683, 691,
    701, 719, 727, 743, 751, 761, 773, 787,
    797, 811, 827, 839, 839, 863, 863, 887,
    887, 911, 919, 929, 947, 953, 971, 983,
    991, 997, 1019, 1031, 1039, 1051, 1063, 1069,
    1091, 1103, 1109, 1123, 1129, 1151, 1163, 1171,
    1187, 1193, 1201, 1223, 1231, 1237, 1259, 1259,
    1283, 1291, 1307, 1319
};


// round(2.^(ceil(log2(primes))+15)./primes)
static const uint16_t iprimes[] =
{
    33825, 44620, 35545, 59075, 50534, 47127, 39199, 37118,
    64035, 60350, 55554, 50231, 46864, 43919, 42154, 39756,
    36954, 35099, 33421, 63792, 61909, 59283, 57260, 53946,
    52925, 50686, 48349, 46733, 45714, 43805, 43129, 41838,
    40041, 38926, 37872, 37366, 35926, 35026, 34169, 33354,
    32961, 64158, 64158, 61343, 59599, 58764, 57163, 56017,
    55279, 54207, 53177, 51862, 50917, 50763, 49128, 48559,
    47867, 46668, 46155, 45161, 44680, 44093, 43408, 42636,
    42101, 41374, 40574, 39993, 39993, 38881, 38881, 37829,
    37829, 36833, 36512, 36119, 35432, 35209, 34557, 34135,
    33859, 33655, 32929, 65091, 64590, 63852, 63132, 62777,
    61511, 60842, 60513, 59759, 59441, 58305, 57703, 57309,
    56537, 56252, 55877, 54872, 54516, 54251, 53303, 53303,
    52306, 51982, 51346, 50879 };

// ceil(log2(x))
static const int ipexp[] =
{
    5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 };

static const int16_t ALIGN(32) mm1[] = { 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 16 * 16 }; //m*(m+1)/2, m=1:15 and 16*16

void zcseq_gen ( int16_t * restrict r, int u, int v, int M )
{
    /* Reference MATLAB code
    % Generate special kind of Zadoff-Chu sequence
    % 3GPP TS 36.211 V8.8.0 (2009-09)
    % 5.5.1.1	Base sequences of length   or larger
    % u is in range 0...29
    % v is 0 or 1
    % M is 12*m where m is 3...110
    function z=zc(u,v,M)
    % compute Nrszc as a biggest prime lesser than M
    t=primes(M-1);
    Nrszc=t(end);
    q=Nrszc*(u+1)/31;
    q=floor(q+0.5)+v*((-1)^floor(2*q));
    % q always equal or less than 1. Maximum value 1 is achieved when
    % u=29,v=1,M=36 !!! Minimum value is always positive and never less than
    % 0.0169
    m=(0:Nrszc-1);
    xq=exp(-j*pi*q*m.*(m+1)/Nrszc);
    n=(0:M-1);
    z=xq(mod(n,Nrszc)+1);   % circular copying
    */

    const xb_vecNx16    * restrict T0 = (xb_vecNx16 *)(sincos_table + 0);
    const xb_vecNx16    * restrict T1 = (xb_vecNx16 *)(sincos_table + 1);
    const xb_vecNx16    * restrict T2 = (xb_vecNx16 *)(sincos_table + 2);
    const xb_vecNx16    * restrict T3 = (xb_vecNx16 *)(sincos_table + 3);
    const xb_vecNx16    * restrict T4 = (xb_vecNx16 *)(sincos_table + 4);
    const xb_vecNx16    * restrict T5 = (xb_vecNx16 *)(sincos_table + 5);
    const xb_vecNx16    * restrict T6 = (xb_vecNx16 *)(sincos_table + 6);
    const xb_vecNx16    * restrict T7 = (xb_vecNx16 *)(sincos_table + 7);
    xb_vecNx16    * restrict R = (xb_vecNx16 *)r;
    xb_vecNx16    * restrict Rs = (xb_vecNx16 *)r;
    xb_vecNx16    inrszc, r0, r1, y0, y1;
    xb_vecNx16    ph, pv, plim, pinc, tv, qv;
    xb_vecNx40    A;
    xb_vecNx16    t0, t1, t2, t3, t4, t5, t6, t7;
    int           Nrszc, INrszc;
    int           i, q, m;
    valign        r_align, rs_align;
    vboolN        b0;
    vsaN          enrszc;
    vselN sel;
    xb_vecNx16 lu, hu, ofs, t;

    NASSERT_ALIGN4(r);

    NASSERT(M % 12 == 0);
    NASSERT(M >= 36);

    i = ((M - 36) * 5461 + 32768) >> 16;

    NASSERT(i == (M - 36) / 12);
    ASSERT((((uintptr_t)r) & 31) == 0);

    /* Table indexes are (M-36)/12 */
    Nrszc = primes[i];
    INrszc = iprimes[i];

    /*
    t = Nrszc*(u+1)/31.;
    q = (int)floor(t+.5)+((((int)floor(2*t))&1)?-v:v);
    */

    /* Beware of incorrect estimating of floor(2*t)
    because of integer ariphmetic (e.g., 1.999 instead
    of 2.000 will case an error in sign of v and incorrect
    result for the whole sequence).
    The code below is expected to give exactly the same
    results for u in 0..29, M in 0..1320, but may
    cause errors if these parameters are out of range.
    */
    {
        int32_t t;
        t = Nrszc*(u + 1);      /* 32-bit ariphmetic is enough for M<=1320, u<=29 */
        t *= 33825;           // t in Q.20 
        t += 128;
        t >>= 8;              // rounding to Q.12 
        q = ((int)(t >> 11)) & 1;
        q = (q ? -v : v);
        t += 0x800;           /* rounding to integer */
        t >>= 12;             //Q.0
        q += (int)t;
        q <<= 1;
    }

    enrszc = BBE_MOVVSA32(ipexp[i]);
    inrszc = BBE_MOVVA16(INrszc);
    pinc = BBE_MOVVA16(q);
    plim = BBE_MOVVA16(Nrszc);
    plim = BBE_SLLINX16(plim, 1);

    {
        xb_vecNx40 A;
        xb_vecNx16 _mm1, dvsr, res, zero, _136;
        dvsr = BBE_MOVVA16(2 * Nrszc);
        /* compute pv=mod(q*m*(m+1),2*Nrszc), m=0...15 and qv=mod(q*16*16,2*Nrszc) */
        _mm1 = BBE_LVNX16_I((const xb_vecNx16*)mm1, 0);
        tv = BBE_MOVVA16(q);
        A = BBE_MULNX16(_mm1, tv);
        BBE_DIVNX32S_5STEP0_0(A, dvsr);
        BBE_DIVNX32S_5STEP0_1(A, dvsr);
        BBE_DIVNX16S_4STEP_0(dvsr);
        BBE_DIVNX16S_4STEP_1(dvsr);
        BBE_DIVNX16S_4STEP_0(dvsr);
        BBE_DIVNX16S_4STEP_1(dvsr);
        res = BBE_DIVNX16S_3STEPN_0(dvsr);
        res = BBE_DIVNX16S_3STEPN_1(dvsr);
        pv = BBE_MOVVREM();
        qv = BBE_REPNX16(pv, BBE_SIMD_WIDTH - 1);
        zero = 0;
        pv = BBE_SELNX16I(pv, zero, BBE_SELI_PACK_1);
        /* compute pinc= mod((136+16*m)*q,2*Nrszc), m=0...15 */
        dvsr = BBE_MOVVA16(2 * Nrszc);
        _mm1 = BBE_SEQNX16();
        _mm1 = BBE_SLLINX16(_mm1, 4);
        _136 = BBE_MOVVA16(136);
        _mm1 = BBE_ADDNX16(_136, _mm1);
        A = BBE_MULNX16(_mm1, tv);

        BBE_DIVNX32S_5STEP0_0(A, dvsr);
        BBE_DIVNX32S_5STEP0_1(A, dvsr);
        BBE_DIVNX16S_4STEP_0(dvsr);
        BBE_DIVNX16S_4STEP_1(dvsr);
        BBE_DIVNX16S_4STEP_0(dvsr);
        BBE_DIVNX16S_4STEP_1(dvsr);
        res = BBE_DIVNX16S_3STEPN_0(dvsr);
        res = BBE_DIVNX16S_3STEPN_1(dvsr);
        pinc = BBE_MOVVREM();
    }
    // main loop
    /* Compute phases.
    Simulate the following code:

    for (i=0; i<16; ++i)
    {
    int64_t p = q*(m+i);
    p *= m+i+1;
    p %= 2*Nrszc;
    p <<= 15;
    p += Nrszc/2;
    p /= Nrszc;
    phases[i] = (int16_t)p;
    }
    */
    A = BBE_MULUSRNX16(inrszc, pv, enrszc);
    A = BBE_SRLNX40(A, enrszc);  // NOTE: can not use PACKV here!!!
    ph = BBE_PACKLNX40(A);
    /* Increment p and update pinc */
    pv = BBE_ADDNX16(pv, pinc);
    tv = BBE_SUBNX16(pv, plim);
    b0 = BBE_LENX16(plim, pv);
    pv = BBE_MOVNX16T(tv, pv, b0);
    pinc = BBE_ADDNX16(pinc, qv);
    tv = BBE_SUBNX16(pinc, plim);
    b0 = BBE_LENX16(plim, pinc);
    pinc = BBE_MOVNX16T(tv, pinc, b0);
    for (m = 0; m<Nrszc - (BBE_SIMD_WIDTH - 1); m += BBE_SIMD_WIDTH)
    {
        /*  exp(-j*pi*q*m.*(m+1)/Nrszc) = cos(phases[m]) - j*sin(phases[m]); */
        /*cosine*/
        BBE_LVNX16_IP(t0, T0, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t1, T1, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t2, T2, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t3, T3, 0 * BBE_SIMD_WIDTH);
        /*sine*/
        BBE_LVNX16_IP(t4, T4, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t5, T5, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t6, T6, 0 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t7, T7, 0 * BBE_SIMD_WIDTH);
        /*POLI_SU*/
        ofs = BBE_POLYNX16_OFF(ph, 12, 2);/* Q.15 */
        sel = BBE_MOVVSELNX16(ph, 12);
        /*------------------------------*/
        /*POLI_STEP - 0*/
        lu = BBE_SHFLNX16(t3, sel);
        hu = BBE_SHFLNX16(t7, sel);
        /*POLI_STEP - 1*/
        y0 = BBE_SHFLNX16(t2, sel);
        y1 = BBE_SHFLNX16(t6, sel);
        t = BBE_MULNX16PACKQ(lu, ofs); lu = BBE_ADDNX16(y0, t);
        t = BBE_MULNX16PACKQ(hu, ofs); hu = BBE_ADDNX16(y1, t);
        /*------------------------------*/
        /*POLI_STEP - 2*/
        y0 = BBE_SHFLNX16(t1, sel);
        y1 = BBE_SHFLNX16(t5, sel);
        t = BBE_MULNX16PACKQ(lu, ofs); lu = BBE_ADDNX16(y0, t);
        t = BBE_MULNX16PACKQ(hu, ofs); hu = BBE_ADDNX16(y1, t);
        /*------------------------------*/
        /*POLI_STEP - 3*/
        y0 = BBE_SHFLNX16(t0, sel);
        y1 = BBE_SHFLNX16(t4, sel);
        //Q.15<- Q.15*Q.15 - 15 w/ rounding
        t = BBE_MULNX16PACKQ(lu, ofs); y0 = BBE_ADDNX16(y0, t);
        t = BBE_MULNX16PACKQ(hu, ofs); y1 = BBE_ADDNX16(y1, t);
        /*------------------------------*/
        y1 = BBE_NEGNX16(y1);
        BBE_DSELNX16I(r1, r0, y1, y0, BBE_DSELI_INTERLEAVE_1);
        BBE_SVNX16_IP(r0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1, R, 2 * BBE_SIMD_WIDTH);

        A = BBE_MULUSRNX16(inrszc, pv, enrszc);
        A = BBE_SRLNX40(A, enrszc);  // NOTE: can not use PACKV here!!!
        ph = BBE_PACKLNX40(A);
        /* Increment p and update pinc */
        pv = BBE_ADDNX16(pv, pinc);
        tv = BBE_SUBNX16(pv, plim);
        b0 = BBE_LENX16(plim, pv);
        pv = BBE_MOVNX16T(tv, pv, b0);
        pinc = BBE_ADDNX16(pinc, qv);
        tv = BBE_SUBNX16(pinc, plim);
        b0 = BBE_LENX16(plim, pinc);
        pinc = BBE_MOVNX16T(tv, pinc, b0);
    }
    r_align = BBE_ZALIGN();

    // finalization for remainder < BBE_SIMD_WIDTH/2
    {
        /*
        exp(-j*pi*q*m.*(m+1)/Nrszc) = cos(phases[m]) - j*sin(phases[m]);
        */
        {
            /*  exp(-j*pi*q*m.*(m+1)/Nrszc) = cos(phases[m]) - j*sin(phases[m]); */
            /*cosine*/
            BBE_LVNX16_IP(t0, T0, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t1, T1, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t2, T2, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t3, T3, 0 * BBE_SIMD_WIDTH);
            /*sine*/
            BBE_LVNX16_IP(t4, T4, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t5, T5, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t6, T6, 0 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(t7, T7, 0 * BBE_SIMD_WIDTH);
            /*POLI_SU*/
            ofs = BBE_POLYNX16_OFF(ph, 12, 2);/* Q.15 */
            sel = BBE_MOVVSELNX16(ph, 12);
            /*------------------------------*/
            /*POLI_STEP - 0*/
            lu = BBE_SHFLNX16(t3, sel);
            hu = BBE_SHFLNX16(t7, sel);
            /*POLI_STEP - 1*/
            y0 = BBE_SHFLNX16(t2, sel);
            y1 = BBE_SHFLNX16(t6, sel);
            t = BBE_MULNX16PACKQ(lu, ofs); lu = BBE_ADDNX16(y0, t);
            t = BBE_MULNX16PACKQ(hu, ofs); hu = BBE_ADDNX16(y1, t);
            /*------------------------------*/
            /*POLI_STEP - 2*/
            y0 = BBE_SHFLNX16(t1, sel);
            y1 = BBE_SHFLNX16(t5, sel);
            t = BBE_MULNX16PACKQ(lu, ofs); lu = BBE_ADDNX16(y0, t);
            t = BBE_MULNX16PACKQ(hu, ofs); hu = BBE_ADDNX16(y1, t);
            /*------------------------------*/
            /*POLI_STEP - 3*/
            y0 = BBE_SHFLNX16(t0, sel);
            y1 = BBE_SHFLNX16(t4, sel);
            //Q.15<- Q.15*Q.15 - 15 w/ rounding
            t = BBE_MULNX16PACKQ(lu, ofs); y0 = BBE_ADDNX16(y0, t);
            t = BBE_MULNX16PACKQ(hu, ofs); y1 = BBE_ADDNX16(y1, t);
            /*------------------------------*/
            y1 = BBE_NEGNX16(y1);
            BBE_DSELNX16I(r1, r0, y1, y0, BBE_DSELI_INTERLEAVE_1);
            if ((Nrszc - m)>8)
            {
                BBE_SAVNX16_XP(r0, r_align, R, 2 * BBE_SIMD_WIDTH);
                BBE_SAVNX16_XP(r1, r_align, R, 2 * 2 * (Nrszc - m - 8));
            }
            else
                BBE_SAVNX16_XP(r0, r_align, R, 2 * 2 * (Nrszc - m));
        }
    }
    /* Generate Ouroboros tail */
    m = 2 * (M - Nrszc);
    rs_align = BBE_LAVNX16_PP(Rs);
    while (m>16)
    {
        BBE_LAVNX16_XP(r0, rs_align, Rs, 2 * BBE_SIMD_WIDTH);
        BBE_SAVNX16_XP(r0, r_align, R, 2 * BBE_SIMD_WIDTH);
        m -= BBE_SIMD_WIDTH;
    }
    BBE_LAVNX16_XP(r0, rs_align, Rs, 2 * m);
    BBE_SAVNX16_XP(r0, r_align, R, 2 * m);
    BBE_SAPOS_FP(r_align, R);
} /* zcseq_gen() */
#endif

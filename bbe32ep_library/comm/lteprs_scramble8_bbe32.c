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
    Scramble by PRS
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"
#include "lteprs_common.h"

/*-------------------------------------------------------------------------
Scramble by PRS

Functions generate pseudo-random sequence defined in the LTE standard with
given length and codeword size, when, make bitwise modulo 2 summation
between generated sequence and input codewords. Supported sizes are 1, 2,
4, 6, 8 10 and 16 bits. 
Two versions of routines are available: regular versions (lteprs_scramble1, 
lteprs_scramble2, lteprs_scramble4, lteprs_scramble6, lteprs_scramble8,
lteprs_scramble10,lteprs_scramble16) work with arbitrary arguments, faster versions 
(lteprs_scramble1_fast, lteprs_scramble2_fast, lteprs_scramble4_fast, 
lteprs_scramble6_fast, lteprs_scramble8_fast, lteprs_scramble10_fast,
lteprs_scramble16_fast) 
apply some restrictions.

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Input/Output:
  r[2]   32-bit LFSR states

Output:
  c[N]   Output codewords

Input:
  b[N]   Input codewords

Return value:
  none

Restrictions:
  For fast versions: 
  N        Multiple of 16
  c[],b[]  Aligned by 32-byte boundary
-------------------------------------------------------------------------*/

void lteprs_scramble8 ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N )
{
#if HAVE_LFSR
    xb_vecNx16 * restrict R;
    xb_vecNx16 * restrict C;
    const xb_vecNx16 * restrict P;
    const xb_vecNx16   * restrict B = (const xb_vecNx16 *)b;
    int32_t _c0, _c1;
    valign c_align, r_align, b_align;

    xb_vecNx16 p00, p01, p02, p03, p10, p11, p12, p13;
    xb_vecNx16 c0, c1, b0;
    xb_vecNx16 y0, y;
    xb_vecNx16  mask;
    vselN      sel0;
    vsaN       shft;

    int i;

    P = (const xb_vecNx16*)lteprs_poly;
    R = (xb_vecNx16*)r;
    C = (xb_vecNx16*)c;

    if (N < 0) return;

    c_align = BBE_ZALIGN();
    b_align = BBE_LAVNX16_PP(B);

    c0 = BBE_LPNX16_I(R, 0); c0 = BBE_REPNX16C(c0, 0);
    c1 = BBE_LPNX16_I(R, 4); c1 = BBE_REPNX16C(c1, 0);
    y = BBE_LVNX16_I((const xb_vecNx16*)lteprs8_tbl, 0 * 2 * BBE_SIMD_WIDTH);
    sel0 = BBE_MOVVSELNX16(y, 0);
    shft = BBE_MOVVSA32(8);

    mask = BBE_MOVVA16(255);
    p00 = BBE_LVNX16_I(P, +0 * 2 * BBE_SIMD_WIDTH);
    p01 = BBE_LVNX16_I(P, +1 * 2 * BBE_SIMD_WIDTH);
    p02 = BBE_LVNX16_I(P, +2 * 2 * BBE_SIMD_WIDTH);
    p03 = BBE_LVNX16_I(P, +3 * 2 * BBE_SIMD_WIDTH);

    p10 = BBE_LVNX16_I(P, +4 * 2 * BBE_SIMD_WIDTH);
    p11 = BBE_LVNX16_I(P, +5 * 2 * BBE_SIMD_WIDTH);
    p12 = BBE_LVNX16_I(P, +6 * 2 * BBE_SIMD_WIDTH);
    p13 = BBE_LVNX16_I(P, +7 * 2 * BBE_SIMD_WIDTH);
    BBE_MOVBMULSTATEV(p01, p00, 0);
    BBE_MOVBMULSTATEV(p03, p02, 1);
    _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
    BBE_MOVBMULACCA(_c0);
    // generate by 256 bits
    for (i = 0; i<(N / BBE_SIMD_WIDTH); ++i)
    {
        BMUL32N(c0, 4);
        BBE_MOVBMULSTATEV(p11, p10, 0);
        BBE_MOVBMULSTATEV(p13, p12, 1);
        _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c1);
        BMUL32N(c1, 4);
        y = BBE_XORNX16(c0, c1);
        BBE_LAVNX16_XP(b0, b_align, B, 2 * BBE_SIMD_WIDTH);
        y0 = BBE_SRLNX16(y, shft);
        y0 = BBE_SELNX16(y0, y, sel0);
        y0 = BBE_ANDNX16(y0, mask);
        y0 = BBE_XORNX16(y0, b0);
        BBE_SAVNX16_XP(y0, c_align, C, 2 * BBE_SIMD_WIDTH);
        BBE_MOVBMULSTATEV(p01, p00, 0);
        BBE_MOVBMULSTATEV(p03, p02, 1);
        _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c0);
    }
    N &= 15;
    // tail
    BMUL32N(c0, 4);
    BBE_MOVBMULSTATEV(p11, p10, 0);
    BBE_MOVBMULSTATEV(p13, p12, 1);
    _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
    BBE_MOVBMULACCA(_c1);
    BMUL32N(c1, 4);
    y = BBE_XORNX16(c0, c1);
    BBE_LAVNX16_XP(b0, b_align, B, 2 * N);
    y0 = BBE_SRLNX16(y, shft);
    y0 = BBE_SELNX16(y0, y, sel0);
    y0 = BBE_ANDNX16(y0, mask);
    y0 = BBE_XORNX16(y0, b0);
    BBE_SAVNX16_XP(y0, c_align, C, 2 * N);
    BBE_SAPOS_FP(c_align, C);

    //update LFSR states: shift right by N*8 bits
    N = N * 8 + 96;
    UPDATE_LFSRL(c0, c1, N);
    R = (xb_vecNx16*)r;
    r_align = BBE_ZALIGN();
    BBE_SAVNX16_XP(c0, r_align, R, 2 * 4);
    BBE_SAPOS_FP(r_align, R);
#else
    int i;
    uint32_t rx1, rx2;

    rx1 = r[0];
    rx2 = r[1];

    for (i = 0; i<N; ++i)
    {
        c[i] = (uint16_t)((rx1 & 255) ^ (rx2 & 255) ^ b[i]);
        rx1 = (rx1 >> 8) | ((((rx1 >> 3) & 510) ^ (rx1 & 510)) << 23);
        rx2 = (rx2 >> 8) | ((((rx2 >> 3) & 510) ^ ((rx2 >> 2) & 510) ^ ((rx2 >> 1) & 510) ^ (rx2 & 510)) << 23);
    }

    r[0] = rx1;
    r[1] = rx2;
#endif
} /* lteprs_scramble8() */

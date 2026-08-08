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

// pack 16x1-bit from pos@y to C: first variant!
#define PACK1(C,y,b0,pos,rsh,intlv) \
{                                           \
    xb_vecNx16 three,y0;                    \
    three = BBE_MOVVA16(3);                 \
    y0 = BBE_REPNX16C(y,pos);               \
    y0 = BBE_SRANX16(y0,rsh);               \
    y0 = BBE_ANDNX16 (y0,three);            \
    y0 = BBE_SELNX16(y0,y0,intlv);          \
    y0 = BBE_XORNX16 (y0,b0);               \
    BBE_SVNX16_IP(y0,C,2*BBE_SIMD_WIDTH);   \
}

void lteprs_scramble2_fast ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N )
{
#if HAVE_LFSR
    int32_t _c0, _c1;
    int i;
    int16_t      * restrict Rs  =   (int16_t *)r;
    xb_vecNx16   * restrict C = (xb_vecNx16 *)c;
    const xb_vecNx16   * restrict B = (const xb_vecNx16 *)b;
    xb_vecNx16    p00, p01, p02, p03, p10, p11, p12, p13;
    xb_vecNx16    c0, c1, b0, y;
    vsaN  rsh;
    vselN intlv;
    const xb_vecNx16   * restrict P = (const xb_vecNx16 *)lteprs_poly;

    NASSERT_ALIGN32(c);
    NASSERT_ALIGN32(B);
    NASSERT((N & 15) == 0);
    if (N <= 0) return;

    y = BBE_LVNX16_I((const xb_vecNx16 *)lteprs2_tbl, 0 * BBE_SIMD_WIDTH); rsh = BBE_MOVVSV(y, 0);
    y = BBE_LVNX16_I((const xb_vecNx16 *)lteprs2_tbl, 2 * BBE_SIMD_WIDTH); intlv = BBE_MOVVSELNX16(y, 0);
    c0 = BBE_LPNX16_I(Rs, 0); c0 = BBE_REPNX16C(c0, 0);
    c1 = BBE_LPNX16_I(Rs, 4); c1 = BBE_REPNX16C(c1, 0);

    p00 = BBE_LVNX16_I(P, 0 * BBE_SIMD_WIDTH);
    p01 = BBE_LVNX16_I(P, 2 * BBE_SIMD_WIDTH);
    p02 = BBE_LVNX16_I(P, 4 * BBE_SIMD_WIDTH);
    p03 = BBE_LVNX16_I(P, 6 * BBE_SIMD_WIDTH);
    p10 = BBE_LVNX16_I(P, 8 * BBE_SIMD_WIDTH);
    p11 = BBE_LVNX16_I(P, 10 * BBE_SIMD_WIDTH);
    p12 = BBE_LVNX16_I(P, 12 * BBE_SIMD_WIDTH);
    p13 = BBE_LVNX16_I(P, 14 * BBE_SIMD_WIDTH);

    BBE_MOVBMULSTATEV(p11, p10, 0);
    BBE_MOVBMULSTATEV(p13, p12, 1);
    _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
    BBE_MOVBMULACCA(_c1);
    N >>= 4;
    for (i = 0; i<(N); ++i)
    {
        //Generates 32 bits sequence
        c1 = BBE_BMUL32A(c1, 0);
        BBE_MOVBMULSTATEV(p01, p00, 0);
        BBE_MOVBMULSTATEV(p03, p02, 1);
        _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c0);
        c0 = BBE_BMUL32A(c0, 0);
        y = BBE_XORNX16(c0, c1);
        BBE_LVNX16_IP(b0, B, 2 * BBE_SIMD_WIDTH);
        PACK1(C, y, b0, 6, rsh, intlv);
        BBE_MOVBMULSTATEV(p11, p10, 0);
        BBE_MOVBMULSTATEV(p13, p12, 1);
        _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c1);
    }
    // save generator state (combine from 2 last pairs of c0,c1)
    r[0] = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
    r[1] = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);

#else
    NASSERT_ALIGN32(c);
    NASSERT_ALIGN32(b);
    NASSERT(N % (BBE_SIMD_WIDTH) == 0);
    lteprs_scramble2(r, c, b, N);
#endif
} /* lteprs_scramble2_fast() */

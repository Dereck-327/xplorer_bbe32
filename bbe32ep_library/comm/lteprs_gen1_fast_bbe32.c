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
    Generate LTE PRS sequence
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
Generate LTE PRS sequence 

Functions generate pseudo-random sequence defined in the LTE standard with
given length and codeword size. Supported sizes are 1, 2, 4, 6, 8 10 and 16 
bits.
Two versions of routines are available: regular versions (lteprs_gen1, 
lteprs_gen2, lteprs_gen4, lteprs_gen6, lteprs_gen8, lteprs_gen10, lteprs_gen16) 
work with arbitrary arguments, faster versions (lteprs_gen1_fast, 
lteprs_gen2_fast, lteprs_gen4_fast, lteprs_gen6_fast, lteprs_gen8_fast,
lteprs_gen10_fast, lteprs_gen16_fast) apply some restrictions.

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Input/Output:
  r[2]   32-bit LFSR states

Output:
  c[N]   Generated sequence

Return value:
  none

Restrictions:
  For fast versions:
  N    Multiple of 16
  c[]  Aligned on 32-byte boundary
-------------------------------------------------------------------------*/

// pack 16x1-bit from pos@y to C: first variant!
#define PACK1(C,y,pos,rsh) \
{                                           \
    xb_vecNx16 one,y0;                      \
    one = BBE_MOVVA16(1);                   \
    y0 = BBE_REPNX16(y,pos);                \
    y0 = BBE_SRANX16(y0,rsh);               \
    y0 = BBE_ANDNX16 (y0,one);              \
    BBE_SVNX16_IP(y0,C,2*BBE_SIMD_WIDTH);   \
}

// pack 16x1-bit from pos@y to C: third variant!
#define PACK3(C,y,pos,rsh) \
{                                           \
    xb_vecNx40 onew,yw,y0w;                 \
    xb_vecNx16 y0;                          \
    onew= BBE_MOVWA32(1);                   \
    onew= BBE_SLLNX40(onew,rsh);            \
    yw  = BBE_UNPKUNX16(y);                 \
    y0w= BBE_REPNX40(yw,pos);               \
    y0w= BBE_ANDNX40 (y0w,onew);            \
    y0 = BBE_PACKVNX40(y0w,rsh);            \
    BBE_SVNX16_XP(y0,C,2*BBE_SIMD_WIDTH);   \
}

void lteprs_gen1_fast ( uint32_t * r, uint16_t * restrict c, int N )
{
#if HAVE_LFSR
    int32_t _c1, _c0;
    int i;
    int16_t      * restrict Rs = (int16_t *)r;
    xb_vecNx16   * restrict C = (xb_vecNx16 *)c;
    xb_vecNx16    p00, p01, p02, p03, p10, p11, p12, p13;
    xb_vecNx16    c0, c1, y;
    vsaN rsh;
    const xb_vecNx16   * restrict P = (const xb_vecNx16 *)lteprs_poly;

    NASSERT_ALIGN32(c);
    NASSERT(N % BBE_SIMD_WIDTH == 0);
    if (N <= 0) return;

    y = BBE_SEQNX16(); rsh = BBE_MOVVSV(y, 0);
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
    for (i = 0; i<(N >> 7); ++i)
    {
        //Generates 128 bits sequence
        BMUL32N(c1, 4);
        BBE_MOVBMULSTATEV(p01, p00, 0);
        BBE_MOVBMULSTATEV(p03, p02, 1);
        _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c0);
        BMUL32N(c0, 4);
        y = BBE_XORNX16(c0, c1);

        PACK1(C, y, 6, rsh);
        PACK3(C, y, 7, rsh);
        PACK1(C, y, 8, rsh);
        PACK3(C, y, 9, rsh);
        PACK1(C, y, 10, rsh);
        PACK3(C, y, 11, rsh);
        PACK1(C, y, 12, rsh);
        PACK1(C, y, 13, rsh);
        BBE_MOVBMULSTATEV(p11, p10, 0);
        BBE_MOVBMULSTATEV(p13, p12, 1);
        _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c1);
    }
    // tail
    for (i = 0; i<((N >> 5) & 3); ++i)
    {
        //Generates 32 bits sequence
        c1 = BBE_BMUL32A(c1, 0);
        BBE_MOVBMULSTATEV(p01, p00, 0);
        BBE_MOVBMULSTATEV(p03, p02, 1);
        _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c0);
        c0 = BBE_BMUL32A(c0, 0);
        y = BBE_XORNX16(c0, c1);
        PACK1(C, y, 12, rsh);
        PACK1(C, y, 13, rsh);
        BBE_MOVBMULSTATEV(p11, p10, 0);
        BBE_MOVBMULSTATEV(p13, p12, 1);
        _c1 = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c1);
    }
    // save generator state (combine from 2 last pairs of c0,c1)
    r[0] = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
    r[1] = BBE_EXTRNX16C(c1, BBE_SIMD_WIDTH / 2 - 1);

    if (N & 16)
    {
        c1 = BBE_BMUL32A(c1, 0);
        BBE_MOVBMULSTATEV(p01, p00, 0);
        BBE_MOVBMULSTATEV(p03, p02, 1);
        _c0 = BBE_EXTRNX16C(c0, BBE_SIMD_WIDTH / 2 - 1);
        BBE_MOVBMULACCA(_c0);
        c0 = BBE_BMUL32A(c0, 0);
        y = BBE_XORNX16(c0, c1);
        PACK1(C, y, 12, rsh);

        N = 16 + 64;
        UPDATE_LFSR(c0, c1, N);
        r[0] = BBE_EXTRNX16C(c0, 0);
        r[1] = BBE_EXTRNX16C(c0, 1);
    }
    
#else
    NASSERT_ALIGN32(c);
    NASSERT(N % (BBE_SIMD_WIDTH) == 0);
    lteprs_gen1(r, c, N);
#endif
} /* lteprs_gen1_fast() */

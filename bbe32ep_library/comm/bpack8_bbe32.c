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
    Bit packing
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
Bit packing

Bitstream packing - convert stream of 16-bit words to stream of 1, 2, 3, 4, 6,
8 bits. 

Input:
s[N]      Bitstream formatted as M bits per word
N         Number of symbols
M         Number of bits per symbol
Output:
b[N*M/8]  Unformatted bitstream

Restrictions:
M         1, 2, 3, 4, 6 or 8
N         Multiple of 16
s[]       Aligned on 32-byte boundary
b[]       Aligned on 32-byte boundary
-------------------------------------------------------------------------*/

void bpack8 ( uint8_t * restrict b, const int16_t* s, int N )
{
#if defined(BBE_PACKLNX2X8)

          xb_vecNx16 * restrict B = (xb_vecNx16 *)b;
    const xb_vecNx16 * restrict S = (xb_vecNx16 *)s;
    int i;
    xb_vecNx16 s0, x0, x1;
    valign b_align;
    x1 = BBE_ZERONX16();

    if (N <= 0) return;
    NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    for (i = 0; i<N / (BBE_SIMD_WIDTH * 2); ++i)
    {
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);
        s0 = BBE_PACKLNX2X8(x1, x0);
        BBE_SVNX16_IP(s0, B, 2 * BBE_SIMD_WIDTH);
    }

    if (N & (2 * BBE_SIMD_WIDTH - 1))
    {
        b_align = BBE_ZALIGN();
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        s0 = BBE_PACKLNX2X8(x1, x0);
        BBE_SAVNX16_XP(s0, b_align, B, BBE_SIMD_WIDTH);
        BBE_SAPOS_FP(b_align, B);
    }

#else

          xb_vecNx16 * restrict B = (      xb_vecNx16 *)b;
    const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
    xb_vecNx16 s0, s1, x0, x1, c8, BitMask;
    int i;
    valign b_align;

    if (N <= 0) return;
    NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    c8 = BBE_MOVVA16(1 << 8);
    BitMask = 0xff;

    for (i = 0; i<N / (BBE_SIMD_WIDTH * 2); ++i)
    {
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, S, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_ANDNX16(x0, BitMask);
        x1 = BBE_ANDNX16(x1, BitMask);
        s0 = BBE_SELNX16I(x1, x0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        s1 = BBE_SELNX16I(x1, x0, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
        s1 = BBE_MULNX16PACKL(s1, c8);
        s0 = BBE_ADDNX16(s0, s1);
        BBE_SVNX16_IP(s0, B, 2 * BBE_SIMD_WIDTH);
    }
    if (N & (2 * BBE_SIMD_WIDTH - 1))
    {
        b_align = BBE_ZALIGN();
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_ANDNX16(x0, BitMask);
        s0 = BBE_SELNX16I(x0, x0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        s1 = BBE_SELNX16I(x0, x0, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
        s1 = BBE_MULNX16PACKL(s1, c8);
        s0 = BBE_ADDNX16(s0, s1);
        BBE_SAVNX16_XP(s0, b_align, B, 16);
        BBE_SAPOS_FP(b_align, B);
    }

#endif
} /* bpack8() */

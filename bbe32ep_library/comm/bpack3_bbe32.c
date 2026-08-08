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

void bpack3 ( uint8_t * restrict b, const int16_t* s, int N )
{
          xb_vecNx16 * restrict B = (      xb_vecNx16 *)b;
    const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
    static const union { int16_t i[16]; _vsaN v; } ALIGN(32) lsh[2] =
    {
        { { 16, 16, 16, 16, 12, 16, 16, 16, 8, 16, 16, 16, 4, 16, 16, 16 } },
        { { 0, 16, 16, 16, 4, 16, 16, 16, 8, 16, 16, 16, 12, 16, 16, 16 } }
    };

    int i;
    valign b_align;
    xb_vecNx16 x0, y0, y1, y, c3, c6, BitMask;
    xb_vecNx40 w;
    vselN squeeze;
    static const int32_t isq = 0x111;
    vboolN bsq = BBE_LBN_I((const vboolN*)&isq, 0);

    if (N <= 0) return;
    NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    BBE_SQZN(squeeze, i, bsq);

    b_align = BBE_ZALIGN();
    c3 = BBE_MOVVINT16(1 << 3);
    c6 = BBE_MOVVINT16(1 << 6);
    BitMask = BBE_MOVVINT16(7);

    for (i = 0; i < (N >> LOG2_BBE_SIMD_WIDTH); ++i)
    {
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        // first, group by 12 bits
        //y0 = BBE_SLLINX16(x0,3);
        x0 = BBE_ANDNX16(x0, BitMask);
        y0 = BBE_MULNX16PACKL(x0, c3);
        y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_1);
        x0 = BBE_ADDNX16(x0, y0);
        //y0 = BBE_SLLINX16(x0,6);
        y0 = BBE_MULNX16PACKL(x0, c6);
        y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_2);
        x0 = BBE_ADDNX16(x0, y0);
        // 4x12bits->3x16
        y0 = BBE_SLLNX16(x0, _V(lsh[0].v));
        //y1 = BBE_SRLNX16(x0,lsh[1].v);
        w = BBE_UNPKUNX16(x0);
        y1 = BBE_PACKVNX40(w, _V(lsh[1].v));
        y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_4);    // rotate by 4 elements
        x0 = BBE_ADDNX16(y0, y1);
        // finally, compress 0,4,8
        y = BBE_SHFLNX16(x0, squeeze);
        BBE_SAVNX16_XP(y, b_align, B, 6);
    }
    BBE_SAPOS_FP(b_align, B);
} /* bpack3() */

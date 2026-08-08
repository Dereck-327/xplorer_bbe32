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
    Bit segmentation
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
Bit segmentation

Convert stream of 1, 2, 4, 6, 8 bits to stream of 16-bit words

Input:
b[N*M/8]  Unformatted bitstream
N         Number of symbols
M         Number of bits per symbol
Output:
s[N]      Bitstream formatted as M bits per word

Restrictions:
M   1, 2, 3, 4, 6 or 8
N   Multiple of 16
s   Aligned on 32-byte boundary
b   Aligned on 32-byte boundary
-------------------------------------------------------------------------*/

void bsegmnt3 ( int16_t * restrict s, const uint8_t* b, int N )
{
    const xb_vecNx16 * restrict B = (const xb_vecNx16 *)b;
          xb_vecNx16 * restrict S = (      xb_vecNx16 *)s;
    static const uint16_t ALIGN(32) tbl[3 * 16] =
    {
        2, 5, 8, 11, 14, 1, 4, 7, 10, 13, 0, 3, 6, 9, 12, 15, //rsh0
        16, 16, 16, 16, 16, 15, 16, 16, 16, 16, 14, 16, 16, 16, 16, 16, //rsh1
        7 << 0, 7 << 3, 7 << 6, 7 << 9, 7 << 12, 7 >> 1, 7 << 2, 7 << 5, 7 << 8, 7 << 11, 7 >> 2, 7 << 1, 7 << 4, 7 << 7, 7 << 10, 7 << 13
    };
    int i;
    xb_vecNx16  b0, r0, r1, z0, c2;
    xb_vecNx40 w0, w1;
    valign balign;
    static const int32_t unsq[2] =
    { 0x8210,//0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2
      0x8420 //0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2
    };
    vboolN unsqueeze;
    vselN rep0, rep1;
    vsaN  rsh0, rsh1;

    if (N <= 0) return;
    NASSERT_ALIGN32(s);
    NASSERT_ALIGN32(b);
    NASSERT((N & (BBE_SIMD_WIDTH - 1)) == 0);

    c2 = BBE_MOVVINT16(1 << 2);
    z0 = BBE_LVNX16_I((const xb_vecNx16*)tbl, 0 * 2 * BBE_SIMD_WIDTH); rsh0 = BBE_MOVVSV(z0, 0);
    z0 = BBE_LVNX16_I((const xb_vecNx16*)tbl, 1 * 2 * BBE_SIMD_WIDTH); rsh1 = BBE_MOVVSV(z0, 0);
    z0 = BBE_LVNX16_I((const xb_vecNx16*)tbl, 2 * 2 * BBE_SIMD_WIDTH);
    unsqueeze = BBE_LBN_I((const vboolN*)unsq, 0); BBE_UNSQZN(rep0, i, unsqueeze);
    unsqueeze = BBE_LBN_I((const vboolN*)unsq, 4); BBE_UNSQZN(rep1, i, unsqueeze);
    balign = BBE_LA_PP(B);

    for (i = 0; i<(N >> LOG2_BBE_SIMD_WIDTH); ++i)
    {
        BBE_LAVNX16_XP(b0, balign, B, 6);
        r0 = BBE_SELNX16(b0, b0, rep0);
        r0 = BBE_ANDNX16(r0, z0);
        r1 = BBE_SELNX16(b0, b0, rep1);
        w0 = BBE_MULUUNX16(r0, c2);
        r0 = BBE_PACKVNX40(w0, rsh0);
        w1 = BBE_UNPKUNX16(r1);
        r1 = BBE_PACKVNX40(w1, rsh1);
        r0 = BBE_ADDNX16(r0, r1);
        BBE_SVNX16_IP(r0, S, 2 * BBE_SIMD_WIDTH);
    }
} /* bsegmnt3() */

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

void bsegmnt2 ( int16_t * restrict s, const uint8_t* b, int N )
{
    const xb_vecNx16 * restrict B = (const xb_vecNx16 *)b;
          xb_vecNx16 * restrict S = (      xb_vecNx16 *)s;
    static const uint16_t ALIGN(32) table[16] = { 0x4000, 0x1000, 0x400, 0x100, 0x40, 0x10, 0x4, 0x1, 0x4000, 0x1000, 0x400, 0x100, 0x40, 0x10, 0x4, 0x1 };
    int i;
    const xb_vecNx16 * restrict R = (xb_vecNx16 *)table;
    xb_vecNx16  b0, r0, z1;
    static const int32_t unsq = 0x8080; // corresponds to 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1
    const vboolN unsqueeze = BBE_LBN_I((const vboolN*)&unsq, 0);
    vselN sel;
    valign balign;
    
    if (N <= 0) return;
    NASSERT_ALIGN32(s);
    NASSERT_ALIGN32(b);
    NASSERT((N & (BBE_SIMD_WIDTH - 1)) == 0);

    BBE_UNSQZN(sel, i, unsqueeze);
    z1 = BBE_LVNX16_I(R, 0);
    balign = BBE_LA_PP(B);

    for (i = 0; i<(N >> LOG2_BBE_SIMD_WIDTH); ++i)
    {
        BBE_LAVNX16_XP(b0, balign, B, 4);
        b0 = BBE_SELNX16(b0, b0, sel);
        r0 = BBE_MULNX16PACKL(b0, z1);
        r0 = BBE_SRLINX16(r0, 14);
        BBE_SVNX16_IP(r0, S, 2 * BBE_SIMD_WIDTH);
    }
} /* bsegmnt2() */

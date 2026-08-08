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

void bpack1 ( uint8_t * restrict b, const int16_t* s, int N )
{
          xb_int16   * restrict B = (      xb_int16   *)b;
    const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
    int i;
    xb_vecNx40 w;
    xb_vecNx16 c, s0;
    xb_int16 z0, _1;

    if (N <= 0) return;
    NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    _1 = BBE_MOVVINT16(1);
    c = BBE_SLLNX16(_1, BBE_MOVVSV(BBE_SEQNX16(), 0)); // form mask 1,2,4....0x8000

    __Pragma("loop_count min=1");
    for (i = 0; i < (N >> LOG2_BBE_SIMD_WIDTH); ++i)
    {
        BBE_LVNX16_XP(s0, S, 2 * BBE_SIMD_WIDTH);
        s0 = BBE_ANDNX16(s0, _1);
        w = BBE_MULUUNX16(s0, c);
        s0 = BBE_PACKLNX40(w);
        z0 = BBE_RADDNX16(s0);
        xb_int16_storeip(z0, B, 2);
    }
} /* bpack1() */

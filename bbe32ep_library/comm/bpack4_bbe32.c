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

void bpack4 ( uint8_t * restrict b, const int16_t* s, int N )
{
        xb_vecNx16 * restrict B = (      xb_vecNx16 *)b;
  const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
  int i;
  xb_vecNx40 w;
  xb_vecNx16 x0, x1, x2, x3, z0, z1, z2, z3;
  xb_vecNx16 c1, c2, c3, BitMask;
  valign b_align, s_align;

  if (N <= 0) return;
  NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
  NASSERT(N % BBE_SIMD_WIDTH == 0);

  c1 = BBE_MOVVINT16(0x10);
  c2 = BBE_MOVVA16(0x100);
  c3 = BBE_MOVVA16(0x1000);
  b_align = BBE_ZALIGN();
  BitMask = BBE_MOVVINT16(15);

  __Pragma("ymemory( S )");
  for (i = 0; i<N / (BBE_SIMD_WIDTH * 4); ++i)
  {
    BBE_LVNX16_IP(x0, S, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x1, S, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x2, S, (2 * BBE_SIMD_WIDTH));
    BBE_LVNX16_IP(x3, S, (2 * BBE_SIMD_WIDTH));
    x0 = BBE_ANDNX16(x0, BitMask);
    x1 = BBE_ANDNX16(x1, BitMask);
    x2 = BBE_ANDNX16(x2, BitMask);
    x3 = BBE_ANDNX16(x3, BitMask);

    BBE_DSELNX16I(z1, z0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(z3, z2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x1, x0, z2, z0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x3, x2, z3, z1, BBE_DSELI_DEINTERLEAVE_1);

    w = BBE_UNPKUNX16(x0);
    BBE_MULUUANX16(w, x1, c1);
    BBE_MULUUANX16(w, x2, c2);
    BBE_MULUUANX16(w, x3, c3);
    x0 = BBE_PACKLNX40(w);
    BBE_SAVNX16_XP(x0, b_align, B, (2 * BBE_SIMD_WIDTH));
  }
  N &= (4 * BBE_SIMD_WIDTH - 1);
  if (N)
  {
    s_align = BBE_LANX16_PP(S);
    BBE_LAVNX16_XP(x0, s_align, S, 2 * (N - 0 * 16));
    BBE_LAVNX16_XP(x1, s_align, S, 2 * (N - 1 * 16));
    BBE_LAVNX16_XP(x2, s_align, S, 2 * (N - 2 * 16));
    BBE_LAVNX16_XP(x3, s_align, S, 2 * (N - 3 * 16));
    x0 = BBE_ANDNX16(x0, BitMask);
    x1 = BBE_ANDNX16(x1, BitMask);
    x2 = BBE_ANDNX16(x2, BitMask);
    x3 = BBE_ANDNX16(x3, BitMask);

    BBE_DSELNX16I(z1, z0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(z3, z2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(x1, x0, z2, z0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(x3, x2, z3, z1, BBE_DSELI_DEINTERLEAVE_1);

    w = BBE_UNPKUNX16(x0);
    BBE_MULUUANX16(w, x1, c1);
    BBE_MULUUANX16(w, x2, c2);
    BBE_MULUUANX16(w, x3, c3);
    x0 = BBE_PACKLNX40(w);
    BBE_SAVNX16_XP(x0, b_align, B, (N >> 1));
  }
  BBE_SAPOS_FP(b_align, B);
} /* bpack4() */

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

void bpack6 ( uint8_t * restrict b, const int16_t* s, int N )
{
        xb_vecNx16 * restrict B = (      xb_vecNx16 *)b;
  const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
  static const int16_t ALIGN(32) tbl[2 * BBE_SIMD_WIDTH] =
  {
      0, 0, 1 << 12, 0, 1 << 8, 0, 1 << 4, 0, 0, 0, 1 << 12, 0, 1 << 8, 0, 1 << 4, 0,
      0, 16, 4, 16, 8, 16, 12, 16, 0, 16, 4, 16, 8, 16, 12, 16
  };
  int i;
  valign b_align;
  xb_vecNx16  x0, y0, y1, y, c6, c0, BitMask;
  vselN squeeze;
  vsaN shr;
  static const int32_t isq = 0x1515;
  vboolN bsq = BBE_LBN_I((const vboolN*)&isq, 0);

  if (N <= 0) return;
  NASSERT_ALIGN(b, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
  NASSERT(N % BBE_SIMD_WIDTH == 0);

  b_align = BBE_ZALIGN();
  c6 = BBE_MOVVINT16(1 << 6);
  c0 = BBE_LVNX16_I((const xb_vecNx16*)tbl, 0);
  y0 = BBE_LVNX16_I((const xb_vecNx16*)tbl, 1 * 2 * BBE_SIMD_WIDTH);
  shr = BBE_MOVVSV(y0, 0);
  BBE_SQZN(squeeze, i, bsq);
  BitMask = BBE_MOVVINT16(63);
  __Pragma("loop_count min=1");
  for (i = 0; i < (N >> LOG2_BBE_SIMD_WIDTH); ++i)
  {
    BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
    x0 = BBE_ANDNX16(x0, BitMask);
    // first group by 8x12bits
    y0 = BBE_MULNX16PACKL(x0, c6);
    y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_1);
    x0 = BBE_ADDNX16(y0, x0);
    // next, combine 8x12->6x16
    y0 = BBE_MULNX16PACKL(x0, c0);
#if 1
    {
      xb_vecNx40  w;
      w = BBE_UNPKUNX16(x0);
      y1 = BBE_PACKVNX40(w, shr);
    }
#else
    y1 = BBE_SRLNX16(x0, shr);
#endif
    y0 = BBE_SELNX16I(y0, y0, BBE_SELI_ROTATE_RIGHT_2);
    x0 = BBE_ADDNX16(y0, y1);
    // compress 0,2,4,6,10,12,  ->  0,1,2,3,4,5
    y = BBE_SELNX16(x0, x0, squeeze);
    BBE_SAVNX16_XP(y, b_align, B, 12);
  }
  BBE_SAPOS_FP(b_align, B);
} /* bpack6() */

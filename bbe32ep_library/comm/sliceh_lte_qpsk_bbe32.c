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
    Hard demapper
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
Hard demapper

Hard demapper algorithm converts equalized complex data into the sequence
decoded bits. It supports QAM/QPSK constellations from LTE standards.

For detailed information on mapping bits onto constellations, please refer
to
[1] ETSI TS 136 211 V10.0.0 (2011-01); LTE; Evolved Universal Terrestrial
    Radio Access (E-UTRA); Physical channels and modulation (3GPP TS
    36.211 version 10.0.0 Release 10);

Methods:
  sliceh_lte_qpsk()       - LTE QPSK [1] 7.1.2
  sliceh_lte_16qam()      - LTE 16QAM [1] 7.1.3
  sliceh_lte_64qam()      - LTE 64QAM [1] 7.1.4
Input:
   s[2*N]  Input complex samples. I and Q components are interleaved
           with the real part going first (at even indices). Fixed point
           position depends on the constellation size, as described
           above
Output:
   b[N]    Bitstream formatted as sequence of 16-bit words each containing 
           K-tuples
Restrictions:
  s[]      Aligned on 32-byte boundary
  b        Aligned on 32-byte boundary
  N        Multiple of 8
Performance restrictions:
  None
---------------------------------------------------------------------------*/

/* LTE QPSK [1] 7.1.2 */
void sliceh_lte_qpsk ( int16_t * restrict b,
                  const int16_t *          s,
                  int N)
{
    int n;
    xb_vecNx16 X0, X1;
    valign balign = BBE_ZALIGN();
          xb_vecNx16 * restrict B = (      xb_vecNx16*)b;
    const xb_vecNx16 * restrict S = (const xb_vecNx16*)s;
   
    if (N <= 0) return;
    NASSERT_ALIGN(b, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(s, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
    
    for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH)); n++)
    {
        BBE_LVNX16_IP(X0, S, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(X1, S, 2 * BBE_SIMD_WIDTH);
        BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_1);
        X0 = BBE_SRLINX16(X0, 15);
        X1 = BBE_SRLINX16(X1, 15);
        X1 = BBE_ADDNX16(X1, X1);
        X0 = BBE_ADDNX16(X1, X0);
        BBE_SVNX16_XP(X0, B, 2 * BBE_SIMD_WIDTH);
    }
    if (N&(BBE_SIMD_WIDTH / 2))
    {
        BBE_LVNX16_IP(X0, S, 2 * BBE_SIMD_WIDTH);
        BBE_DSELNX16I(X1, X0, X0, X0, BBE_DSELI_DEINTERLEAVE_1);
        X0 = BBE_SRLINX16(X0, 15);
        X1 = BBE_SRLINX16(X1, 15);
        X1 = BBE_ADDNX16(X1, X1);
        X0 = BBE_ADDNX16(X1, X0);
        BBE_SAVNX16_XP(X0, balign, B, BBE_SIMD_WIDTH);
        BBE_SAPOS_FP(balign, B);
    }
} /* sliceh_lte_qpsk() */

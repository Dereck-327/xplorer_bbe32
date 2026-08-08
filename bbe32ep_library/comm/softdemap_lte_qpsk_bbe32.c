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
    Soft demapper
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
Soft demapper

  Soft demapper algorithm converts equalized complex data into the sequence
  of log-likelihood ratios (LLRs). It supports QAM/QPSK constellations from
  LTE standards.

  For detailed information on mapping bits onto constellations, please refer
  to
  [1] ETSI TS 136 211 V10.0.0 (2011-01); LTE; Evolved Universal Terrestrial
      Radio Access (E-UTRA); Physical channels and modulation (3GPP TS
      36.211 version 10.0.0 Release 10);

  Soft decoding algorithm is run over a block of I/Q samples. For each input
  sample s, the demapper computes approximations of sub-optimal simplified
  log-likelihood ratios (LLR) for K bit positions, where K equals the base-2 
  logarithm of the constellation size. Let S(0,k,s) denote the constellation
  point with 0 at k-th bit position which is closest to s, while S(1,k,s) is
  the closest point with 1 at k-th bit position, k=0..K-1. Then the value of
  sub-optimal simplified LLR is defined as:

    LLR(s,k) = (|s-S(0,k,s)|^2 - |s-S(1,k,s)|^2),

          |...|^2 - squared absolute value of a complex number

  I/Q components are represented in Q10

  Methods:
    softdemap_lte_qpsk()       - LTE QPSK [1] 7.1.2
    softdemap_lte_16qam()      - LTE 16QAM [1] 7.1.3
    softdemap_lte_64qam()      - LTE 64QAM [1] 7.1.4
  Input:
     s[2*N]     Input complex samples. I and Q components are interleaved
                with the real part going first (at even indices). Fixed point
                position depends on the constellation size, as described
                above
  Output:
     llr[K*N]   LLR estimations, Q0. K equals the base-2 logarithm of the
                constellation size. llr[0] holds the LLR value for the left-most
                bit position of a K-tuple that matches s[0]+j*s[1].
  Restrictions:
    s[], llr[]  Must be aligned on 32-byte boundary
    N           Must be a multiple of 8
  Performance restrictions:
    None
---------------------------------------------------------------------------*/

/* LTE QPSK [2] 7.1.2 */
void softdemap_lte_qpsk ( int16_t * restrict llr,
                    const int16_t *          s,
                    int N )
{
#if ( HAVE_SDMAP && 1 )
          xb_vecNx16 * restrict L = (xb_vecNx16 *)llr;
    const xb_vecNx16 *          S = (xb_vecNx16 *)s;
    xb_vecNx16 x, y0, y1, z, scl;
    int n;

    NASSERT(llr != NULL && s != NULL);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(llr, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

    scl = BBE_MOVVA16C(0x00000008); // -1.0

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        /* Load 8 CQ5.10 samples */
        BBE_LVNX16_IP(x, S, 2 * BBE_SIMD_WIDTH);
        /* Phase 0, demap lower BBE_SIMD_WIDTH/4 complex inputs to BBE_SIMD_WIDTH/2 soft bits */
        y0 = BBE_SDMAPQPSKNX16C(scl, x, 0, 1);
        /* Phase 1, demap higher BBE_SIMD_WIDTH/4 complex inputs to BBE_SIMD_WIDTH/2 soft bits */
        y1 = BBE_SDMAPQPSKNX16C(scl, x, 1, 1);
        /* Concatenate soft bits from phase 0 to soft bits from phase 1. */
        z = BBE_SELNX16I(y1, y0, BBE_SELI_INTERLEAVE_4_LO);
        /* Convert signed 8-bit values to 16-bit. */
        z = BBE_UNPKSNX2X8L(z);
        /* Save BBE_SIMD_WIDTH soft bits */
        BBE_SVNX16_IP(z, L, BBE_SIMD_WIDTH * 2);
    }
#else
          xb_vecNx16 * restrict L = (      xb_vecNx16 *)llr;
    const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
    int n;
    xb_vecNx16 _1_256, x0;
    xb_vecNx16 _8, __8;

    NASSERT(llr != NULL && s != NULL);
    NASSERT_ALIGN32(s);
    NASSERT_ALIGN32(llr);
    NASSERT((N % (BBE_SIMD_WIDTH / 2) == 0));

    _1_256 = BBE_MOVVA16(-128);   // 1/256,Q15
    _8 = BBE_MOVVINT16(8);
    __8 = BBE_NEGNX16(_8);

    // Apply the recursive form of the Tosato-Bisaglia algorithm.
    for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
    {
        // Load 8 complex samples
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_MULNX16PACKQ(x0, _1_256);
        x0 = BBE_MAXNX16(__8, x0);
        x0 = BBE_MINNX16(_8, x0);
        BBE_SVNX16_IP(x0, L, 2 * BBE_SIMD_WIDTH);
    }
#endif
} /* softdemap_lte_qpsk() */

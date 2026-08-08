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

/* LTE 64QAM [2] 7.1.4 */
void softdemap_lte_64qam ( int16_t * restrict llr,
                     const int16_t *          s,
                     int N )
{
#if ( HAVE_SDMAP && 1 )
          xb_vecNx16 * restrict L = (      xb_vecNx16 *)llr;
    const xb_vecNx16 *          S = (const xb_vecNx16 *)s;
    xb_vecNx16 x, y0, y1, z0, z1, z2, scl;
    vsaN sel;
    int n;

    static const ALIGN(32) int16_t sel_i[BBE_SIMD_WIDTH] = { 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };

    NASSERT(llr != NULL && s != NULL);
    NASSERT_ALIGN(s, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(llr, 2 * BBE_SIMD_WIDTH);
    NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

    scl = BBE_MOVVA16C(0x00000008); // -1.0
    sel = BBE_MOVVSELNX16(BBE_LVNX16_I((xb_vecNx16*)sel_i, 0), 0);

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        /* Load 8 CQ5.10 samples */
        BBE_LVNX16_IP(x, S, BBE_SIMD_WIDTH * 2);
        /* Phase 0, demap lower BBE_SIMD_WIDTH/4 complex inputs to 6*BBE_SIMD_WIDTH/4 soft bits */
        y0 = BBE_SDMAP64QAMNX16C(scl, x, 0, 1, 1);
        /* Phase 1, demap higher BBE_SIMD_WIDTH/4 complex inputs to 6*BBE_SIMD_WIDTH/4 soft bits */
        y1 = BBE_SDMAP64QAMNX16C(scl, x, 1, 1, 1);
        /* Remove the hole between phase 0/1 soft bits. */
        y1 = BBE_SELNX16(y1, y0, sel);
        /* Convert signed 8-bit values to 16-bit. */
        z0 = BBE_UNPKSNX2X8L(y0);
        z1 = BBE_UNPKSNX2X8L(y1);
        z2 = BBE_UNPKSNX2X8H(y1);
        /* Save 3*BBE_SIMD_WIDTH soft bits */
        BBE_SVNX16_IP(z0, L, BBE_SIMD_WIDTH * 2);
        BBE_SVNX16_IP(z1, L, BBE_SIMD_WIDTH * 2);
        BBE_SVNX16_IP(z2, L, BBE_SIMD_WIDTH * 2);
    }
#else
          xb_vecNx16 * restrict L = (      xb_vecNx16 *)llr;
    const xb_vecNx16 * restrict S = (const xb_vecNx16 *)s;
    int n;
    xb_vecNx16    x1;
    xb_vecNx16    x0;
    xb_vecNx16 Y1, Y0;
    xb_vecNx16    l1, l0, l2;
    xb_vecNx16    _1_256;
    xb_vecNx16    _4Q10, _2Q10, _min, _max;

    NASSERT(llr != NULL && s != NULL);
    NASSERT_ALIGN32(s);
    NASSERT_ALIGN32(llr);
    NASSERT((N % (BBE_SIMD_WIDTH / 2)) == 0);

    _1_256 = BBE_MOVVA16(-128);   // 1/256,Q15
    _2Q10 = BBE_MOVPINT16(2);   // 2,Q10
    _4Q10 = BBE_MOVPINT16(4);   // 4,Q10
    _max = BBE_MOVPINT16(8);
    _min = BBE_NEGNX16(_max);

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        // Load 8 complex samples
        BBE_LVNX16_IP(x0, S, 2 * BBE_SIMD_WIDTH);
        x0 = BBE_MINNX16(_max, x0);
        x0 = BBE_MAXNX16(_min, x0);
        x1 = BBE_ABSSNX16(x0);
        l1 = BBE_SUBSNX16(_4Q10, x1);
        l2 = BBE_ABSSNX16(l1);
        l2 = BBE_SUBSNX16(_2Q10, l2);
        // multiplex 3 outputs
        BBE_DSELNX16I(Y1, Y0, l1, x0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(l1, l0, l2, Y0, BBE_DSELI_INTERLEAVE_C3_STEP_0);
        BBE_DSELNX16I_H(l1, l2, l2, Y1, BBE_DSELI_INTERLEAVE_C3_STEP_1);

        l0 = BBE_MULNX16PACKQ(l0, _1_256);
        l1 = BBE_MULNX16PACKQ(l1, _1_256);
        l2 = BBE_MULNX16PACKQ(l2, _1_256);
        BBE_SVNX16_IP(l0, L, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(l1, L, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(l2, L, 2 * BBE_SIMD_WIDTH);
    }
#endif
} /* softdemap_lte_64qam() */

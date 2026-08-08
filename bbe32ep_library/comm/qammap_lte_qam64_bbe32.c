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
    LTE QAM64 mapper
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
  LTE QAM64 mapper

  Translate binary 6-tuples into LTE QAM64 points.

  Input:
    d[N]     Binary 6-tuples, each one stored at the least significant bits
             of a 16-bit word (0..5). 10 MSBs of the word (6..15) must hold
             zero
  Output:
   r[2*N]    QAM64 points, Q5.10. In-phase (I) and quadrature (Q) phase 
             components are interleaved with the I-component going first
             (at even indices).
  Restrictions:
    N        Multiple of 16
    d[],r[]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/

void qammap_lte_qam64 ( int16_t * restrict r,
                  const int16_t *          d,
                  int N )
{
    /*    1/4 LTE QAM64 constellation (I and Q components ): */
    static const int16_t ALIGN(32) qam64_mapIQ[32] =
    {                             //type C    
        4096 / 4 * 3, 4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 1, 4096 / 4 * 3, 4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 1, 4096 / 4 * 5, 4096 / 4 * 5, 4096 / 4 * 7, 4096 / 4 * 7, 4096 / 4 * 5, 4096 / 4 * 5, 4096 / 4 * 7, 4096 / 4 * 7,
        4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 5, 4096 / 4 * 7, 4096 / 4 * 5, 4096 / 4 * 7, 4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 3, 4096 / 4 * 1, 4096 / 4 * 5, 4096 / 4 * 7, 4096 / 4 * 5, 4096 / 4 * 7
    };
    int n;
    const xb_vecNx16 * restrict D = (const xb_vecNx16 *)d;
    xb_vecNx16 * restrict R = (xb_vecNx16 *)r;
    xb_vecNx16  t_re, t_im;
    xb_vecNx16  c10, c11, d0, d1, r0, r1, s0, s1, mask15, zero;
    vselN       idx;
    vboolN b0, b1;

    if (N <= 0) return;
    NASSERT(r != NULL && d != NULL && N >= 0);
    NASSERT_ALIGN32(r);
    NASSERT_ALIGN32(d);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    mask15 = BBE_MOVVA16(15);
    c10 = BBE_MOVVA16(1 << 10);
    c11 = BBE_MOVVA16(1 << 11);
    zero = BBE_ZERONX16();
    t_re = BBE_LVNX16_I((const xb_vecNx16*)qam64_mapIQ, 0);
    t_im = BBE_LVNX16_I((const xb_vecNx16*)qam64_mapIQ, (2 * BBE_SIMD_WIDTH));

    for (n = 0; n<N / (BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(d0, D, 2 * BBE_SIMD_WIDTH);
        idx = BBE_MOVVSELNX16(d0, 0);
        r0 = BBE_SELNX16(t_re, t_re, idx);
        r1 = BBE_SELNX16(t_im, t_im, idx);
        s0 = BBE_MULNX16PACKL(d0, c10);
        s1 = BBE_MULNX16PACKL(d0, c11);
        b0 = BBE_LTNX16(s0, zero);
        b1 = BBE_LTNX16(s1, zero);
        BBE_NEGNX16T(r0, r0, b0);
        BBE_NEGNX16T(r1, r1, b1);
        BBE_DSELNX16I(d1, d0, r1, r0, BBE_DSELI_INTERLEAVE_1);
        BBE_SVNX16_IP(d0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(d1, R, 2 * BBE_SIMD_WIDTH);
    }
} /* qammap_lte_qam64() */

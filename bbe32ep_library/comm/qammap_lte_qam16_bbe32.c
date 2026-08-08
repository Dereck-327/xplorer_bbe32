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
    Fast LTE QAM16 mapper
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
  Fast LTE QAM16 mapper

  Translate binary 4-tuples into LTE QAM16 points.

  Input:
    d[N]     Binary 4-tuples, each one stored at the least significant bits
             of a 16-bit word (0..3). 12 MSBs of the word (4..15) must hold
             zero
  Output:
   r[2*N]    QAM16 points, Q5.10. In-phase (I) and quadrature (Q) phase 
             components are interleaved with the I-component going first
             (at even indices).
  Restrictions:
    N        Multiple of 16
    d[],r[]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/

void qammap_lte_qam16 ( int16_t * restrict r,
                  const int16_t *          d,
                  int N )
{
    /*
    LTE QAM16 constellation (I and Q components ):
    0   4096   4096
    1   4096  12288
    2  12288   4096
    3  12288  12288
    4   4096  -4096
    5   4096 -12288
    6  12288  -4096
    7  12288 -12288
    8  -4096   4096
    9  -4096  12288
    10 -12288   4096
    11 -12288  12288
    12  -4096  -4096
    13  -4096 -12288
    14 -12288  -4096
    15 -12288 -12288
    */
    static const int16_t ALIGN(32) qam16_mapIQ[32] =
    {                             //type C    
        4096 / 4, 4096 / 4, 12288 / 4, 12288 / 4, 4096 / 4, 4096 / 4, 12288 / 4, 12288 / 4, -4096 / 4, -4096 / 4, -12288 / 4, -12288 / 4, -4096 / 4, -4096 / 4, -12288 / 4, -12288 / 4,
        4096 / 4, 12288 / 4, 4096 / 4, 12288 / 4, -4096 / 4, -12288 / 4, -4096 / 4, -12288 / 4, 4096 / 4, 12288 / 4, 4096 / 4, 12288 / 4, -4096 / 4, -12288 / 4, -4096 / 4, -12288 / 4
    };
    int n;
    const xb_vecNx16 * restrict D = (const xb_vecNx16 *)d;
          xb_vecNx16 * restrict R = (      xb_vecNx16 *)r;
    xb_vecNx16  t_re, t_im;
    xb_vecNx16  d0, d1, r0, r1;
    vselN       idx;

    if (N <= 0) return;
    NASSERT(r != NULL && d != NULL && N >= 0);
    NASSERT_ALIGN32(r);
    NASSERT_ALIGN32(d);
    NASSERT(N % BBE_SIMD_WIDTH == 0);

    t_re = BBE_LVNX16_I((const xb_vecNx16*)qam16_mapIQ, 0);
    t_im = BBE_LVNX16_I((const xb_vecNx16*)qam16_mapIQ, (2 * BBE_SIMD_WIDTH));

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(d0, D, 2 * BBE_SIMD_WIDTH);
        idx = BBE_MOVVSELNX16(d0, 0);
        r0 = BBE_SELNX16(t_re, t_re, idx);
        r1 = BBE_SELNX16(t_im, t_im, idx);
        BBE_DSELNX16I(d1, d0, r1, r0, BBE_DSELI_INTERLEAVE_1);
        BBE_SVNX16_IP(d0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(d1, R, 2 * BBE_SIMD_WIDTH);
    }
} /* qammap_lte_qam16() */

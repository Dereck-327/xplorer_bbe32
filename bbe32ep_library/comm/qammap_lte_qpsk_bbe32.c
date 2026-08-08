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
    LTE QPSK mapper
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
  LTE QPSK mapper

  Translate binary 2-tuples into QPSK points. QPSK LTE constellation: 
                          
                          [10]    [00]
                          -1+j     1+j
                     
                          -1-j     1-j
                          [11]    [01]
  Input:
    d[N]         Binary 2-tuples, each one stored at the least significant bits
                 of a 16-bit word (0..1). 14 MSBs of the word (2..15) must hold
                 zero
  Output:
    r[2*N]       QPSK points, Q5.10. In-phase (I) and quadrature (Q) phase 
                 components are interleaved with the I-component going first
                 (at even indices).
  Restrictions:
    N            Multiple of 16
    d[N],r[2*N]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/

/* Conversion of a floating-point complex number into Q3.12 format. */
#define CQ10( re, im )      (int16_t)( (re)*(1<<10) ), \
                            (int16_t)( (im)*(1<<10) )

void qammap_lte_qpsk ( int16_t * restrict r,
                 const int16_t *          d,
                 int N )
{
    // QPSK type C constellation (LTE, see [2] 7.1.2)
    static const int16_t ALIGN(16) qpsk_c_map[2 * 4] =
    { //      00               01               10               11
        CQ10(+1.0, +1.0), CQ10(+1.0, -1.0), CQ10(-1.0, +1.0), CQ10(-1.0, -1.0)
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

    d0 = BBE_LVNX16_I((const xb_vecNx16*)qpsk_c_map, 0);
    t_re = BBE_SELNX16I(d0, d0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    t_im = BBE_SELNX16I(d0, d0, BBE_SELI_EXTRACT_1_OF_2_OFF_1);

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(d0, D, 2 * BBE_SIMD_WIDTH);
        d0 = BBE_ANDNX16(d0, 3);
        idx = BBE_MOVVSELNX16(d0, 0);
        r0 = BBE_SELNX16(t_re, t_re, idx);
        r1 = BBE_SELNX16(t_im, t_im, idx);
        BBE_DSELNX16I(d1, d0, r1, r0, BBE_DSELI_INTERLEAVE_1);
        BBE_SVNX16_IP(d0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(d1, R, 2 * BBE_SIMD_WIDTH);
    }
} /* qammap_lte_qpsk() */

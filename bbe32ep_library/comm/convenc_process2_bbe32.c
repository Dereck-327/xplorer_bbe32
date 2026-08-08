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
    Encoding
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

#if !(HAVE_INTLV && HAVE_LFSR && 1)
DISCARD_FUN(int16_t, convenc_process2, (convenc_handle_t handle, 
                         int16_t* e, const int16_t* s, int N))
#else
/*---------------------------------------------------------------------------
Encoding

Input:
s[N]            Unformatted bitstream (16 bits per word)
N               Number of input 16-bit words
Output:
e[N*R]          Unformatted bitstream
Returned value  Updated state of encoder
---------------------------------------------------------------------------*/

#define MAGIC 0x756bb623
typedef struct
{
    uint32_t magic;       // Instance pointer validation number
    uint16_t  polyMask[4];// polynomials 
    int16_t  R;           // code rate (2…4)
    int16_t  K;           // constraint length (3…9)
    int16_t  state;       // state of encoder
}
tConvenc_;

int16_t convenc_process2 ( convenc_handle_t handle, 
                           int16_t * restrict e, const int16_t * s, int N )
{
    tConvenc_ *pEnc = (tConvenc_ *)handle;
    const xb_vecNx16 * restrict S = (xb_vecNx16 *)s;
    xb_vecNx16 * restrict E = (xb_vecNx16 *)e;
    int n, M, K;
    uint16_t p0, p1;
    valign s_align, e_align;
    xb_vecNx16  s0, s1, st, y0, y1, s3;
    vsaN        sh15;
    NASSERT(pEnc != NULL);
    NASSERT(pEnc->R == 2);
    NASSERT(pEnc->magic == MAGIC);
    if (N <= 0) return pEnc->state;
    p0 = pEnc->polyMask[0] * 2;
    p1 = pEnc->polyMask[1] * 2;
    K = pEnc->K;
    st = BBE_MOVVA16(pEnc->state << (16 - K));
    s_align = BBE_LAVNX16_PP(S);
    e_align = BBE_ZALIGN();
    sh15 = BBE_MOVVSA32(15);
    y0 = y1 = BBE_ZERONX16();
    M = ((int16_t)N * 34953) >> 19; // M==N/15
    //M = ((int16_t)N * 37450) >> 18; // M==N/7
    N = N - (int16_t)M * 15;
    NASSERT(N<15);
    for (n = 0; n<M; n++)
    {
        BBE_LAVNX16_XP(s0, s_align, S, (BBE_SIMD_WIDTH - 1) * 2);
        s1 = BBE_SELNX16I(s0, st, BBE_SELI_PACK_1);
        // update status register
        st = BBE_REPNX16(s0, BBE_SIMD_WIDTH - 2);
        s1 = BBE_SRLINX16(s1, 1);
        s0 = BBE_SLLNX16(s0, sh15);
        s1 = BBE_ORNX16(s1, s0);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_4);
        //generate 256 bits of a convolution encoder
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_8);
        s1 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_12);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        s1 = BBE_INTLVNX16X1H(y0, y1);
        s0 = BBE_INTLVNX16X1L(y0, y1);
        BBE_SAVNX16_XP(s0, e_align, E, BBE_SIMD_WIDTH * 2);
        BBE_SAVNX16_XP(s1, e_align, E, (BBE_SIMD_WIDTH - 2) * 2);
    }
    if (N)  // last unprocessed data words
    {
        vselN sel0;
        BBE_LAVNX16_XP(s0, s_align, S, 2 * N);
        s1 = BBE_SELNX16I(s0, st, BBE_SELI_PACK_1);
        // update status register
        sel0 = BBE_MOVVSELNX16(N - 1, 0);
        st = BBE_SHFLNX16(s0, sel0);
        s1 = BBE_SRLINX16(s1, 1);
        s0 = BBE_SLLNX16(s0, sh15);
        s1 = BBE_ORNX16(s1, s0);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_4);
        //generate 256 bits of a convolution encoder
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_8);
        s1 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_12);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        s1 = BBE_INTLVNX16X1H(y0, y1);
        s0 = BBE_INTLVNX16X1L(y0, y1);
        N <<= 2;
        BBE_SAVNX16_XP(s0, e_align, E, N); N -= BBE_SIMD_WIDTH * 2;
        BBE_SAVNX16_XP(s1, e_align, E, N);
    }
    BBE_SAPOS_FP(e_align, E);
    return pEnc->state = (uint16_t)BBE_MOVAV16(st) >> (16 - K);
} /* convenc_process2() */
#endif

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
DISCARD_FUN(int16_t, convenc_process3, (convenc_handle_t handle, 
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

// interleave 64x3->192 bits: several variants
// workable variant
#define INTERLEAVE0(s,sel0,sel1,sh0,sh1) \
{                                  \
    xb_vecNx16 m0,m1,x0,x1;        \
    /* group bits 4x3->12: */      \
    /* fed_ba9_765_321_ -> */      \
    /* 0000fedba9765321    */      \
    m0=BBE_MOVVA16(0x0E0E);        \
    m1=BBE_MOVVA16(0xE0E0);        \
    x0=BBE_ANDNX16(m0,s);          \
    x1=BBE_ANDNX16(m1,s);          \
    x0=BBE_SRLINX16(x0,1);         \
    x1=BBE_SRLINX16(x1,2);         \
    s =BBE_ADDNX16(x0,x1);         \
    m0=BBE_MOVVA16(0x003F);        \
    m1=BBE_MOVVA16(0x3F00);        \
    x0=BBE_ANDNX16(m0,s );         \
    x1=BBE_ANDNX16(m1,s );         \
    x1=BBE_SRLINX16(x1,2);         \
    s =BBE_ADDNX16(x0,x1);         \
    /* group by 8x12->96 bits */   \
    x0=BBE_SELNX16(s ,s ,sel0.v);  \
    x1=BBE_SELNX16(s ,s ,sel1.v);  \
    x0=BBE_SRLNX16(x0,sh0.v);      \
    x1=BBE_SLLNX16(x1,sh1.v);      \
    s =BBE_ORNX16(x0,x1);          \
}

// more instructions but better parallelization
#define INTERLEAVE1(s,sel0,sel1,sh0,sh1) \
{                                  \
    xb_vecNx40 w;                  \
    xb_vecNx16 m0,m1,x0,x1,x2,x3;  \
    xb_vecNx16 c0,c1,c2,c3;        \
    /* group bits 4x3->12: */      \
    /* fed_ba9_765_321_ -> */      \
    /* 0000fedba9765321    */      \
    m0=BBE_MOVVA16(0x000E);        \
    m1=BBE_MOVVA16(0x00E0);        \
    x0=BBE_ANDNX16(s,m0);          \
    x1=BBE_ANDNX16(s,m1);          \
    m0=BBE_MOVVA16(0x0E00);        \
    m1=BBE_MOVVA16(0xE000);        \
    x2=BBE_ANDNX16(s,m0);          \
    x3=BBE_ANDNX16(s,m1);          \
    c0=BBE_MOVVA16(0x4000);        \
    c1=BBE_MOVVA16(0x2000);        \
    c2=BBE_MOVVA16(0x1000);        \
    c3=BBE_MOVVA16(0x0800);        \
    w=BBE_MULUUNX16(x0,c0);        \
    BBE_MULUUANX16(w,x1,c1);       \
    BBE_MULUUANX16(w,x2,c2);       \
    BBE_MULUUANX16(w,x3,c3);       \
    s=BBE_PACKQNX40(w);            \
    /* group by 8x12->96 bits */   \
    x0=BBE_SELNX16(s ,s ,sel0.v);  \
    x1=BBE_SELNX16(s ,s ,sel1.v);  \
    x0=BBE_SRLNX16(x0,sh0.v);      \
    x1=BBE_SLLNX16(x1,sh1.v);      \
    s =BBE_ORNX16(x0,x1);          \
}

#define INTERLEAVE2(s,sel0,sel1,sh0,sh1) \
{                                  \
    xb_vecNx40 w0 ;                \
    xb_vecNx16 m0,m1,x0,x1,x2,x3;  \
    xb_vecNx16 c0,c1,c2,c3;        \
    /* group bits 4x3->12: */      \
    /* fed_ba9_765_321_ -> */      \
    /* 0000fedba9765321    */      \
    m0=BBE_MOVVA16(0x000E);        \
    m1=BBE_MOVVA16(0x00E0);        \
    x0=BBE_ANDNX16(s,m0);          \
    x1=BBE_ANDNX16(s,m1);          \
    m0=BBE_MOVVA16(0x0E00);        \
    m1=BBE_MOVVA16(0xE000);        \
    x2=BBE_ANDNX16(s,m0);          \
    x3=BBE_ANDNX16(s,m1);          \
    c0=BBE_MOVVA16(0x4000);        \
    c1=BBE_MOVVA16(0x2000);        \
    c2=BBE_MOVVA16(0x1000);        \
    c3=BBE_MOVVA16(0x0800);        \
    w0=BBE_MULUUNX16(x0,c0);       \
    BBE_MULUUANX16(w0,x1,c1);      \
    BBE_MULUUANX16(w0,x2,c2);      \
    BBE_MULUUANX16(w0,x3,c3);      \
    s=BBE_PACKQNX40(w0);           \
    /* group by 8x12->96 bits */   \
    x0=BBE_SELNX16(s ,s ,sel0.v);  \
    x1=BBE_SELNX16(s ,s ,sel1.v);  \
    w0=BBE_UNPKUNX16(x0);          \
    x0=BBE_PACKVNX40(w0,sh0.v);    \
    x1=BBE_SLLNX16(x1,sh1.v);      \
    s =BBE_ORNX16(x0,x1);          \
}

#define INTERLEAVE3(s,sel0,sel1,sh0,sh1) \
{                                  \
    xb_vecNx40 w0 ;                \
    xb_vecNx16 m0,m1,x0,x1,x2,x3;  \
    xb_vecNx16 c0,c1,c2,c3,csh1;   \
    csh1=BBE_MOVVA16(1);           \
    csh1=BBE_SLLNX16(csh1,sh1 );   \
    /* group bits 4x3->12: */      \
    /* fed_ba9_765_321_ -> */      \
    /* 0000fedba9765321    */      \
    m0=BBE_MOVVA16(0x000E);        \
    m1=BBE_MOVVA16(0x00E0);        \
    x0=BBE_ANDNX16(s,m0);          \
    x1=BBE_ANDNX16(s,m1);          \
    m0=BBE_MOVVA16(0x0E00);        \
    m1=BBE_MOVVA16(0xE000);        \
    x2=BBE_ANDNX16(s,m0);          \
    x3=BBE_ANDNX16(s,m1);          \
    c0=BBE_MOVVA16(0x4000);        \
    c1=BBE_MOVVA16(0x2000);        \
    c2=BBE_MOVVA16(0x1000);        \
    c3=BBE_MOVVA16(0x0800);        \
    w0=BBE_MULUUNX16(x0,c0);       \
    BBE_MULUUANX16(w0,x1,c1);      \
    BBE_MULUUANX16(w0,x2,c2);      \
    BBE_MULUUANX16(w0,x3,c3);      \
    s=BBE_PACKQNX40(w0);           \
    /* group by 8x12->96 bits */   \
    x0=BBE_SELNX16(s ,s ,sel0);    \
    x1=BBE_SELNX16(s ,s ,sel1);    \
    w0=BBE_UNPKUNX16(x0);          \
    w0=BBE_SRLNX40(w0,sh0  );      \
    BBE_MULUUANX16(w0,x1,csh1);    \
    s =BBE_PACKLNX40(w0);          \
}

int16_t convenc_process3 ( convenc_handle_t handle, 
                           int16_t * restrict e, const int16_t * s, int N )
{
    static const union { int16_t i[16]; _vselN s; } ALIGN(32) sel0 = { { 0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14 } };
    static const union { int16_t i[16]; _vsaN  v; } ALIGN(32) sh0 = { { 0, 4, 8, 0, 4, 8, 0, 4, 8, 0, 4, 8, 16, 16, 16, 16 } };
    static const union { int16_t i[16]; _vselN s; } ALIGN(32) sel1 = { { 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15 } };
    static const union { int16_t i[16]; _vsaN  v; } ALIGN(32) sh1 = { { 12, 8, 4, 12, 8, 4, 12, 8, 4, 12, 8, 4, 16, 16, 16, 16 } };
    tConvenc_ *pEnc = (tConvenc_ *)handle;
    const xb_vecNx16   * restrict S = (xb_vecNx16 *)s;
    xb_vecNx16   * restrict E = (xb_vecNx16 *)e;
    int n, K, M;
    uint16_t p0, p1, p2;
    valign s_align, e_align;
    xb_vecNx16    s0, s1, s3, s4, st, y0, y1, y2;
    xb_vecNx16    z1;
    vselN sel2;

    NASSERT(pEnc != NULL);
    NASSERT(pEnc->R == 3);
    NASSERT(pEnc->magic == MAGIC);
    p0 = pEnc->polyMask[0] * 2;
    p1 = pEnc->polyMask[1] * 2;
    p2 = pEnc->polyMask[2] * 2;
    K = pEnc->K;
    st = BBE_MOVVA16(pEnc->state << (16 - K));
    s_align = BBE_LAVNX16_PP(S);
    e_align = BBE_ZALIGN();
    z1 = BBE_ZERONX16();
    y0 = y1 = y2 = BBE_ZERONX16();

    M = ((int16_t)N * 43691) >> 19; // M==N/12
    N = N - (int16_t)(M * 12);
    NASSERT(N<12);

    for (n = 0; n<M; n++)
    {
        BBE_LAVNX16_XP(s0, s_align, S, 2 * 12);
        s1 = BBE_SELNX16I(s0, st, BBE_SELI_PACK_1);
        s1 = BBE_SRLINX16(s1, 1);
        // update status register
        st = BBE_REPNX16(s0, 11);
        s0 = BBE_SLLINX16(s0, 15);
        s1 = BBE_ORNX16(s1, s0);
        s0 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_4);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_8);
        s4 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_12);
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        BBE_CC64(y2, s1, p2);
        BBE_CC64(y0, s0, p0);
        BBE_CC64(y1, s0, p1);
        BBE_CC64(y2, s0, p2);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        BBE_CC64(y2, s3, p2);
        BBE_CC64(y0, s4, p0);
        BBE_CC64(y1, s4, p1);
        BBE_CC64(y2, s4, p2);

        s0 = BBE_INTLVNX16X1H(y0, y2);
        s1 = BBE_INTLVNX16X1H(y1, z1);
        s0 = BBE_INTLVNX16X1L(s0, s1);
        y0 = BBE_INTLVNX16X1L(y0, y2);
        y1 = BBE_INTLVNX16X1L(y1, z1);
        s3 = BBE_INTLVNX16X1H(y0, y1);
        s4 = BBE_INTLVNX16X1L(y0, y1);
        INTERLEAVE3(s4, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));
        INTERLEAVE3(s3, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));
        INTERLEAVE3(s0, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));

        BBE_SAVNX16_XP(s4, e_align, E, 2 * 12);
        BBE_SAVNX16_XP(s3, e_align, E, 2 * 12);
        BBE_SAVNX16_XP(s0, e_align, E, 2 * 12);
    }
    //generate last bits
    if (N)
    {
        BBE_LAVNX16_XP(s0, s_align, S, 2 * N);
        s1 = BBE_SELNX16I(s0, st, BBE_SELI_PACK_1);
        s1 = BBE_SRLINX16(s1, 1);
        // update status register
        sel2 = BBE_MOVVSELNX16(N - 1, 0);
        st = BBE_SHFLNX16(s0, sel2);
        s0 = BBE_SLLINX16(s0, 15);
        s1 = BBE_ORNX16(s1, s0);
        s0 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_4);
        s3 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_8);
        s4 = BBE_SELNX16I(s1, s1, BBE_SELI_ROTATE_RIGHT_12);
        BBE_CC64(y0, s1, p0);
        BBE_CC64(y1, s1, p1);
        BBE_CC64(y2, s1, p2);
        BBE_CC64(y0, s0, p0);
        BBE_CC64(y1, s0, p1);
        BBE_CC64(y2, s0, p2);
        BBE_CC64(y0, s3, p0);
        BBE_CC64(y1, s3, p1);
        BBE_CC64(y2, s3, p2);
        BBE_CC64(y0, s4, p0);
        BBE_CC64(y1, s4, p1);
        BBE_CC64(y2, s4, p2);

        s0 = BBE_INTLVNX16X1H(y0, y2);
        s1 = BBE_INTLVNX16X1H(y1, z1);
        s0 = BBE_INTLVNX16X1L(s0, s1);
        y0 = BBE_INTLVNX16X1L(y0, y2);
        y1 = BBE_INTLVNX16X1L(y1, z1);
        s3 = BBE_INTLVNX16X1H(y0, y1);
        s4 = BBE_INTLVNX16X1L(y0, y1);

        INTERLEAVE3(s4, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));
        INTERLEAVE3(s3, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));
        INTERLEAVE3(s0, _S(sel0.s), _S(sel1.s), _V(sh0.v), _V(sh1.v));

        N = 6 * N;
        M = XT_MIN(24, N);
        BBE_SAVNX16_XP(s4, e_align, E, M); N -= M;
        M = XT_MIN(24, N);
        BBE_SAVNX16_XP(s3, e_align, E, M); N -= M;
        M = XT_MIN(24, N);
        BBE_SAVNX16_XP(s0, e_align, E, M);
    }
    BBE_SAPOS_FP(e_align, E);
    pEnc->state = (uint16_t)BBE_MOVAV16(st) >> (16 - K);
    return pEnc->state;
} /* convenc_process3() */
#endif

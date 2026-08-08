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
    CRC
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

#if !(HAVE_LFSR && 1)
DISCARD_FUN(size_t, crc_alloc,(int order, uint32_t poly))
DISCARD_FUN(crc_handle_t, crc_init,(void * objmem, 
                      int order, uint32_t poly))
DISCARD_FUN(uint32_t, crc_process,(crc_handle_t handle, uint32_t reg, const uint8_t* bitstream, int N))
#else
/*-------------------------------------------------------------------------
CRC 

Description: This function calculates 8-, 16-, 24- or 32-bit CRC checksum

Parameters:
Input:
objmem          Allocated memory block
bitstream[N]    Bitstream
order           Order of polynomial (8,16,24 or 32)
poly            Polynomial without higher degree, for example: 
                x^16 + x^12 + x^5 + 1 corresponds to 0x1021
reg             Original crc register value

Returned value  Updated crc register 

Restrictions:
bistream must be aligned on 2-byte boundary
-------------------------------------------------------------------------*/

#define MAGIC 0xe94b0fb2

#define BBE_MOVVW(hi,lo,w) { hi=BBE_MOVVWH(w);lo=BBE_MOVVWL(w);}

typedef struct
{
    uint32_t magic;
    int order;
    uint32_t poly;
    uint32_t * multbl;
}
tCrc_;

size_t crc_alloc ( int order, uint32_t poly )
{
    return sizeof(tCrc_) + 7 + 31 + sizeof(uint32_t) * 32;
} /* crc_alloc() */

crc_handle_t crc_init ( void * objmem, int order, uint32_t poly )
{
    xb_int16U  * restrict pTbl;
    xb_vecNx40 wt, wrez0, wix0, wpoly;
    xb_vecNx40 wrez1, wix1;
    xb_vecNx16 t, t1, t0;
    vsaN sh;
    static const union { int16_t i[16]; _vsaN v; }
    ALIGN(32) shw[] = { { { 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23 } },
    { { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 } }
    };
    tCrc_ *pCrc;
    uintptr_t a;
    int k;
    NASSERT((order == 8) || (order == 16) || (order == 24) || (order == 32));
    if (objmem == NULL) return NULL;
    a = (uintptr_t)objmem;
    a = (a + 7)&(~7);
    pCrc = (tCrc_*)a;
    a = (a + sizeof(tCrc_) + 31)&(~31);
    pCrc->multbl = (uint32_t *)a;
    pCrc->magic = MAGIC;
    pCrc->order = order;
    poly <<= (32 - order);
    pCrc->poly = poly;
    pTbl = (xb_int16U   *)pCrc->multbl;

    /*
    multbl is a table with remainders of division 2^k, k=0...31 to the given polynomial
    this table is written in the transposed form to be directly used by BBE_BMUL32
    instruction
    */
    wpoly = BBE_MOVWAU32((poly));
    wix0 = BBE_MOVWA32(1);
    wix0 = BBE_SLLNX40(wix0, _V(shw[0].v));
    wrez0 = BBE_ZERONX40();

    for (k = 0; k<32; k++)
    {
        wt = BBE_XORNX40(wrez0, wix0);
        wrez0 = BBE_SLLINX40(wrez0, 1);
        wix0 = BBE_SLLINX40(wix0, 1);
        wt = BBE_SLLINX40(wt, 8);
        wt = BBE_SRAINX40(wt, 39);
        wt = BBE_ANDNX40(wt, wpoly);
        wrez0 = BBE_XORNX40(wrez0, wt);
    }

    wix1 = BBE_MOVWA32(1);
    wix1 = BBE_SLLNX40(wix1, _V(shw[1].v));
    wrez1 = BBE_ZERONX40();

    for (k = 0; k<32; k++)
    {
        wt = BBE_XORNX40(wrez1, wix1);
        wrez1 = BBE_SLLINX40(wrez1, 1);
        wix1 = BBE_SLLINX40(wix1, 1);
        wt = BBE_SLLINX40(wt, 8);
        wt = BBE_SRAINX40(wt, 39);
        wt = BBE_ANDNX40(wt, wpoly);
        wrez1 = BBE_XORNX40(wrez1, wt);
    }

    // byte reverse 
    BBE_MOVVW(t1, t0, wrez0);
    t = BBE_SLLINX16(t0, 8);
    t0 = BBE_SRLINX16(t0, 8);
    t0 = BBE_ORNX16(t0, t);
    t = BBE_SLLINX16(t1, 8);
    t1 = BBE_SRLINX16(t1, 8);
    t1 = BBE_ORNX16(t1, t);
    t0 = BBE_SHFLNX16I(t0, BBE_SHFLI_SWAP_1);   //1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14
    t1 = BBE_SHFLNX16I(t1, BBE_SHFLI_SWAP_1);
    wrez0 = BBE_MOVWV(t1, t0);

    BBE_MOVVW(t1, t0, wrez1);
    t = BBE_SLLINX16(t0, 8);
    t0 = BBE_SRLINX16(t0, 8);
    t0 = BBE_ORNX16(t0, t);
    t = BBE_SLLINX16(t1, 8);
    t1 = BBE_SRLINX16(t1, 8);
    t1 = BBE_ORNX16(t1, t);
    t0 = BBE_SHFLNX16I(t0, BBE_SHFLI_SWAP_1);   //1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14
    t1 = BBE_SHFLNX16I(t1, BBE_SHFLI_SWAP_1);
    wrez1 = BBE_MOVWV(t1, t0);

    // transpose ans save resulted table
    t = BBE_SEQNX16();
    t0 = BBE_MOVVINT16(39);
    t = BBE_SUBNX16(t0, t);
    sh = BBE_MOVVSV(t, 0);
    for (k = 0; k<32; k++)
    {
        xb_int40 i0;
        xb_int16U x;

        wt = BBE_SLLINX40(wrez0, 39);
        wt = BBE_SRLNX40(wt, sh);
        i0 = BBE_RADDNX40(wt);
        x = xb_int40_rtor_xb_int16U(i0);
        xb_int16U_storeip(x, pTbl, 2);

        wt = BBE_SLLINX40(wrez1, 39);
        wt = BBE_SRLNX40(wt, sh);
        i0 = BBE_RADDNX40(wt);
        x = xb_int40_rtor_xb_int16U(i0);
        xb_int16U_storeip(x, pTbl, 2);

        wrez0 = BBE_SRLINX40(wrez0, 1);
        wrez1 = BBE_SRLINX40(wrez1, 1);
    }

    return (crc_handle_t)pCrc;
} /* crc_init() */

uint32_t crc_process ( crc_handle_t handle, uint32_t reg, const uint8_t * bitstream, int N )
{
    tCrc_ *pCrc = (tCrc_ *)handle;
    uint32_t a, b;
    int n;
    vsaN shN;
    xb_vecNx40 w;

    xb_vecNx16    p00, p01;
    xb_vecNx16    c, y, t, z;
    const xb_vecNx16* restrict B;
    valign b_align;
    const xb_vecNx16   * restrict P = (const xb_vecNx16 *)pCrc->multbl;
    z = BBE_ZERONX16();

    // load polynomial to BMUL32 unit
    p00 = BBE_LVNX16_I(P, 0 * 2 * BBE_SIMD_WIDTH);
    p01 = BBE_LVNX16_I(P, 1 * 2 * BBE_SIMD_WIDTH);
    BBE_MOVBMULSTATEV(p01, p00, 0);
    p00 = BBE_LVNX16_I(P, 2 * 2 * BBE_SIMD_WIDTH);
    p01 = BBE_LVNX16_I(P, 3 * 2 * BBE_SIMD_WIDTH);
    BBE_MOVBMULSTATEV(p01, p00, 1);

    B = (const xb_vecNx16*)bitstream;
    b_align = BBE_LA_PP(B);
    // load status register, shift left by (32-pCrc->order), reverse bytes and put it into 14,15 positions 
    w = BBE_MOVWAU32(reg);
    shN = BBE_MOVVSA32(32 - pCrc->order);
    w = BBE_SLLNX40(w, shN);
    BBE_MOVVW(t, y, w);
    t = BBE_SLLINX16(y, 8);
    y = BBE_SRLINX16(y, 8);
    y = BBE_ORNX16(y, t);
    y = BBE_SHFLNX16I(y, BBE_SHFLI_REVERSE_1); //15,14,13,12...
    // main loop: unrolled 4 times: 16-byte per iteration
    BBE_LAVNX16_XP(c, b_align, B, 4);
    b = BBE_EXTRNX16C(c, 0);
    b ^= BBE_EXTRNX16C(y, 0);
    BBE_MOVBMULACCA(b);
    for (n = 0; n<((N - 4)&(~15)); n += 16)
    {
        uint32_t b;
        BBE_LAVNX16_XP(c, b_align, B, 16);
        b = BBE_EXTRNX16C(c, 0); BBE_BMUL32A(y, b);
        b = BBE_EXTRNX16C(c, 1); BBE_BMUL32A(y, b);
        b = BBE_EXTRNX16C(c, 2); BBE_BMUL32A(y, b);
        b = BBE_EXTRNX16C(c, 3); BBE_BMUL32A(y, b);
    }
    for (; n<(N&(~3)); n += 4)
    {
        uint32_t b;
        BBE_LAVNX16_XP(c, b_align, B, 4);
        b = BBE_EXTRNX16C(c, 0);
        BBE_BMUL32A(y, b);
    }
    // process tail
    a = BBE_MOVABMULACC();
    b = BBE_EXTRNX16C(c, 0);
    b ^= a;
    // rotate y by N bits right and c M bits right
    XT_SSA8L(N & 3);
    a = XT_SLL(a);
    b = XT_SRL(b);
    // XT_MOVEQZ(a,0,N);
    BBE_MOVBMULACCA(a);
    y = BBE_BMUL32A(y, b);

    // take from last 2 words, byte reverse and shift right by (32-pCrc->order)
    y = BBE_SHFLNX16I(y, BBE_SHFLI_REVERSE_1); //15,14,13,12...
    t = BBE_SLLINX16(y, 8);
    y = BBE_SRLINX16(y, 8);
    y = BBE_ORNX16(y, t);
    y = BBE_SELNX16I(z, y, BBE_SELI_PACK_2); //0, 1,16,17,18,19,20...
    reg = BBE_EXTRNX16C(y, 0);
    reg >>= (32 - pCrc->order);
    return reg;
} /* crc_process() */
#endif

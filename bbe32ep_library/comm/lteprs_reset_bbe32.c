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
    Reset LTE PRS generator
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
Reset LTE PRS generator

Function lteprs_reset() prepares shift registers and fills them with proper
values derived from initial values (1 for first register and cinit for second
one).

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Output:
  r[2]   32-bit LFSR states

Input:
  cinit  Initialization value for second shift register
  M      Codeword size (1, 2, 4, 6, 8 10 and 16 bits)

Return value:
  none
-------------------------------------------------------------------------*/

void lteprs_reset ( uint32_t * r, uint32_t cinit, int16_t M )
{
#if !HAVE_LFSR
    int i;
    uint32_t rx2, rx1;

    rx2 = cinit & 0x7fffffff;

    if (M == 16)
    {
        rx1 = 1;
        for (i = 0; i<58; ++i)
        {
            rx1 = ((rx1 >> 27) | ((((rx1 >> 3) ^ rx1)) << 4));
        }
        /* rx1 contains bits 1566..1597. Calculate last 18 bits: */
        rx1 = ((rx1 >> 18) | ((((rx1 >> 3) ^ rx1)) << 13));
        r[0] = rx1;

        for (i = 0; i<58; ++i)
        {
            rx2 = ((rx2 >> 27) | ((((rx2 >> 3) ^ (rx2 >> 2) ^ (rx2 >> 1) ^ rx2)) << 4));
        }
        /* rx1 contains bits 1566..1597. Calculate last 18 bits: */
        rx2 = ((rx2 >> 18) | ((((rx2 >> 3) ^ (rx2 >> 2) ^ (rx2 >> 1) ^ rx2)) << 13));
    }
    else
    {
        for (i = 0; i<59; ++i)
        {
            rx2 = ((rx2 >> 27) | ((((rx2 >> 3) ^ (rx2 >> 2) ^ (rx2 >> 1) ^ rx2)) << 4)); //NOTE: use 31 bit as 0 bit of next state (to simplify BBE calculations)
        }
        /* rx2 contains bits 1593..1624. Calculate last 7 bits: */
        rx2 = ((rx2 >> 7) | ((((rx2 >> 3) ^ (rx2 >> 2) ^ (rx2 >> 1) ^ rx2)) << 24)); //NOTE: use 31 bit as 0 bit of next state (to simplify BBE calculations)

        r[0] = 0x5e485840; /* precalculated */
    }
    r[1] = rx2;
#else
    int i;
    uint32_t cinitx2, StateLo, StateHi, State;
    xb_vecNx16 r1;
    xb_vecNx16 p00, p01, p02, p03, p10, p11, p12, p13;

    static const uint32_t ALIGN(32) poly[][8] = {
        { 0x0000001E, 0x0000003C, 0x00000078, 0x000000F0, 0x000001E0, 0x000003C0, 0x00000780, 0x00000F00 },
        { 0x00001E00, 0x00003C00, 0x00007800, 0x0000F000, 0x0001E000, 0x0003C000, 0x00078000, 0x000F0000 },
        { 0x001E0000, 0x003C0000, 0x00780000, 0x00F00000, 0x01E00000, 0x03C00000, 0x07800000, 0x0F000000 },
        { 0x1E000000, 0x3C000000, 0x78000000, 0xF0000000, 0xe000001e, 0xc0000022, 0x8000005a, 0x000000aa },
        { 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100 },
        { 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000 },
        { 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000 },
        { 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000, 0x0000001E },
    };
    const xb_vecNx16 * restrict P = (const xb_vecNx16 *)poly;
    cinitx2 = cinit * 2;


    p00 = BBE_LVNX16_I(P, 0 * BBE_SIMD_WIDTH);
    p01 = BBE_LVNX16_I(P, 2 * BBE_SIMD_WIDTH);
    p02 = BBE_LVNX16_I(P, 4 * BBE_SIMD_WIDTH);
    p03 = BBE_LVNX16_I(P, 6 * BBE_SIMD_WIDTH);
    p10 = BBE_LVNX16_I(P, 8 * BBE_SIMD_WIDTH);
    p11 = BBE_LVNX16_I(P, 10 * BBE_SIMD_WIDTH);
    p12 = BBE_LVNX16_I(P, 12 * BBE_SIMD_WIDTH);
    p13 = BBE_LVNX16_I(P, 14 * BBE_SIMD_WIDTH);

    BBE_MOVBMULSTATEV(p01, p00, 0);
    BBE_MOVBMULSTATEV(p03, p02, 1);
    BBE_MOVBMULACCA(cinitx2);
    if (M == 16)
    {
        r1 = BBE_LPNX16_I(&cinitx2, 0);
        r1 = BBE_REPNX16C(r1, 0);

        //calculate 1599 bits
        for (i = 0; i<50 - 1; ++i)
        {
            BBE_BMUL32A(r1, 0);
        }
        StateLo = BBE_MOVABMULACC();
        // need to calculate 16 bits (summary 1568+16=1584 bits):
        BBE_BMUL32A(r1, 0); // calculate 32 bits
        StateHi = BBE_MOVABMULACC();
        State = (StateHi << 16) | (StateLo >> 16);
        BBE_MOVBMULACCA(State);

        //OLD:r[0] = 0xA9A43648UL; /* 0xA9A43648 precalculated */
        r[0] = 0x5840A9A4UL; /*0x5840A9A4 precalculated */
    }
    else
    {
        r1 = 0;

        //calculate 1599 bits
        for (i = 0; i<50; ++i)
        {
            BBE_BMUL32A(r1, 0);
        }
        
        r[0] = 0x5e485840UL; /* 0x5e485840 precalculated */
    }

    //calculate last bit
    BBE_MOVBMULSTATEV(p11, p10, 0);
    BBE_MOVBMULSTATEV(p13, p12, 1);
    BBE_BMUL32A(r1, 0);
    r[1] = BBE_MOVABMULACC();
#endif
} /* lteprs_reset() */

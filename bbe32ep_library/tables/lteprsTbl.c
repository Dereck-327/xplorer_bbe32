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
	NatureDSP_Baseband library. Communication part
    LTE PRS generator
	Integrit, 2006-2017
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "lteprs_common.h"

// NOTE: Bits 1...31 contain the polynomial state for the NEXT output bit

// NOTE: In standard algorithm bit 31 is unused.
//       In BBE-optimized algorithm bit 31 used to store bit 0 of next state
const uint32_t ALIGN(32) lteprs_poly[8][8] = {
    //poly0 = (1<<3 | 1);  
    { 0x00000012, 0x00000024, 0x00000048, 0x00000090, 0x00000120, 0x00000240, 0x00000480, 0x00000900 },
    { 0x00001200, 0x00002400, 0x00004800, 0x00009000, 0x00012000, 0x00024000, 0x00048000, 0x00090000 },
    { 0x00120000, 0x00240000, 0x00480000, 0x00900000, 0x01200000, 0x02400000, 0x04800000, 0x09000000 },
    { 0x12000000, 0x24000000, 0x48000000, 0x90000000, 0x20000012, 0x40000024, 0x80000048, 0x00000082 },
    //poly1 = (1<<3 | 1<<2 | 1<<1 | 1); 
    { 0x0000001E, 0x0000003C, 0x00000078, 0x000000F0, 0x000001E0, 0x000003C0, 0x00000780, 0x00000F00 },
    { 0x00001E00, 0x00003C00, 0x00007800, 0x0000F000, 0x0001E000, 0x0003C000, 0x00078000, 0x000F0000 },
    { 0x001E0000, 0x003C0000, 0x00780000, 0x00F00000, 0x01E00000, 0x03C00000, 0x07800000, 0x0F000000 },
    { 0x1E000000, 0x3C000000, 0x78000000, 0xF0000000, 0xE000001E, 0xC0000022, 0x8000005A, 0x000000AA }
};

const tlteprssel ALIGN(32) lteprs2_tbl[2] = {
    { { 0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14 } },
    { { 0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15 } }
};

const tlteprssel ALIGN(32) lteprs4_tbl[2] = {
    { { 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13 } },
    { { 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12 } }
};

const tlteprssel ALIGN(32) lteprs6_tbl[4] =
{
    { { 8, 8, 8, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13 } },
    { { 0, 0, 9, 0, 0, 10, 0, 0, 0, 0, 12, 0, 0, 13, 0, 0 } },
    { { 0, 6, 12, 2, 8, 14, 4, 10, 0, 6, 12, 2, 8, 14, 4, 10 } },
    { { 0, 0, 16, 0, 0, 4, 0, 0, 0, 0, 16, 0, 0, 4, 0, 0 } }
};

const tlteprssel ALIGN(32) lteprs8_tbl[1] =
{
    { { 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29 } }
};

const tlteprssel ALIGN(32) lteprs10_tbl[4] =
{
    { { 4, 4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13 } }, // sel0
    { { 0, 5, 0, 6, 7, 0, 8, 0, 0, 10, 0, 11, 12, 0, 13, 0 } }, // sel1
    { { 0, 10, 4, 14, 8, 2, 12, 6, 0, 10, 4, 14, 8, 2, 12, 6 } }, // shft
    { { 0, 64, 0, 4, 256, 0, 16, 0, 0, 64, 0, 4, 256, 0, 16, 0 } }  // prod
};

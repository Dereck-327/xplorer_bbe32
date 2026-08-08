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
    Universal Grey coding/QAM mapping
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
  Universal Grey coding/QAM mapping

  Map data n-tuples (n=1..8) onto a user-defined QAM<2^n> constellation.

  User provides a constellation by passing an array S of 16-bit complex numbers,
  where the correspondence between a QAM point s and its related data item m
  is established through the array index: s_I = S[2*m+0] (in-phase component), 
  s_Q = S[2*m+1] (quadrature phase component).

  Input:
  S[2*(2^n)]  Complex constellation
  d[N]        Input n-tuples

  Output:
  r[2*N]      Coded complex points

  Scratch:
  pScr        Scratch memory for storing temporary data. Allocated space must 
              be at least GCODE_SCRATCH() bytes

  Restrictions:
  N           Multiple of 16. 
  pScr,r,d    Aligned on 32-byte boundary.
---------------------------------------------------------------------------*/

void gcode ( void * pScr,
             int16_t * restrict r,
       const int16_t *          d,
       const int16_t *          S,
       int N )
{
    int n;
    int id0, id1;
    const int16_t    * restrict pS = (int16_t *)S;
    const xb_vecNx16 * restrict vD = (const xb_vecNx16*)d;
    xb_vecNx16 * restrict R = (xb_vecNx16 *)r;
    xb_vecNx16  d0;
    xb_vecNx16  s0, y0, s1, y1;
    xb_vecNx16  z0;
    xb_vecNx16  id_vec0, id_vec1;

    if (N <= 0) return;
    ASSERT(r && d && S && N >= 0);
    NASSERT_ALIGN32(r);
    NASSERT_ALIGN32(d);
    NASSERT((N & 0xf) == 0);

    z0 = BBE_MOVVA16(0x4);

#ifdef COMPILER_XTENSA
#pragma ymemory(R)
#pragma ymemory(vD)
#pragma ymemory(pS)
#endif
    for (n = 0; n<N / (BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(d0, vD, (2 * BBE_SIMD_WIDTH));
        s0 = BBE_SELNX16I(BBE_ZERONX16(), d0, BBE_SELI_INTERLEAVE_1_LO);
        s1 = BBE_SELNX16I(BBE_ZERONX16(), d0, BBE_SELI_INTERLEAVE_1_HI);
        s0 = BBE_MULNX16PACKL(s0, z0);
        s1 = BBE_MULNX16PACKL(s1, z0);
        id_vec0 = s0;
        id_vec1 = s1;
        //7 Point
        id0 = BBE_EXTRNX16C(id_vec0, 7);
        y0 = BBE_LPNX16_X(pS, id0);
        id1 = BBE_EXTRNX16C(id_vec1, 7);
        y1 = BBE_LPNX16_X(pS, id1);
        //6 Point
        id0 = BBE_EXTRNX16C(id_vec0, 6);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 6);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //5 Point
        id0 = BBE_EXTRNX16C(id_vec0, 5);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 5);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //4 Point
        id0 = BBE_EXTRNX16C(id_vec0, 4);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 4);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //3 Point
        id0 = BBE_EXTRNX16C(id_vec0, 3);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 3);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //2 Point
        id0 = BBE_EXTRNX16C(id_vec0, 2);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 2);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //1 Point
        id0 = BBE_EXTRNX16C(id_vec0, 1);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 1);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);
        //0 Point
        id0 = BBE_EXTRNX16C(id_vec0, 0);
        s0 = BBE_LPNX16_X(pS, id0);
        y0 = BBE_SELNX16I(y0, s0, BBE_SELI_PACK_2);
        id1 = BBE_EXTRNX16C(id_vec1, 0);
        s1 = BBE_LPNX16_X(pS, id1);
        y1 = BBE_SELNX16I(y1, s1, BBE_SELI_PACK_2);

        BBE_SVNX16_IP(y0, R, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(y1, R, 2 * BBE_SIMD_WIDTH);
    }
} /* gcode() */

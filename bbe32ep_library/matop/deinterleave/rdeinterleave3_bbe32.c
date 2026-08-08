/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Matrix Operations
    1-to-M complex/real streams deinterleave
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"
/*-------------------------------------------------------------------------
1-to-M complex/real streams deinterleave

Description: decompose the input data stream into M=2,3,4 or 8 output 
streams, element by element. Use rdeinterleave<M>() functions for real
data, and cdeinterleave<M>() functions - for complex data.

Representation:
<r|c>deinterleave<M>   16-bit fixed-point data
<r|c>deinterleave<M>f  IEEE-754 Std single precision floating-point data

Parameters:
Input:
M                    Number of streams
N                    Number of elements per each output stream
x[M*N]               Input data stream
Output:
y0[N],...,y<M-1>[N]  Deinterleaved data streams
y[M]                 M pointers to output data streams

Restrictions:
x, y0,...,y<M-1>,
y[0..M-1]            Must not overlap and must be aligned on 32-byte boundary
N                    Must be a multiple of:
                       16 for real fixed-point data 
                        8 for complex fixed-point data and real floating-point data
                        4 for complex floating-point data
-------------------------------------------------------------------------*/

/* M=3 */
void rdeinterleave3 ( int16_t * restrict y0,
                      int16_t * restrict y1,
                      int16_t * restrict y2, 
                const int16_t * restrict x, int N )
{
  xb_vecNx16* restrict w0 = (xb_vecNx16*)y0;
  xb_vecNx16* restrict w1 = (xb_vecNx16*)y1;
  xb_vecNx16* restrict w2 = (xb_vecNx16*)y2;
  const xb_vecNx16* restrict z = (const xb_vecNx16*)x;
  xb_vecNx16 X0, X1, X2, Y0, Y1, Y2;
  int n;
  vselN sel;  //BBE_SELI_INTERLEAVE_1_EVENODD but with exchanged neighbour elements

  X0 = BBE_SEQNX16(); X1 = BBE_MOVVA16(16); X1 = BBE_ADDNX16(X0, X1);
  X0 = BBE_SELNX16I(X1, X0, BBE_SELI_INTERLEAVE_1_EVENODD);
  X0 = BBE_SHFLNX16I(X0, BBE_SHFLI_SWAP_1);
  sel = BBE_MOVVSELNX16(X0, 0);

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y0, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y1, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y2, 2 * BBE_SIMD_WIDTH);
  NASSERT((N&(BBE_SIMD_WIDTH- 1)) == 0);
  for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    xb_vecNx16 Z0, Z1;
    BBE_LVNX16_IP(X0, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, z, 2 * BBE_SIMD_WIDTH);
    /* make complex deinterleaving first */
    BBE_DSELNX16I(Y2, Z0, X1, X0, BBE_DSELI_DEINTERLEAVE_C3_STEP_0);
    BBE_DSELNX16I_H(Y2, Z1, X1, X2, BBE_DSELI_DEINTERLEAVE_C3_STEP_1);
    BBE_DSELNX16I(Y1, Y0, Z1, Z0, BBE_DSELI_DEINTERLEAVE_2);
    /* convert complex deinterleaved pairs to real deinterleved data */
    X0 = BBE_SELNX16I(Y1, Y0, BBE_SELI_INTERLEAVE_1_EVENODD);/* 0,17,2,19,4,21,6,23... */
    X1 = BBE_SELNX16(Y0, Y2, sel);/*1,32,3,34,5,36,7,38,9,40...*/
    X2 = BBE_SELNX16I(Y2, Y1, BBE_SELI_INTERLEAVE_1_EVENODD);/* 16,33,18,35,20,37,22,39... */

    BBE_SVNX16_IP(X0, w0, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1, w1, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2, w2, 2 * BBE_SIMD_WIDTH);
  }
} /* rdeinterleave3() */

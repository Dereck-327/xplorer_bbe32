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

/* M=8 */
void rdeinterleave8 ( int16_t * restrict y0,
                      int16_t * restrict y1,
                      int16_t * restrict y2,
                      int16_t * restrict y3,
                      int16_t * restrict y4,
                      int16_t * restrict y5,
                      int16_t * restrict y6,
                      int16_t * restrict y7, 
                const int16_t * restrict x, int N )
{
  xb_vecNx16* restrict w0 = (xb_vecNx16*)y0;
  xb_vecNx16* restrict w1 = (xb_vecNx16*)y1;
  xb_vecNx16* restrict w2 = (xb_vecNx16*)y2;
  xb_vecNx16* restrict w3 = (xb_vecNx16*)y3;
  xb_vecNx16* restrict w4 = (xb_vecNx16*)y4;
  xb_vecNx16* restrict w5 = (xb_vecNx16*)y5;
  xb_vecNx16* restrict w6 = (xb_vecNx16*)y6;
  xb_vecNx16* restrict w7 = (xb_vecNx16*)y7;
  const xb_vecNx16* restrict z = (const xb_vecNx16*)x;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  int n;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y0, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y1, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y2, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y3, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y4, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y5, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y6, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y7, 2 * BBE_SIMD_WIDTH);
  NASSERT((N&(BBE_SIMD_WIDTH- 1)) == 0);

  for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    BBE_LVNX16_IP(Y0, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y2, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y3, z, 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(X1, X0, Y1, Y0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X3, X2, Y3, Y2, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y2, Y0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y3, Y1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
    BBE_LVNX16_IP(Y4, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y5, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y6, z, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y7, z, 2 * BBE_SIMD_WIDTH);

    BBE_DSELNX16I(X5, X4, Y5, Y4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7, X6, Y7, Y6, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y6, Y4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y7, Y5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(X1, X0, Y4, Y0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X3, X2, Y5, Y1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X5, X4, Y6, Y2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(X7, X6, Y7, Y3, BBE_DSELI_DEINTERLEAVE_1);

    BBE_SVNX16_IP(X0, w0, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1, w1, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2, w2, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X3, w3, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X4, w4, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X5, w5, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X6, w6, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X7, w7, 2 * BBE_SIMD_WIDTH);
  }
} /* rdeinterleave8() */

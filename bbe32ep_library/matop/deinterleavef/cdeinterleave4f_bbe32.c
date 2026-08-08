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

/* M=4 */
void cdeinterleave4f ( complex_float * restrict y0,
                       complex_float * restrict y1,
                       complex_float * restrict y2,
                       complex_float * restrict y3, 
                 const complex_float * restrict x, int N )
{
  const xb_vecNx16 * restrict px;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
        xb_vecNx16 * restrict py2;
        xb_vecNx16 * restrict py3;
  xb_vecNx16 V0, V1, V2, V3;
  int n;

  NASSERT_ALIGN(x , 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y0, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y1, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y2, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y3, 2*BBE_SIMD_WIDTH);
  NASSERT((N & (BBE_SIMD_WIDTH/4-1)) == 0);

  px  = (const xb_vecNx16 *)(x);
  py0 = (      xb_vecNx16 *)(y0);
  py1 = (      xb_vecNx16 *)(y1);
  py2 = (      xb_vecNx16 *)(y2);
  py3 = (      xb_vecNx16 *)(y3);

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++)
  {
    BBE_LVNX16_IP(V0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(V1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(V2, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(V3, px, 2*BBE_SIMD_WIDTH);

    BBE_DSELNX16I(V2, V0, V2, V0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(V3, V1, V3, V1, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(V1, V0, V1, V0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(V3, V2, V3, V2, BBE_DSELI_INTERLEAVE_4);

    BBE_SVNX16_IP(V0, py0, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(V1, py1, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(V2, py2, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(V3, py3, 2*BBE_SIMD_WIDTH);
  }
} /* cdeinterleave4f() */

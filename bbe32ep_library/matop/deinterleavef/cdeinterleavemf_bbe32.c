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

/* any M */
void cdeinterleavemf ( complex_float * restrict y[],
                 const complex_float * restrict x, int M, int N )
{
#if 0
  int n, m;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT((N % 4) == 0); /* multiple of 4 */
  for (m=0; m<M; m++)
  {
    NASSERT_ALIGN(y[m], (2 * BBE_SIMD_WIDTH));
    for (n=0; n<N; n++)
    {
      y[m][n]=x[M*n+m];
    }
  }
#else
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
  const xb_vecNx16 * restrict px2;
  const xb_vecNx16 * restrict px3;
  const int64_t    * restrict px;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
        xb_vecNx16 * restrict py2;
        xb_vecNx16 * restrict py3;
  xb_vecNx16 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  valign alx1, alx2, alx3;
  int Mtail, m, n;

  /* Check restrictions */
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((N & (BBE_SIMD_WIDTH/4-1)) == 0);
  if (M<=0 || N<=0) return;

  Mtail = M & (BBE_SIMD_WIDTH/4 - 1);
  
  /* Process by 4 streams */
  for ( m=0; m<(M>>2); m++ )
  {
    px0 = (const xb_vecNx16 *)(x+m*4);
    py0 = (      xb_vecNx16 *)(y[m*4+0]);
    py1 = (      xb_vecNx16 *)(y[m*4+1]);
    py2 = (      xb_vecNx16 *)(y[m*4+2]);
    py3 = (      xb_vecNx16 *)(y[m*4+3]);

    /* Check restrictions */
    NASSERT_ALIGN(py0, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(py1, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(py2, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(py3, 2*BBE_SIMD_WIDTH);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++ )
    {
      px1 = (const xb_vecNx16 *)((complex_float *)px0+M);
      px2 = (const xb_vecNx16 *)((complex_float *)px1+M);
      px3 = (const xb_vecNx16 *)((complex_float *)px2+M);
      alx1 = BBE_LANX16_PP(px1);
      alx2 = BBE_LANX16_PP(px2);
      alx3 = BBE_LANX16_PP(px3);
      BBE_LVNX16_XP(X0, px0, M*sizeof(complex_float)*(BBE_SIMD_WIDTH/4));
      BBE_LANX16_IP(X1, alx1, px1);
      BBE_LANX16_IP(X2, alx2, px2);
      BBE_LANX16_IP(X3, alx3, px3);

      BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

      BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y1, py1, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y3, py3, 2*BBE_SIMD_WIDTH);
    }
  }

  /* Interleave last (M%4) streams */
  for ( m=(M-Mtail); m<M; m++ )
  {
    py0 = (xb_vecNx16 *)(y[m]);
    px  = (const int64_t *)(x+m);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++ )
    {
#ifdef BBE_LSN_4X64_XP
      xb_vecN_4x64 vtmp64;

      BBE_LSN_4X64_XP(vtmp64, px, M*sizeof(complex_float)); X0 = BBE_MOVNX16_FROMN_4X64(vtmp64);
      BBE_LSN_4X64_XP(vtmp64, px, M*sizeof(complex_float)); X1 = BBE_MOVNX16_FROMN_4X64(vtmp64);
      BBE_LSN_4X64_XP(vtmp64, px, M*sizeof(complex_float)); X2 = BBE_MOVNX16_FROMN_4X64(vtmp64);
      BBE_LSN_4X64_XP(vtmp64, px, M*sizeof(complex_float)); X3 = BBE_MOVNX16_FROMN_4X64(vtmp64);
#else
      X0 = BBE_LV4X16_I(px, 0); px = (const int64_t *)((uintptr_t)px + M*sizeof(complex_float));
      X1 = BBE_LV4X16_I(px, 0); px = (const int64_t *)((uintptr_t)px + M*sizeof(complex_float));
      X2 = BBE_LV4X16_I(px, 0); px = (const int64_t *)((uintptr_t)px + M*sizeof(complex_float));
      X3 = BBE_LV4X16_I(px, 0); px = (const int64_t *)((uintptr_t)px + M*sizeof(complex_float));
#endif

      Y0 = BBE_SELNX16I(X1, X0, BBE_SELI_INTERLEAVE_4_EVEN);
      Y2 = BBE_SELNX16I(X3, X2, BBE_SELI_INTERLEAVE_4_EVEN);
      Y0 = BBE_SELNX16I(Y2, Y0, BBE_SELI_EXTRACT_LO_HALVES);

      BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
    }
  }
#endif
} /* cdeinterleavemf() */

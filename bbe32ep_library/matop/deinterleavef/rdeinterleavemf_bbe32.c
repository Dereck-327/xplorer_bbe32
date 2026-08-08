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
void rdeinterleavemf ( float32_t * restrict y[], 
                 const float32_t * restrict x, int M, int N )
{
#if 0
  int n, m;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT((N % 8) == 0); /* multiple of 8 */
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
  const int16_t    * restrict px;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
        xb_vecNx16 * restrict py2;
        xb_vecNx16 * restrict py3;
        xb_vecNx16 * restrict py4;
        xb_vecNx16 * restrict py5;
        xb_vecNx16 * restrict py6;
        xb_vecNx16 * restrict py7;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  valign alx0, alx1, alx2, alx3;
  int Mtail, m, n;

  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((N & (BBE_SIMD_WIDTH/2-1)) == 0);
  if (M<=0 || N<=0) return;

  Mtail = M & (BBE_SIMD_WIDTH/2 - 1);

  /* Process by 8 streams */
  for ( m=0; m<(M>>3); m++ )
  {
    py0 = (xb_vecNx16 *)(y[m*8+0]);
    py1 = (xb_vecNx16 *)(y[m*8+1]);
    py2 = (xb_vecNx16 *)(y[m*8+2]);
    py3 = (xb_vecNx16 *)(y[m*8+3]);
    py4 = (xb_vecNx16 *)(y[m*8+4]);
    py5 = (xb_vecNx16 *)(y[m*8+5]);
    py6 = (xb_vecNx16 *)(y[m*8+6]);
    py7 = (xb_vecNx16 *)(y[m*8+7]);
    px  = (const int16_t *)(x+m*8);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++ )
    {
      px0 = (const xb_vecNx16 *)(px);
      px1 = (const xb_vecNx16 *)((float32_t *)px0+M);
      px2 = (const xb_vecNx16 *)((float32_t *)px1+M);
      px3 = (const xb_vecNx16 *)((float32_t *)px2+M);
      BBE_LVNX16_XP(X0, px0, sizeof(float32_t)*M*(BBE_SIMD_WIDTH/4));
      alx1 = BBE_LANX16_PP(px1);
      alx2 = BBE_LANX16_PP(px2);
      alx3 = BBE_LANX16_PP(px3);
      BBE_LANX16_IP(X1, alx1, px1);
      BBE_LANX16_IP(X2, alx2, px2);
      BBE_LANX16_IP(X3, alx3, px3);

      px = (const int16_t *)((intptr_t)px0 + sizeof(float32_t)*M*(BBE_SIMD_WIDTH/4));
      px1 = (xb_vecNx16 *)((float32_t *)px0+M);
      px2 = (xb_vecNx16 *)((float32_t *)px1+M);
      px3 = (xb_vecNx16 *)((float32_t *)px2+M);
      alx0 = BBE_LANX16_PP(px0);
      alx1 = BBE_LANX16_PP(px1);
      alx2 = BBE_LANX16_PP(px2);
      alx3 = BBE_LANX16_PP(px3);
      BBE_LANX16_IP(X4, alx0, px0);
      BBE_LANX16_IP(X5, alx1, px1);
      BBE_LANX16_IP(X6, alx2, px2);
      BBE_LANX16_IP(X7, alx3, px3);

      BBE_DSELNX16I(X4, X0, X4, X0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X5, X1, X5, X1, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X6, X2, X6, X2, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X7, X3, X7, X3, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y5, Y4, X5, X4, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELNX16I(Y7, Y6, X7, X6, BBE_DSELI_INTERLEAVE_2);

      BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y1, py1, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y3, py3, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y4, py4, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y5, py5, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y6, py6, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(Y7, py7, 2*BBE_SIMD_WIDTH);
    }
  }

  /* Interleave last (M%8) streams */
  for ( m=(M-Mtail); m<M; m++ )
  {
    px  = (const int16_t *)(x+m);
    py0 = (xb_vecNx16 *)(y[m]);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++ )
    {
      BBE_LPNX16_XP(X0, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X1, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X2, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X3, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X4, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X5, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X6, px, M*sizeof(float32_t));
      BBE_LPNX16_XP(X7, px, M*sizeof(float32_t));

      Y0 = BBE_SELNX16I(X1, X0, BBE_SELI_INTERLEAVE_2_EVEN);
      Y2 = BBE_SELNX16I(X3, X2, BBE_SELI_INTERLEAVE_2_EVEN);
      Y4 = BBE_SELNX16I(X5, X4, BBE_SELI_INTERLEAVE_2_EVEN);
      Y6 = BBE_SELNX16I(X7, X6, BBE_SELI_INTERLEAVE_2_EVEN);
      Y0 = BBE_SELNX16I(Y2, Y0, BBE_SELI_INTERLEAVE_4_EVEN);
      Y4 = BBE_SELNX16I(Y6, Y4, BBE_SELI_INTERLEAVE_4_EVEN);
      Y0 = BBE_SELNX16I(Y4, Y0, BBE_SELI_EXTRACT_LO_HALVES);

      BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
    }
  }
#endif
} /* rdeinterleavemf() */

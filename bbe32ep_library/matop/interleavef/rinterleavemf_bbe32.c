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
    M-to-1 complex/real streams interleave
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
M-to-1 complex/real streams interleave

Description: interleave element by element M=2,3,4 or 8 streams into a
single stream. Use rinterleave<M>() functions for real data, and
cinterleave<M>() functions - for complex data.

Representation:
<r|c>interleave<M>   16-bit fixed-point data
<r|c>interleave<M>f  IEEE-754 Std single precision floating-point data

Parameters:
Input:
M                    Number of streams
N                    Number of elements per each input stream
x0[N],...,x<M-1>[N]  Input data streams
x[M]                 M pointers to input data streams
Output:
y[M*N]               Interleaved data streams

Restrictions:
x0,...,x<M-1>,
x[0..M-1],y          Must not overlap and must be aligned on 32-byte boundary
N                    Must be a multiple of:
                       16 for real fixed-point data 
                        8 for complex fixed-point data and real floating-point data
                        4 for complex floating-point data
-------------------------------------------------------------------------*/

/* any M */
void rinterleavemf ( float32_t * restrict y,
                     float32_t * restrict x[], int M, int N )
{
#if 0
  int n, m;
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT((N%8)==0); /* multiple of 8 */
  for (m=0; m<M; m++)
  {
      NASSERT_ALIGN(x[m], (2 * BBE_SIMD_WIDTH));
      for (n=0; n<N; n++)
      {
          y[M*n+m]=x[m][n];
      }
  }
#else
  const xb_vecNx16 * restrict px0;
  const xb_vecNx16 * restrict px1;
  const xb_vecNx16 * restrict px2;
  const xb_vecNx16 * restrict px3;
  const xb_vecNx16 * restrict px4;
  const xb_vecNx16 * restrict px5;
  const xb_vecNx16 * restrict px6;
  const xb_vecNx16 * restrict px7;
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
        xb_vecNx16 * restrict py2;
        xb_vecNx16 * restrict py3;
        int16_t    * restrict py;
  xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
  valign aly0, aly1, aly2, aly3;
  int Mtail, m, n;

  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((N & (BBE_SIMD_WIDTH/2-1)) == 0);
  if (M<=0 || N<=0) return;

  Mtail = M & (BBE_SIMD_WIDTH/2 - 1);
  aly0 = aly1 = aly2 = aly3 = BBE_ZALIGN();

  /* Process by 8 streams */
  for ( m=0; m<(M>>3); m++ )
  {
    px0 = (const xb_vecNx16 *)(x[m*8+0]);
    px1 = (const xb_vecNx16 *)(x[m*8+1]);
    px2 = (const xb_vecNx16 *)(x[m*8+2]);
    px3 = (const xb_vecNx16 *)(x[m*8+3]);
    px4 = (const xb_vecNx16 *)(x[m*8+4]);
    px5 = (const xb_vecNx16 *)(x[m*8+5]);
    px6 = (const xb_vecNx16 *)(x[m*8+6]);
    px7 = (const xb_vecNx16 *)(x[m*8+7]);
    py  = (int16_t *)(y+m*8);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++ )
    {
      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1, px1, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X2, px2, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X3, px3, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X4, px4, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X5, px5, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X6, px6, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X7, px7, 2*BBE_SIMD_WIDTH);

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

      py0 = (xb_vecNx16 *)(py);
      py1 = (xb_vecNx16 *)((float32_t *)py0+M);
      py2 = (xb_vecNx16 *)((float32_t *)py1+M);
      py3 = (xb_vecNx16 *)((float32_t *)py2+M);
      BBE_SVNX16_XP(Y0, py0, sizeof(float32_t)*M*(BBE_SIMD_WIDTH/4));
      BBE_SANX16_IP(Y1, aly1, py1);
      BBE_SANX16_IP(Y2, aly2, py2);
      BBE_SANX16_IP(Y3, aly3, py3);
      BBE_SANX16POS_FP(aly1, py1);
      BBE_SANX16POS_FP(aly2, py2);
      BBE_SANX16POS_FP(aly3, py3);

      py = (int16_t *)((intptr_t)py0 + sizeof(float32_t)*M*(BBE_SIMD_WIDTH/4));
      py1 = (xb_vecNx16 *)((float32_t *)py0+M);
      py2 = (xb_vecNx16 *)((float32_t *)py1+M);
      py3 = (xb_vecNx16 *)((float32_t *)py2+M);
      BBE_SANX16_IP(Y4, aly0, py0);
      BBE_SANX16_IP(Y5, aly1, py1);
      BBE_SANX16_IP(Y6, aly2, py2);
      BBE_SANX16_IP(Y7, aly3, py3);
      BBE_SANX16POS_FP(aly0, py0);
      BBE_SANX16POS_FP(aly1, py1);
      BBE_SANX16POS_FP(aly2, py2);
      BBE_SANX16POS_FP(aly3, py3);
    }
  }

  /* Interleave last (M%8) streams */
  for ( m=(M-Mtail); m<M; m++ )
  {
    px0 = (const xb_vecNx16 *)(x[m]);
    py  = (int16_t *)(y+m);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++ )
    {
      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);

      Y0 = BBE_REPNX16C(X0, 0);
      Y1 = BBE_REPNX16C(X0, 1);
      Y2 = BBE_REPNX16C(X0, 2);
      Y3 = BBE_REPNX16C(X0, 3);
      Y4 = BBE_REPNX16C(X0, 4);
      Y5 = BBE_REPNX16C(X0, 5);
      Y6 = BBE_REPNX16C(X0, 6);
      Y7 = BBE_REPNX16C(X0, 7);

      BBE_SPNX16_XP(Y0, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y1, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y2, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y3, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y4, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y5, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y6, py, M*sizeof(float32_t));
      BBE_SPNX16_XP(Y7, py, M*sizeof(float32_t));
    }
  }
#endif
} /* rinterleavmf() */

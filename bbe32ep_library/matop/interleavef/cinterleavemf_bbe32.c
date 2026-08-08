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
void cinterleavemf ( complex_float * restrict y,
                     complex_float * restrict x[], int M, int N )
{
#if 0
  int n, m;
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT((N%4)==0); /* multiple of 4 */
  for (m=0; m<M; m++)
  {
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
        xb_vecNx16 * restrict py0;
        xb_vecNx16 * restrict py1;
        xb_vecNx16 * restrict py2;
        xb_vecNx16 * restrict py3;
        int64_t    * restrict py;
  xb_vecNx16 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
  valign aly1, aly2, aly3;
  int Mtail, m, n;

  /* Check restrictions */
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT((N & (BBE_SIMD_WIDTH/4-1)) == 0);
  if (M<=0 || N<=0) return;

  Mtail = M & (BBE_SIMD_WIDTH/4 - 1);
  aly1 = aly2 = aly3 = BBE_ZALIGN();
  
  /* Process by 4 streams */
  for ( m=0; m<(M>>2); m++ )
  {
    px0 = (const xb_vecNx16 *)(x[m*4+0]);
    px1 = (const xb_vecNx16 *)(x[m*4+1]);
    px2 = (const xb_vecNx16 *)(x[m*4+2]);
    px3 = (const xb_vecNx16 *)(x[m*4+3]);
    py0 = (      xb_vecNx16 *)(y+m*4);

    /* Check restrictions */
    NASSERT_ALIGN(px0, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(px1, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(px2, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(px3, 2*BBE_SIMD_WIDTH);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++ )
    {
      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1, px1, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X2, px2, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X3, px3, 2*BBE_SIMD_WIDTH);

      BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_4);

      py1 = (xb_vecNx16 *)((complex_float *)py0+M);
      py2 = (xb_vecNx16 *)((complex_float *)py1+M);
      py3 = (xb_vecNx16 *)((complex_float *)py2+M);
      BBE_SVNX16_XP(Y0, py0, M*sizeof(complex_float)*(BBE_SIMD_WIDTH/4));
      BBE_SANX16_IP(Y1, aly1, py1);
      BBE_SANX16_IP(Y2, aly2, py2);
      BBE_SANX16_IP(Y3, aly3, py3);
      BBE_SANX16POS_FP(aly1, py1);
      BBE_SANX16POS_FP(aly2, py2);
      BBE_SANX16POS_FP(aly3, py3);
    }
  }

  /* Interleave last (M%4) streams */
  for ( m=(M-Mtail); m<M; m++ )
  {
    px0 = (const xb_vecNx16 *)(x[m]);
    py  = (int64_t *)(y+m);

    __Pragma("loop_count min=1");
    for ( n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++ )
    {
#ifdef BBE_SSN_4X64_XP
      xb_vecN_4x64 vtmp64, Y0_64, Y1_64, Y2_64, Y3_64;

      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);
      vtmp64 = BBE_MOVN_4X64_FROMNX16(X0);

      Y0_64 = BBE_REPN_4X64(vtmp64, 0);
      Y1_64 = BBE_REPN_4X64(vtmp64, 1);
      Y2_64 = BBE_REPN_4X64(vtmp64, 2);
      Y3_64 = BBE_REPN_4X64(vtmp64, 3);
      BBE_SSN_4X64_XP(Y0_64, py, M*sizeof(complex_float));
      BBE_SSN_4X64_XP(Y1_64, py, M*sizeof(complex_float));
      BBE_SSN_4X64_XP(Y2_64, py, M*sizeof(complex_float));
      BBE_SSN_4X64_XP(Y3_64, py, M*sizeof(complex_float));
#else
      BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);

      Y0 = BBE_SHFLNX16I(X0, BBE_SHFLI_REP_0X4);
      Y1 = BBE_SHFLNX16I(X0, BBE_SHFLI_REP_1X4);
      Y2 = BBE_SHFLNX16I(X0, BBE_SHFLI_REP_2X4);
      Y3 = BBE_SHFLNX16I(X0, BBE_SHFLI_REP_3X4);
      BBE_SV4X16_I(Y0, py, 0); py = (int64_t *)((uintptr_t)py + M*sizeof(complex_float));
      BBE_SV4X16_I(Y1, py, 0); py = (int64_t *)((uintptr_t)py + M*sizeof(complex_float));
      BBE_SV4X16_I(Y2, py, 0); py = (int64_t *)((uintptr_t)py + M*sizeof(complex_float));
      BBE_SV4X16_I(Y3, py, 0); py = (int64_t *)((uintptr_t)py + M*sizeof(complex_float));
#endif
    }
  }
#endif
} /* cinterleavemf() */

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
#include "bs_common.h"
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

/* M=3 */
void rinterleave3 ( int16_t * restrict y, const int16_t * restrict x0,
                                          const int16_t * restrict x1,
                                          const int16_t * restrict x2, int N )
{
  const xb_vecNx16* restrict w0 = (const xb_vecNx16*)x0;
  const xb_vecNx16* restrict w1 = (const xb_vecNx16*)x1;
  const xb_vecNx16* restrict w2 = (const xb_vecNx16*)x2;
  xb_vecNx16* restrict z = (xb_vecNx16*)y;
  xb_vecNx16 X0, X1, X2, Y0, Y1, Z0, Z1, Z2;
  int n;
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x0, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x1, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x2, 2 * BBE_SIMD_WIDTH);
  NASSERT((N&(BBE_SIMD_WIDTH / 2 - 1)) == 0);

  for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
  {
    BBE_LVNX16_IP(X0, w0, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, w1, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, w2, 2 * BBE_SIMD_WIDTH);
    /* convert real data in such a way to apply complex interleaving process */
    Z0 = BBE_SELNX16I(X1, X0, BBE_SELI_INTERLEAVE_1_EVEN);   /*0,32,2,34,4,36,....*/
    Z1 = BBE_SELNX16I(X0, X2, BBE_SELI_INTERLEAVE_1_EVENODD);/*64,1,66,3...       */
    X2 = BBE_SELNX16I(X2, X1, BBE_SELI_INTERLEAVE_1_ODD);    /*33,65,35,37...     */
    X0 = Z0; X1 = Z1;
    /* complex interleaving */
    BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_2);
    /*Y0=X(0,1,32,33,2,3,34,35,4,5,36,37,6,7,38,39,8,9,40,41,10,11,42,43,12,13,44,45,14,15,46,47)
      Y1=X(16,17,48,49,18,19,50,51,20,21,52,53,22,23,54,55,24,25,56,57,26,27,58,59,28,29,60,61,30,31,62,63)*/
    BBE_DSELNX16I(Z1, Z0, X2, Y0, 0);
    /*Z0=0,1,32,33,64,65,2,3,34,35,66,67,4,5,36,37,68,69,6,7,38,39,70,71,8,9,40,41,72,73,10,11
      Z1=42,43,74,75,12,13,44,45,76,77,14,15,46,47,78,79,0,1,32,33,80,81,2,3,34,35,82,83,4,5,36,37*/
    BBE_DSELNX16I_H(Z1, Z2, X2, Y1, 0);
    /*Z2=84,85,22,23,54,55,86,87,24,25,56,57,88,89,26,27,58,59,90,91,28,29,60,61,92,93,30,31,62,63,94,95
      Z1=  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,  ,16,17,48,49,80,81,18,19,50,51,82,83,20,21,52,53*/

    BBE_SVNX16_IP(Z0, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Z1, z, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Z2, z, 2 * BBE_SIMD_WIDTH);
  }
} /* rinterleave3() */

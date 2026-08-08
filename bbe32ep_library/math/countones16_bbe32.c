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
  NatureDSP_Baseband library. Math functions
    Count One Bits in a Word
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* Tables for 16-bit fixed-point countones(x) functions. */
#include "countones_16b_tbl.h"
/*-------------------------------------------------------------------------
Count One Bits in a Word

Description: Functions count the number of one bits in the number or each number of a vector.

Data format: 16-bit/32-bit fixed-point format

Accuracy: exact

Parameters:
Input:
x[N]   Input data, 16-bit/32-bit
N      Length of input/output data vectors
Output:
z[N]   Results, 16-bit/32-bit

Restrictions:
z,x,y  Aligned on 32-byte boundary
z,x,y  Must not overlap
N      Multiple of 16 (vcountones16) or 8 (vcountones32)
-------------------------------------------------------------------------*/

int countones16 (int16_t x)
{
  /*int16_t tab[16] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
  int n;
  uint16_t x0,x1,x2,x3;
  for (n=0;n<N;n++)
  {
  x0 = x[n] >> 12;
  x1 = x[n] >> 8;
  x2 = x[n] >> 4;
  x3 = x[n];
  x0 &= 15; x1 &= 15;
  x2 &= 15; x3 &= 15;
  z[n] = tab[x0] + tab[x1] + tab[x2] + tab[x3];*/
  uint16_t u0, u1, u2, u3;
  int32_t z0;
  u0 = (uint16_t)x;
  u1 = XT_SRLI(u0, 12);
  u2 = XT_SRLI(u0, 8);
  u3 = XT_SRLI(u0, 4);
  u0 = XT_AND(u0, 15);
  u1 = XT_AND(u1, 15);
  u2 = XT_AND(u2, 15);
  u3 = XT_AND(u3, 15);
  z0 = (int16_t)countones_16b_tbl[u0] + countones_16b_tbl[u1] + countones_16b_tbl[u2] + countones_16b_tbl[u3];
  return z0;
} /* countones16() */

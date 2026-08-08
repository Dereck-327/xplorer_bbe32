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
    Complex Matrix Conjugate Transpose
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, cmattran2x2nf, ( complex_float * restrict y, 
                             const complex_float * restrict x, 
                                   int L ))
#else

#ifndef BBE_SHFLN_4XCF32
#define BBE_SHFLN_4XCF32(a, b) (BBE_MOVN_4XCF32_FROMNX16(BBE_SHFLNX16(BBE_MOVNX16_FROMN_4XCF32(a), (b) )))
#endif

/*-------------------------------------------------------------------------
Complex Matrix Conjugate Transpose

Description: These functions perform transposition and then take the complex
conjugate for each matrix from input sequence. Results are stored to output
sequence. Both the block order and streaming order are allowed for input/output
matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
1. Complex conjugation of fixed-point data may involve 16-bit saturation of
   imaginary components
2. The functions cmattrannxnn(), cmattrannxmn() and cmattrannxmnf() (conjugate
   transpose for the block order) may distort the input matrices sequence x[L*S].

Parameters:
Input:
x[L*S]  Sequence of input matrices.
N,M     Matrix dimensions 
L       Number of matrices
Output:
y[L*S]  Sequence of output matrices

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, 2x2->2x2, S=4
   Restrictions:
    None
*/
void cmattran2x2nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L )
{
  static const uint16_t ALIGN(32) tblSel[BBE_SIMD_WIDTH]=
  {
    0x0, 0x1, 0x2, 0x3, 0x8, 0x9, 0xA, 0xB, 0x4, 0x5, 0x6, 0x7, 0xC, 0xD, 0xE, 0xF
  };
  const xb_vecN_4xcf32 * restrict px;
        xb_vecN_4xcf32 * restrict py;
  int l;

  xb_vecN_4xcf32 X, Y;
  xb_vecNx16 vTmp;
  vselN sel;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel, 0);
  sel = BBE_MOVVSELNX16(vTmp,0);
  px = (const xb_vecN_4xcf32 *)x;
  py = (      xb_vecN_4xcf32 *)y;

  for (l=0; l<L; l++)
  {
    /* Load input matrix X */
    BBE_LVN_4XCF32_IP(X, px, 2*BBE_SIMD_WIDTH);

    /* Transpose and conjugate */
    Y = BBE_SHFLN_4XCF32(X, sel);
    Y = BBE_CONJN_4XCF32(Y);

    /* Save results */
    BBE_SVN_4XCF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
  }
} /* cmattran2x2nf() */
#endif

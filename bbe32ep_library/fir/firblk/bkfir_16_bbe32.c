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
    NatureDSP_Baseband library. FIR filters and Related Functions
    Block Real FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "bkfir_common.h"

/*-------------------------------------------------------------------------
Block Real FIR Filter

Computes a real FIR filter (direct-form) using IR stored in vector h. The
real data input is stored in vector x. The filter output result is stored
in vector y. The filter calculates N output samples using M coefficients
and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by bkfir[f]_algDelay(M).

Representation:
bkfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
bkfirf  IEEE-754 Std. single precision floating-point format for filter 
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[M]    Filter coefficients; h[0] is to be multiplied by the newest sample
N       Length of sample block
M       Length of filter
x[N]    Input samples
Output:
y[N]    Output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (bkfir) or 8 (bkfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, bkfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/*processing function for M==16*/
void bkfir_process_16(int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N)
{

  int n;
  const xb_vecNx16 * restrict X;
        xb_vecNx16 * restrict Y;
        xb_vecNx16 * restrict D;

  xb_vecNx16 d0, d1;
  xb_vecNx16 x0, x1, y0;
  uint32_t h0, h1, h2, h3, h4, h5, h6, h7;
  xb_vecNx40 w0; 

  NASSERT(N>0 && !(N&15));
  NASSERT(M==16);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(delay);
  NASSERT_ALIGN32(h);

  X = (const xb_vecNx16 *) x;
  Y = (      xb_vecNx16 *) y;
  D = (      xb_vecNx16 *) delay;

  h0 = XT_L32I((int*)h, 0);
  h1 = XT_L32I((int*)h, 4);
  h2 = XT_L32I((int*)h, 8);
  h3 = XT_L32I((int*)h, 12);
  h4 = XT_L32I((int*)h, 16);
  h5 = XT_L32I((int*)h, 20); 
  h6 = XT_L32I((int*)h, 24);
  h7 = XT_L32I((int*)h, 28);

  __Pragma("ymemory(X)")
  __Pragma("ymemory(D)")
  __Pragma("loop_count min=1")
  for (n=0; n<(N/BBE_SIMD_WIDTH); n++)
  {
    BBE_LVNX16_IP(d0, D, 0 * 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(d1, X, 1 * 2 * BBE_SIMD_WIDTH);
   
    BBE_SELPRNX16I(x1, x0, d1, d0, 1); w0 = BBE_MULNX16PR(x0, x1, h0);
    BBE_SELPRNX16I(x1, x0, d1, d0, 3); BBE_MULANX16PR(w0, x0, x1, h1);
    BBE_SELPRNX16I(x1, x0, d1, d0, 5); BBE_MULANX16PR(w0, x0, x1, h2);
    BBE_SELPRNX16I(x1, x0, d1, d0, 7); BBE_MULANX16PR(w0, x0, x1, h3);
    BBE_SELPRNX16I(x1, x0, d1, d0, 9); BBE_MULANX16PR(w0, x0, x1, h4);
    BBE_SELPRNX16I(x1, x0, d1, d0, 11); BBE_MULANX16PR(w0, x0, x1, h5);
    BBE_SELPRNX16I(x1, x0, d1, d0, 13); BBE_MULANX16PR(w0, x0, x1, h6);
    BBE_SELPRNX16I(x1, x0, d1, d0, 15); BBE_MULANX16PR(w0, x0, x1, h7);

    //Save
    y0 = BBE_PACKQNX40(w0);
    BBE_SVNX16_IP(y0, Y, 1 * 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(d1, D, 0 * 2 * BBE_SIMD_WIDTH);
    //BBE_SVNX16_I(d1, D, 0 * 2 * BBE_SIMD_WIDTH);
  }

}

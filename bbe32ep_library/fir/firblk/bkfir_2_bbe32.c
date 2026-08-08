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

/*processing function for M==2*/
void bkfir_process_2(int16_t* restrict delay, int16_t* restrict y, const int16_t * restrict x, const int16_t* restrict h, int M, int N)
{
  int n;
  const xb_vecNx16 * restrict pX = (const xb_vecNx16 *) x;
        xb_vecNx16 * restrict pY = (      xb_vecNx16 *) y;
        xb_vecNx16 * restrict pD = (      xb_vecNx16 *)(delay);
  const int        * restrict pH = (const int        *) h;

  xb_vecNx16 d0, d1;
  xb_vecNx16 x0, x1, t0;
  uint32_t h0;
  xb_vecNx40 A0; 

  NASSERT(N>0 && !(N&15));
  NASSERT(M==2);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(delay);
  NASSERT_ALIGN32(h);

  h0 = XT_L32I(pH,0);
  d1 = BBE_LVNX16_I(pD, 0);

  __Pragma("ymemory(pX)");
  __Pragma("loop_count min=1");
  for (n=0; n<(N/BBE_SIMD_WIDTH); n++)
  {
    d0 = d1;
    BBE_LVNX16_IP(d1, pX, 2*BBE_SIMD_WIDTH);
    
    BBE_SELPRNX16I(x1,x0,d1,d0,15);
    A0 = BBE_MULNX16PR(x0,x1,h0);
  
    t0 = BBE_PACKQNX40(A0);
    BBE_SVNX16_IP(t0, pY, 2*BBE_SIMD_WIDTH);
  }
  BBE_SVNX16_I(d1, pD, 0);
}

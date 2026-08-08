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
    Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "cxfir_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR Filter

Computes a complex FIR filter (direct-form) using complex IR stored in 
vector h. The complex data input is stored in vector x. The filter output
result is stored in vector y. The filter calculates N output samples using
M coefficients and requires last M+N-1 samples in the delay line. 

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by cxfir[f]_algDelay(M).

Representation:
cxfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
cxfirf  IEEE-754 Std. single precision floating-point format for filter
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[2*M]  Filter coefficients; h[0]+j*h[1] is to be multiplied by the newest
        sample,
N       Length of sample block
M       Length of filter
x[N]    Complex input samples
Output:
y[N]    Complex output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cxfir) or 4 (cxfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, cxfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/
/*processing function for M==2*/
void cxfir_process_2(complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N)
{
  int n;
  const xb_vecNx16 *          C;
  const xb_vecNx16 *          X;
        xb_vecNx16 * restrict Y;
  
  xb_vecNx16 d0, d1;
  xb_vecNx16 p0;
  xb_vecNx16 c0, c1;
  xb_vecNx40 w0;
  xb_vecNx16 y0;

  NASSERT(N>0 && !(N&7)); 
  NASSERT(M==2);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(delayLine);
  NASSERT_ALIGN32(coef);

  C = (const xb_vecNx16*)coef;
  X = (const xb_vecNx16*)x;
  Y = (      xb_vecNx16*)y;
  
  //Load delay line
  d0 = BBE_LVNX16_I( (const xb_vecNx16*)delayLine, 0*4*BBE_SIMD_WIDTH/2 );
  d0 = BBE_REPNX16C(d0, 0);

  // Load filter coefficients
  BBE_LPNX16_IP(c0, C, 4);
  BBE_LPNX16_IP(c1, C, 4);

  c0 = BBE_REPNX16C(c0, 0);
  c1 = BBE_REPNX16C(c1, 0);
  
  __Pragma("ymemory(X)")
  __Pragma("ymemory(Y)")
  __Pragma("loop_count min=1")
  for ( n=0; n<N/(BBE_SIMD_WIDTH/2); n++ )
  {
      //Load input matrix  
      BBE_LVNX16_IP(d1, X, 4*BBE_SIMD_WIDTH/2);

      p0 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_LEFT_2);
      d0 = d1;

      //Multiple
      w0 = BBE_MULNX16C(p0, c0);
      BBE_MULANX16C(w0, d1, c1);

      //Pack and save result
      y0 = BBE_PACKQNX40(w0);
      BBE_SVNX16_IP(y0, Y, +4*BBE_SIMD_WIDTH/2);
  }

  //Save delay line
  d0 = BBE_REPNX16C(d0, BBE_SIMD_WIDTH/2-1);
  BBE_SPNX16_I(d0, delayLine, 0);
  
} //cxfir_process_2()

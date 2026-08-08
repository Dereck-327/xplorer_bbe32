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

/*processing function for M==8*/
void cxfir_process_8(complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N)
{

  int n;
  const xb_vecNx16 * restrict X = (const xb_vecNx16 *) x;
        xb_vecNx16 * restrict Y = (      xb_vecNx16 *) y;
        xb_vecNx16 * restrict D;
  const xb_vecNx16 * restrict H;

  xb_vecNx16 d0, d1, c0;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
  xb_vecNx16 t0;
  xb_vecNx16 h0, h1, h2, h3, h4, h5, h6, h7;
  xb_vecNx40 w0; 

  NASSERT(N>0 && !(N&7)); 
  NASSERT(M==8);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);
  NASSERT_ALIGN32(delayLine);
  NASSERT_ALIGN32(coef);

  D = (xb_vecNx16 *)(delayLine);

  d1 = BBE_LVNX16_I(D, 0);

  __Pragma("ymemory(X)")
  __Pragma("ymemory(Y)")
  __Pragma("ymemory(D)")
  __Pragma("loop_count min=1")
  for (n=0; n<N/(BBE_SIMD_WIDTH/2); n++)
  {
    H = (const xb_vecNx16 *)coef;    

    d0 = d1;

    BBE_LVNX16_IP(d1, X, 2 * BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(c0, H, 2 * BBE_SIMD_WIDTH);
    h0 = BBE_REPNX16C(c0, 0);
    h1 = BBE_REPNX16C(c0, 1);
    h2 = BBE_REPNX16C(c0, 2);
    h3 = BBE_REPNX16C(c0, 3);
    h4 = BBE_REPNX16C(c0, 4);
    h5 = BBE_REPNX16C(c0, 5);
    h6 = BBE_REPNX16C(c0, 6);
    h7 = BBE_REPNX16C(c0, 7);

    x0 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_2);
    x1 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_4);
    x2 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_6);
    x3 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_8);
    x4 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_10);
    x5 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_12);
    x6 = BBE_SELNX16I(d1, d0, BBE_SELI_ROTATE_RIGHT_14);
    x7  = d1;

    w0 = BBE_MULNX16C(x0, h0);
    BBE_MULANX16C(w0, x1, h1);
    BBE_MULANX16C(w0, x2, h2);
    BBE_MULANX16C(w0, x3, h3);
    BBE_MULANX16C(w0, x4, h4);
    BBE_MULANX16C(w0, x5, h5);
    BBE_MULANX16C(w0, x6, h6);
    BBE_MULANX16C(w0, x7, h7);

    // Save
    t0 = BBE_PACKQNX40(w0);
    BBE_SVNX16_IP(t0, Y, 2 * BBE_SIMD_WIDTH);
  }

  // Save delayLine after calculate
  BBE_SVNX16_IP( d1, D, 2*BBE_SIMD_WIDTH);

} //cxfir_process_8()

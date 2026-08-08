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
    Decimating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "firdec_common.h"

/*-------------------------------------------------------------------------
Decimating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with decimation using real IR 
stored in vector h. The complex data input is stored in vector x. The filter
output result is stored in vector y. The filter calculates N output samples
using M coefficients and requires last D*N+M-1 samples in the delay line.

NOTE:
To avoid aliasing, the IR should be synthesized in such a way that filter pass
band is limited by input sample frequency divided by 2*D.

Representation:
firdec   16-bit signed fixed-point format
         Filter coefficients are Q15
         Number of fractional bits for input/output samples is user-difined
firdecf  IEEE-754 Std. single precision floating-point format for filter 
         coefficients and input/output samples

Parameters:
Input:
D        Decimation factor
N        Length of output sample block
M        Length of filter
h[M]     Filter coefficients; h[0] is to be multiplied by the newest 
         sample
x[N*D]   Input complex samples
Output:
y[N]     Output complex samples

Restrictions:
x,y      Must not overlap
x,y      Aligned on 32-byte boundary
N        Multiple of 8 (firdec) or 4 (firdecf)
M        2,4,8 or a positive multiple of 16 for D=2,3,4; or 
         a positive multiple of 16 for D>4
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
filter lengths M=2,4,8,16 and 32 and decimation factors D=2,3 and 4, in
any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firdec[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* processing function D==4, M==8 */
static void firdec_proc_D4_M8( int16_t * restrict y,
                        const int16_t * restrict x,
                        const int16_t * restrict coef,
                              int16_t * restrict delayLine,
                        int M, int N, int D )
{

        xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X;

  xb_vecNx40 w0;

  xb_vecNx16 cf;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 d0, d1, d2, d3;
  uint32_t   q0, q1, q2, q3;
  xb_vecNx16 y0;

  xb_vecNx16 p00, p01, p10, p11;
  xb_vecNx16 p20, p21, p30, p31;

  int n;

  NASSERT( N>0  && !(N&7) );
  NASSERT( D==4 && M==8 );
  NASSERT_ALIGN32( y         );
  NASSERT_ALIGN32( x         );
  NASSERT_ALIGN32( coef      );
  NASSERT_ALIGN32( delayLine );

  Y = (      xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;

  //
  // Load the delay line state.
  //

  x0 = BBE_LVNX16_I((const xb_vecNx16*)delayLine, 0);

  d0 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_LEFT_12);
  d1 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_LEFT_8);
  d2 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_LEFT_4);
  d3 = x0;

  //
  // Process data.
  //

  cf = BBE_LVNX16_I((const xb_vecNx16*)coef, 0);

  for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
  {
      // Load 8x4 input samples, CQ15
      BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(x1, X, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(x2, X, 2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(x3, X, 2 * BBE_SIMD_WIDTH);

      // Transposition 8x4 -> 4x8 (4 banks, each of 8 samples).
      BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);

      //
      // Coefficients bank 3
      //

      BBE_SELPCNX16I(p01, p00, x0, d0, 7);

      d0 = x0;

      // Load 4 real filter coefficients, Q15.
      q0 = BBE_EXTRNX16C(cf, 3);

      w0 = BBE_MULNX16PR(p01, p00, q0);

      //                            
      // Coefficients bank 0
      //                            

      BBE_SELPCNX16I(p11, p10, x1, d1, 6);

      d1 = x1;

      q1 = BBE_EXTRNX16C(cf, 0);

      BBE_MULANX16PR(w0, p11, p10, q1);

      //
      // Coefficients bank 1
      //

      BBE_SELPCNX16I(p21, p20, x2, d2, 6);

      d2 = x2;

      q2 = BBE_EXTRNX16C(cf, 1);

      BBE_MULANX16PR(w0, p21, p20, q2);

      //
      // Coefficients bank 2
      //

      BBE_SELPCNX16I(p31, p30, x3, d3, 6);

      d3 = x3;

      q3 = BBE_EXTRNX16C(cf, 2);

      BBE_MULANX16PR(w0, p31, p30, q3);

      //
      // Save 8 output samples.
      //

      // CQ15 <- CQ30 - 15 w/ rounding and saturation.
      y0 = BBE_PACKQNX40(w0);

      BBE_SVNX16_IP(y0, Y, 2 * BBE_SIMD_WIDTH);
  }

  //
  // Save the delay line state.
  //
  
  d0 = BBE_SELNX16I(d1, d0, BBE_SELI_INTERLEAVE_4_HI);
  d1 = BBE_SELNX16I(d3, d2, BBE_SELI_INTERLEAVE_4_HI);
  d0 = BBE_SELNX16I(d1, d0, BBE_SELI_EXTRACT_HI_HALVES);
  BBE_SVNX16_I(d0, (xb_vecNx16*)delayLine, 0 * 2 * BBE_SIMD_WIDTH);
  //d0 = BBE_SELNX16I(d0, d0, BBE_SELI_ROTATE_RIGHT_12);
  //d1 = BBE_SELNX16I(d1, d1, BBE_SELI_ROTATE_RIGHT_12);
  //d2 = BBE_SELNX16I(d2, d2, BBE_SELI_ROTATE_RIGHT_12);
  //d3 = BBE_SELNX16I(d3, d3, BBE_SELI_ROTATE_RIGHT_12);

  //BBE_SV4X16_I(d0, delayLine, 0 * 2 * 4);
  //BBE_SV4X16_I(d1, delayLine, 1 * 2 * 4);
  //BBE_SV4X16_I(d2, delayLine, 2 * 2 * 4);
  //BBE_SV4X16_I(d3, delayLine, 3 * 2 * 4);
} // firdec_proc_D4_M8()

const tFirFxdxns firdec_4d_8_8n  ={&firdec_alloc_gen,firdec_proc_D4_M8 };

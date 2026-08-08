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

/* processing function D==3, M==8 */
static void firdec_proc_D3_M8 (  int16_t * restrict y,
                          const int16_t * restrict x,
                          const int16_t * restrict coef,
                                int16_t * restrict delayLine,
                            int M, int N, int D )
{

        xb_vecNx16 * restrict Y;
  const xb_vecNx16 *          X;
        xb_vecNx16 *          S;

  xb_vecNx40 w0;

  xb_vecNx16 cf;
  xb_vecNx16 x0, x1, x2;
  xb_vecNx16 s0, s1, s2;
  xb_vecNx16 d0, d1, d2;
  xb_vecNx16 y0;
  xb_vecNx16 p00, p01, p02;
  xb_vecNx16 p10, p11;
  xb_vecNx16 p20, p21, p22;
  uint32_t   q0, q1, q2, q3;

  int n;

  NASSERT( N>0  && N%8==0 );
  NASSERT( M==8 && D==3 );
  NASSERT_ALIGN32( y         );
  NASSERT_ALIGN32( x         );
  NASSERT_ALIGN32( coef      );
  NASSERT_ALIGN32( delayLine );

  Y = (      xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;
  S = (      xb_vecNx16*)delayLine;

  //
  // Load the delay line state.
  //

  x0 = BBE_LVNX16_I( S, 0*4*BBE_SIMD_WIDTH/2 );

  d0 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_4);
  d1 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_8);
  d2 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_14);

  BBE_SVNX16_I( d2, S, 0 );

  //
  // Process data.
  //

  cf = BBE_LVNX16_I( (const xb_vecNx16*)coef, 0 );

  __Pragma("ymemory(X)")
  __Pragma("ymemory(S)")
  for ( n=0; n<N/(BBE_SIMD_WIDTH/2); n++ )
  {
    // Load 8x3 input samples, CQ15
    BBE_LVNX16_IP( x0, X, 2 * BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x1, X, 2 * BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x2, X, 2 * BBE_SIMD_WIDTH );

    // Transposition 8x3 -> 3x8 (3 banks, each of 8 samples).
    BBE_DSELNX16I  ( s2, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_C3_STEP_0 );
    BBE_DSELNX16I_H( s2, x1, x1, x2, BBE_DSELI_DEINTERLEAVE_C3_STEP_1 );
    BBE_DSELNX16I  ( s1, s0, x1, x0, BBE_DSELI_DEINTERLEAVE_2        );

    d2 = BBE_LVNX16_I( S, 0 );

    BBE_SVNX16_I( s2, S, 0 );

    //
    // Coefficients bank 1, offsets 0 and 1
    //

    BBE_SELPCNX16I( p01, p00, s0, d0, 6 );

    q0 = BBE_EXTRNX16C( cf, 1 );

    w0 = BBE_MULNX16PR( p01, p00, q0 );

    //
    // Coefficients bank 2, offsets 0 and 1
    //

    BBE_SELPCNX16I( p11, p10, s1, d1, 6 );

    q1 = BBE_EXTRNX16C( cf, 2 );

    BBE_MULANX16PR( w0, p11, p10, q1 );

    //
    // Coefficients bank 0, offsets 0 and 1
    //

    BBE_SELPCNX16I( p21, p20, s2, d2, 5 );

    q2 = BBE_EXTRNX16C( cf, 0 );

    BBE_MULANX16PR( w0, p21, p20, q2 );

    //
    // Coefficient banks 0 and 1, offset 2
    //

    p02 = s0;                                              // 1
    p22 = BBE_SELNX16I( s2, d2, BBE_SELI_ROTATE_RIGHT_14 ); // 0

    q3 = BBE_EXTRNX16C( cf, 3 );

    BBE_MULANX16PR( w0, p02, p22, q3 );

    //
    // Update 2 delay line banks, the 3rd is updated through memory.
    //

    d0 = s0;
    d1 = s1;

    //
    // Save 16 output samples.
    //

    // CQ15 <- CQ30 - 15 w/ rounding and saturation.
    y0 = BBE_PACKQNX40( w0 );

    BBE_SVNX16_IP( y0, Y, 2 * BBE_SIMD_WIDTH );
  }

  //
  // Save the delay line state.
  //

  d2 = BBE_LVNX16_I(S, 0);

  x0 = BBE_SELNX16I(d2, d2, BBE_SELI_ROTATE_RIGHT_10);
  x0 = BBE_SELNX16I(x0, d1, BBE_SELI_ROTATE_LEFT_4);
  x0 = BBE_SELNX16I(x0, d0, BBE_SELI_ROTATE_LEFT_4);

  BBE_SVNX16_I(x0, S, 0);

} // firdec_proc_D3_M8()
const tFirFxdxns firdec_3d_8_8n  ={&firdec_alloc_d3_m8,firdec_proc_D3_M8};

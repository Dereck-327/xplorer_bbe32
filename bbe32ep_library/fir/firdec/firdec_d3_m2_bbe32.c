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

/* processing function D==3, M==2 */
#if !(HAVE_DIPACK && HAVE_MULPC && 1)
DISCARD_FUN(void, firdec_proc_D3_M2,(  int16_t * restrict y,
                         const int16_t * restrict x,
                         const int16_t * restrict coef,
                               int16_t * restrict delayLine,
                         int M, int N, int D ))
#else
void firdec_proc_D3_M2( int16_t * restrict y,
                 const int16_t * restrict x,
                 const int16_t * restrict coef,
                       int16_t * restrict delayLine,
                 int M, int N, int D )
{
        xb_vecNx16 * restrict Y;
  const xb_vecNx16 * X;

  valign X_va;

  xb_vecNx40 w0;

  xb_vecNx16 x0, x1, x2;
  xb_vecNx16 d0, d1, d2;
  xb_vecNx16 y0;

  xb_vecNx16 h0;

  int n;

  NASSERT( N>0 && !(N&7) );

  NASSERT( D==3 && M==BBE_SIMD_WIDTH );

  NASSERT_ALIGN32( y         );
  NASSERT_ALIGN32( x         );
  NASSERT_ALIGN32( coef      );
  NASSERT_ALIGN32( delayLine );

  Y = (      xb_vecNx16*)y;
  X = (const xb_vecNx16*)x;

  X_va = BBE_LANX16_PP( X );

  x0 = BBE_LVNX16_I( (const xb_vecNx16*)delayLine, 0 );
  BBE_LAVNX16_XP( x1, X_va, X, 4*7 );
  d0 = BBE_SELNX16I( x1, x0, BBE_SELI_ROTATE_LEFT_2 );

  h0 = BBE_LVNX16_I( (const xb_vecNx16*)coef, 0 );

  for ( n=0; n<N/(BBE_SIMD_WIDTH/2); n++ )
  {
    BBE_LAVNX16_XP( d1, X_va, X, 2*BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( d2, X_va, X, 2*BBE_SIMD_WIDTH );
  
    BBE_DSELNX16I( x1, x0, d1, d0, BBE_DSELI_DEINTERLEAVE_C3_STEP_0 );
    BBE_DSELNX16I_H( x1, x2, d1, d2, BBE_DSELI_DEINTERLEAVE_C3_STEP_1 );

    w0 = BBE_MULNX16PC_0( x0, h0 );
    BBE_MULNX16PC_1( w0, x2, h0 );

    // CQ15 <- CQ30 - 15 w/ rounding and saturation.
    y0 = BBE_DIPACKQNX40C(w0);

    BBE_SVNX16_IP( y0, Y, 2 * BBE_SIMD_WIDTH );

    BBE_LAVNX16_XP( d0, X_va, X, 2*BBE_SIMD_WIDTH );
  }

  d0 = BBE_SELNX16I( d0, d0, BBE_SELI_ROTATE_RIGHT_2 );
  BBE_SVNX16_I( d0, (xb_vecNx16*)delayLine, 0 );

}
#endif

const tFirFxdxns firdec_3d_2_8n  ={&firdec_alloc_d2_m2,firdec_proc_D3_M2};

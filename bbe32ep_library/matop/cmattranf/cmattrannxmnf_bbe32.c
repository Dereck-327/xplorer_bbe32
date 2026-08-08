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
DISCARD_FUN(void, cmattrannxmnf, ( complex_float * restrict y, 
                                   complex_float * restrict x, 
                                   int N, int M, int L ))
#else

#define sz_cf32 ((int)sizeof(complex_float))

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

/* Block Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void cmattrannxmnf ( complex_float * restrict y, 
                     complex_float * restrict x, 
                     int N, int M, int L )
{
    const xb_vecN_4xcf32 * restrict px0;
    const xb_vecN_4xcf32 * restrict px1;
          xb_vecN_4xcf32 * restrict py0;
          xb_vecN_4xcf32 * restrict py1;
    int strideX, strideY;
    int backX, backY;
    int l, p, q, P, Q;

    xb_vecN_4xcf32 X0, X1, X2, X3;
    xb_vecN_4xcf32 Y0, Y1, Y2, Y3;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);

    if (N<=0 || M<=0 || L<=0) return;

    /* At first select the direction of traverse of matrices */
    if (M>N)
    {
      /* traverse input matrix by columns, *
       * output - by rows                  */
      P = N;
      Q = M;
      strideX = 3*N*sz_cf32;
      strideY = -2*M*sz_cf32+4*sz_cf32;
      backX = 4*sz_cf32 - N*M*sz_cf32;
      backY = 3*M*sz_cf32;
    }
    else
    {
      /* traverse input matrix by rows, *
       * output - by columns            */
      P = M;
      Q = N;
      strideX = -N*sz_cf32+4*sz_cf32;
      strideY = 2*M*sz_cf32;
      backX = 3*N*sz_cf32;
      backY = 4*sz_cf32 - N*M*sz_cf32;
    }
    
    __Pragma("loop_count min=1");
    for ( l=0; l<L; l++ )
    {
        px0 = (const xb_vecN_4xcf32 *)(x+l*M*N);
        px1 = (const xb_vecN_4xcf32 *)(x+l*M*N+2*N);
        py0 = (      xb_vecN_4xcf32 *)(y+l*M*N);
        py1 = (      xb_vecN_4xcf32 *)(y+l*M*N+M);

        /* Transpose matrices by 4x4 blocks */
        __Pragma("loop_count min=1");
        for ( p=0; p<(P>>2); p++ )
        {
            __Pragma("loop_count min=1");
            for ( q=0; q<(Q>>2); q++ )
            {
                /* Load input matrix X */
                BBE_LVN_4XCF32_XP(X0, px0, N*sz_cf32);
                BBE_LVN_4XCF32_XP(X1, px0, strideX);
                BBE_LVN_4XCF32_XP(X2, px1, N*sz_cf32);
                BBE_LVN_4XCF32_XP(X3, px1, strideX);
                /* Transpose and conjugate */
                BBE_DSELN_4XCF32I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
                BBE_DSELN_4XCF32I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
                BBE_DSELN_4XCF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_4);
                BBE_DSELN_4XCF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_4);
                Y0 = BBE_CONJN_4XCF32(X0);
                Y1 = BBE_CONJN_4XCF32(X1);
                Y2 = BBE_CONJN_4XCF32(X2);
                Y3 = BBE_CONJN_4XCF32(X3);
                /* Save results */
                BBE_SVN_4XCF32_XP(Y0, py0, 2*M*sz_cf32);
                BBE_SVN_4XCF32_XP(Y1, py1, 2*M*sz_cf32);
                BBE_SVN_4XCF32_XP(Y2, py0, strideY);
                BBE_SVN_4XCF32_XP(Y3, py1, strideY);
            }
            /* Jump to the next block */
            px0 = (const xb_vecN_4xcf32 *)((intptr_t)px0+backX);
            px1 = (const xb_vecN_4xcf32 *)((intptr_t)px1+backX);
            py0 = (      xb_vecN_4xcf32 *)((intptr_t)py0+backY);
            py1 = (      xb_vecN_4xcf32 *)((intptr_t)py1+backY);
        }
    }
} /* cmattrannxmnf() */
#endif

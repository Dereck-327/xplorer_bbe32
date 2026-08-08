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
    Matrix Hermitian Product
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
DISCARD_FUN(void, cmathermnxmnf,( void * pScr,
                         complex_float * restrict y, 
                   const complex_float * restrict x, 
                   int N, int M, int L ))
#else
/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     M,N must be multiples of 4
*/
void cmathermnxmnf ( void * pScr,
                     complex_float * restrict y, 
               const complex_float * restrict x, 
               int N, int M, int L )
{
  int nrow, ncol, m, l;

  const xb_vecN_4xcf32 * restrict px_col;
  const xb_vecN_4xcf32 * restrict px_row;
        xb_vecN_4xcf32 * restrict py0;
        xb_vecN_4xcf32 * restrict py1;
        xb_vecN_4xcf32 * restrict py2;
        xb_vecN_4xcf32 * restrict py3;
  xb_vecN_4xcf32 x0, x1;
  xb_vecN_4xcf32 x_sel00, x_sel01;
  xb_vecN_4xcf32 x_sel10, x_sel11;
  xb_vecN_4xcf32 x_sel20, x_sel21;
  xb_vecN_4xcf32 x_sel30, x_sel31;
  xb_vecN_4xcf32 y0, y1, y2, y3;
  xb_vecN_4xcf32 tmp0, tmp1;

  /* check restrictions */
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(M % 4 == 0);
  NASSERT(N % 4 == 0);

  if ( N<=0 || L<=0 ) return;
  if ( M<=0 )
  {
      int k;
      xb_vecN_4xcf32 zero;

      py0 = (xb_vecN_4xcf32 *)y;
      zero = BBE_ZERON_4XCF32();
      __Pragma("loop_count min=2, factor=2");
      for (k=0; k<((N*N*L)>>(LOG2_BBE_SIMD_WIDTH-2)); k++)
      {
          BBE_SVN_4XCF32_IP(zero, py0, 2*BBE_SIMD_WIDTH);
      }

      return;
  }
  
  py0 = (xb_vecN_4xcf32 *)(y+0*N);
  py1 = (xb_vecN_4xcf32 *)(y+1*N);
  py2 = (xb_vecN_4xcf32 *)(y+2*N);
  py3 = (xb_vecN_4xcf32 *)(y+3*N);
  
  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      px_row = (const xb_vecN_4xcf32 *)(x+l*N*M);
      /* Compute matrices by 4x(BBE_SIMD_WIDTH/4) pieces */
      __Pragma("loop_count min=1");
      for (nrow=0; nrow<(N>>2); nrow++)
      {
          px_col = (const xb_vecN_4xcf32 *)(x+l*N*M);
          __Pragma("loop_count min=1");
          for (ncol=0; ncol<(N>>(LOG2_BBE_SIMD_WIDTH-2)); ncol++)
          {
              y0 = y1 = y2 = y3 = BBE_ZERON_4XCF32();

              __Pragma("loop_count min=2, factor=2");
              for (m=0; m<(M>>1); m++)
              {
                  /* Load 2x(BBE_SIMD_WIDTH/4) elements of matrix X */
                  BBE_LVN_4XCF32_XP(x0, px_col, sizeof(complex_float)*N);
                  BBE_LVN_4XCF32_XP(x1, px_col, sizeof(complex_float)*N);

                  /* Load 2x4 elements of matrix X */
                  BBE_LVN_4XCF32_XP(tmp0, px_row, sizeof(complex_float)*N);
                  BBE_LVN_4XCF32_XP(tmp1, px_row, sizeof(complex_float)*N);
                  x_sel00 = BBE_REPN_4XCF32(tmp0, 0);
                  x_sel10 = BBE_REPN_4XCF32(tmp0, 1);
                  x_sel20 = BBE_REPN_4XCF32(tmp0, 2);
                  x_sel30 = BBE_REPN_4XCF32(tmp0, 3);
                  x_sel01 = BBE_REPN_4XCF32(tmp1, 0);
                  x_sel11 = BBE_REPN_4XCF32(tmp1, 1);
                  x_sel21 = BBE_REPN_4XCF32(tmp1, 2);
                  x_sel31 = BBE_REPN_4XCF32(tmp1, 3);

                  /* Y[4x(BBE_SIMD_WIDTH/4)] = X[4x2]' * X[2x(BBE_SIMD_WIDTH/4)] */
                  BBE_MULJAN_4XCF32(y0, x0, x_sel00);
                  BBE_MULJAN_4XCF32(y1, x0, x_sel10);
                  BBE_MULJAN_4XCF32(y2, x0, x_sel20);
                  BBE_MULJAN_4XCF32(y3, x0, x_sel30);
                  BBE_MULJAN_4XCF32(y0, x1, x_sel01);
                  BBE_MULJAN_4XCF32(y1, x1, x_sel11);
                  BBE_MULJAN_4XCF32(y2, x1, x_sel21);
                  BBE_MULJAN_4XCF32(y3, x1, x_sel31);
              }
              /* Save outputs */
              BBE_SVN_4XCF32_IP(y0, py0, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(y1, py1, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(y2, py2, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(y3, py3, 2*BBE_SIMD_WIDTH);

              /* Prepare pointers for next (BBE_SIMD_WIDTH/4) columns */
              px_col = (const xb_vecN_4xcf32 *)((complex_float*)px_col-N*M+(BBE_SIMD_WIDTH/4));
              px_row = (const xb_vecN_4xcf32 *)((complex_float*)px_row-N*M);
          }
          /* Jump to next 4 rows */
          py0 = py3;
          py1 = (xb_vecN_4xcf32 *)((complex_float*)py0+N);
          py2 = (xb_vecN_4xcf32 *)((complex_float*)py1+N);
          py3 = (xb_vecN_4xcf32 *)((complex_float*)py2+N);

          px_row  = (const xb_vecN_4xcf32 *)((complex_float*)px_row+4);
      }
  }
} /* cmathermnxmnf() */
#endif

/* Return the scratch area size, in bytes. */
size_t cmathermnxmnf_getScratchSize ( int N, int M, int L )
{
  return 0;
} /* cmathermnxmnf_getScratchSize() */

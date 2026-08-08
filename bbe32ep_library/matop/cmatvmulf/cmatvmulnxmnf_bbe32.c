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
    Complex Matrix-Matrix/Matrix-Vector Multiply
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
DISCARD_FUN(void, cmatvmulnxmnf, ( void * pScr,
                          complex_float * restrict z, 
                    const complex_float * restrict x, 
                    const complex_float * restrict y, 
                    int N, int M, int L ))
#else

/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     N,M must be multiples of 4
*/
void cmatvmulnxmnf ( void * pScr,
                     complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int N, int M, int L )
{
  const xb_vecN_4xcf32 * restrict px0;
  const xb_vecN_4xcf32 * restrict px1;
  const xb_vecN_4xcf32 * restrict px2;
  const xb_vecN_4xcf32 * restrict px3;
  const xb_vecN_4xcf32 * restrict py0;
  const xb_vecN_4xcf32 * restrict py1;
        xb_vecN_4xcf32 * restrict pz;
  int l, m, n;

  xb_vecN_4xcf32 X0, X1, X2, X3;
  xb_vecN_4xcf32 Y0, Y1;
  xb_vecN_4xcf32 Z, Z0, Z1, Z2, Z3;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(M % 4 == 0);
  NASSERT(N % 4 == 0);

  if ( M<=0 || L<=0 ) return;
  if ( N<=0 )
  {
      xb_vecN_4xcf32 zero;

      pz = (xb_vecN_4xcf32 *)z;
      zero = BBE_ZERON_4XCF32();
      __Pragma("loop_count min=1");
      for (m=0; m<((M*L)>>(LOG2_BBE_SIMD_WIDTH-2)); m++)
      {
          BBE_SVN_4XCF32_IP(zero, pz, 2*BBE_SIMD_WIDTH);
      }

      return;
  }

  px0 = (const xb_vecN_4xcf32 *)(x+0*N);
  px1 = (const xb_vecN_4xcf32 *)(x+1*N);
  px2 = (const xb_vecN_4xcf32 *)(x+2*N);
  px3 = (const xb_vecN_4xcf32 *)(x+3*N);
  py0 = (const xb_vecN_4xcf32 *)(y);
  py1 = (const xb_vecN_4xcf32 *)(y);
  pz  = (      xb_vecN_4xcf32 *)(z);
  
  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      __Pragma("loop_count min=1");
      for (m=0; m<(M>>2); m++)
      {
          Z0 = Z1 = Z2 = Z3 = BBE_ZERON_4XCF32();
          __Pragma("loop_count min=1");
          for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-2)); n++)
          {
              /* Load part of input matrices X and Y */
              BBE_LVN_4XCF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
              BBE_LVN_4XCF32_IP(X1, px1, 2*BBE_SIMD_WIDTH);
              BBE_LVN_4XCF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
              BBE_LVN_4XCF32_IP(X3, px3, 2*BBE_SIMD_WIDTH);
              BBE_LVN_4XCF32_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
              BBE_LVN_4XCF32_IP(Y1, py1, 2*BBE_SIMD_WIDTH);

              /* Multiply input matrices X and Y */
              BBE_MULAN_4XCF32(Z0, X0, Y0);
              BBE_MULAN_4XCF32(Z1, X1, Y0);
              BBE_MULAN_4XCF32(Z2, X2, Y1);
              BBE_MULAN_4XCF32(Z3, X3, Y1);
          }

          BBE_DSELN_4XCF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELN_4XCF32I(Z3, Z1, Z3, Z1, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELN_4XCF32I(Z1, Z0, Z1, Z0, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELN_4XCF32I(Z3, Z2, Z3, Z2, BBE_DSELI_INTERLEAVE_4);

          Z0 = BBE_ADDN_4XCF32(Z0, Z1);
          Z2 = BBE_ADDN_4XCF32(Z2, Z3);
          Z  = BBE_ADDN_4XCF32(Z0, Z2);

          /* Save results */
          BBE_SVN_4XCF32_IP(Z, pz, 2*BBE_SIMD_WIDTH);

          /* Prepare pointers for next 4 rows */
          px0 = (const xb_vecN_4xcf32 *)((complex_float*)px0+3*N);
          px1 = (const xb_vecN_4xcf32 *)((complex_float*)px1+3*N);
          px2 = (const xb_vecN_4xcf32 *)((complex_float*)px2+3*N);
          px3 = (const xb_vecN_4xcf32 *)((complex_float*)px3+3*N);
          py0 = (const xb_vecN_4xcf32 *)((complex_float*)py0-N);
          py1 = (const xb_vecN_4xcf32 *)((complex_float*)py1-N);
      }
      py0 = (const xb_vecN_4xcf32 *)((complex_float*)py0+N);
      py1 = (const xb_vecN_4xcf32 *)((complex_float*)py1+N);
  }
} /* cmatvmulnxmnf() */
#endif

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmnf_getScratchSize ( int N, int M, int L )
{
  return 0;
} /* cmatvmulnxmnf_getScratchSize() */

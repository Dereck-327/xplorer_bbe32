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
DISCARD_FUN(void, cmatmulnxmnf,( void * pScr,
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

/* Block Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     N,M must be multiples of 4
*/
void cmatmulnxmnf ( void * pScr,
                    complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int N, int M, int L )
{
  int mrow, mcol, n, l;

  const complex_float  * restrict px0;
  const complex_float  * restrict px1;
  const complex_float  * restrict px2;
  const complex_float  * restrict px3;
  const xb_vecN_4xcf32 * restrict py;
        xb_vecN_4xcf32 * restrict pz0;
        xb_vecN_4xcf32 * restrict pz1;
        xb_vecN_4xcf32 * restrict pz2;
        xb_vecN_4xcf32 * restrict pz3;
  xb_vecN_4xcf32 x_sel00, x_sel01;
  xb_vecN_4xcf32 x_sel10, x_sel11;
  xb_vecN_4xcf32 x_sel20, x_sel21;
  xb_vecN_4xcf32 x_sel30, x_sel31;
  xb_vecN_4xcf32 y0, y1;
  xb_vecN_4xcf32 z0, z1, z2, z3;
  xb_vecN_4xcf32 tmp0, tmp1, tmp2, tmp3;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(M % 4 == 0);
  NASSERT(N % 4 == 0);

  if ( M<=0 || L<=0 ) return;
  if ( N<=0 )
  {
      int k;
      xb_vecN_4xcf32 zero;

      pz0 = (xb_vecN_4xcf32 *)z;
      zero = BBE_ZERON_4XCF32();
      __Pragma("loop_count min=2, factor=2");
      for (k=0; k<((M*M*L)>>(LOG2_BBE_SIMD_WIDTH-2)); k++)
      {
          BBE_SVN_4XCF32_IP(zero, pz0, 2*BBE_SIMD_WIDTH);
      }

      return;
  }
  
  pz0 = (xb_vecN_4xcf32 *)(z+0*M);
  pz1 = (xb_vecN_4xcf32 *)(z+1*M);
  pz2 = (xb_vecN_4xcf32 *)(z+2*M);
  pz3 = (xb_vecN_4xcf32 *)(z+3*M);
  px0 = (x+0*N);
  px1 = (x+1*N);
  px2 = (x+2*N);
  px3 = (x+3*N);
  
  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      /* Compute matrices by 4x(BBE_SIMD_WIDTH/4) pieces */
      __Pragma("loop_count min=1");
      for (mrow=0; mrow<(M>>2); mrow++)
      {
          py  = (const xb_vecN_4xcf32 *)(y+l*N*M);
          __Pragma("loop_count min=1");
          for (mcol=0; mcol<(M>>(LOG2_BBE_SIMD_WIDTH-2)); mcol++)
          {
              z0 = z1 = z2 = z3 = BBE_ZERON_4XCF32();

              __Pragma("loop_count min=2, factor=2");
              for (n=0; n<(N>>1); n++)
              {
                  /* Load 2x(BBE_SIMD_WIDTH/4) elements of matrix Y */
                  BBE_LVN_4XCF32_XP(y0, py, sizeof(complex_float)*M);
                  BBE_LVN_4XCF32_XP(y1, py, sizeof(complex_float)*M);

                  /* Load 4x4 elements of matrix X */
                  BBE_LSN_4XCF32_IP(tmp0, px0, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp1, px1, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp2, px2, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp3, px3, sizeof(complex_float));
                  x_sel00 = BBE_REPN_4XCF32(tmp0, 0);
                  x_sel10 = BBE_REPN_4XCF32(tmp1, 0);
                  x_sel20 = BBE_REPN_4XCF32(tmp2, 0);
                  x_sel30 = BBE_REPN_4XCF32(tmp3, 0);
                  BBE_LSN_4XCF32_IP(tmp0, px0, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp1, px1, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp2, px2, sizeof(complex_float));
                  BBE_LSN_4XCF32_IP(tmp3, px3, sizeof(complex_float));
                  x_sel01 = BBE_REPN_4XCF32(tmp0, 0);
                  x_sel11 = BBE_REPN_4XCF32(tmp1, 0);
                  x_sel21 = BBE_REPN_4XCF32(tmp2, 0);
                  x_sel31 = BBE_REPN_4XCF32(tmp3, 0);

                  /* Z[4x(BBE_SIMD_WIDTH/4)] = X[4x2]*Y[2x(BBE_SIMD_WIDTH/4)] */
                  BBE_MULAN_4XCF32(z0, x_sel00, y0);
                  BBE_MULAN_4XCF32(z0, x_sel01, y1);
                  BBE_MULAN_4XCF32(z1, x_sel10, y0);
                  BBE_MULAN_4XCF32(z1, x_sel11, y1);
                  BBE_MULAN_4XCF32(z2, x_sel20, y0);
                  BBE_MULAN_4XCF32(z2, x_sel21, y1);
                  BBE_MULAN_4XCF32(z3, x_sel30, y0);
                  BBE_MULAN_4XCF32(z3, x_sel31, y1);
              }
              /* Save outputs */
              BBE_SVN_4XCF32_IP(z0, pz0, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(z1, pz1, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(z2, pz2, 2*BBE_SIMD_WIDTH);
              BBE_SVN_4XCF32_IP(z3, pz3, 2*BBE_SIMD_WIDTH);

              /* Prepare pointers for next (BBE_SIMD_WIDTH/4) columns */
              px0 = px0-N;
              px1 = px1-N;
              px2 = px2-N;
              px3 = px3-N;
              py  = (const xb_vecN_4xcf32 *)((complex_float*)py-N*M+(BBE_SIMD_WIDTH/4));
          }
          /* Jump to next 4 rows */
          pz0 = pz3;
          pz1 = (xb_vecN_4xcf32 *)((complex_float*)pz0+M);
          pz2 = (xb_vecN_4xcf32 *)((complex_float*)pz1+M);
          pz3 = (xb_vecN_4xcf32 *)((complex_float*)pz2+M);
          
          px0 = px0+4*N;
          px1 = px1+4*N;
          px2 = px2+4*N;
          px3 = px3+4*N;

          py  = (const xb_vecN_4xcf32 *)((complex_float*)py-M);
      }
  }
} /* cmatmulnxmnf() */
#endif

/* Return the scratch area size, in bytes. */
size_t cmatmulnxmnf_getScratchSize ( int N, int M, int L )
{
  return 0;
} /* cmatmulnxmnf_getScratchSize() */

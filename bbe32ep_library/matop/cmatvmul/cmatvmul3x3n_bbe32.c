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
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatvmul3x3n,(complex_fract16 * restrict z, 
            const complex_fract16 * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
#else

/* pack 3 5-bit into 15 bit: a - low, c - high */
#define PCK555(a,b,c) (a)+((b)<<5)+((c)<<10)
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

/* Block Order, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void cmatvmul3x3n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  /* packed selections */
  static
  const int16_t ALIGN(32) aSel0[BBE_SIMD_WIDTH]= {
      PCK555(0,2,4),      PCK555(1,3,5),
      PCK555(6,8,10),     PCK555(7,9,11),
      PCK555(12,14,0),    PCK555(13,15,0),
      PCK555(0,0,0),      PCK555(0,0,0),
      PCK555(16,18,20),   PCK555(17,19,21),
      PCK555(22,24,26),   PCK555(23,25,27),
      PCK555(28,30,0),    PCK555(29,31,0),
      PCK555(0,0,0),      PCK555(0,0,0)
  };
  int l;
  xb_vecNx16 vTmp;
  xb_vecNx16 vX00, vX01, vX10, vX11,vY0,vZ; 
  xb_vecNx16 vXsel, vYsel;
  xb_vecNx40 wvZ;

  vselN vsX0, vsX1, vsX2;
  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;

  const vsaN  q= BBE_MOVVSA32(Q);

  /* check restrictions */
  NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));;
  NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));;
  NASSERT_ALIGN(z,(2*BBE_SIMD_WIDTH));;
  NASSERT(L%2==0);
  NASSERT(Q>=0 && Q<=16);
  if (L<=0) return;

  vTmp= BBE_LVNX16_I((const xb_vecNx16*)aSel0,0);
  vsX0= BBE_MOVVSELNX16(vTmp,0);
  vsX1= BBE_MOVVSELNX16(vTmp,5);
  vsX2= BBE_MOVVSELNX16(vTmp,10);

  __Pragma( "loop_count min=1" )
  for ( l=0; l<L; l+=2 )
  {
    /* load 2 X matrices */
    BBE_LVNX16_IP(vX00,pXrd,2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(vX01,pXrd,2*BBE_SIMD_WIDTH);

    BBE_LVNX16_IP(vX10,pXrd,2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(vX11,pXrd,2*BBE_SIMD_WIDTH);

    /* load 2 Y vectors */
    BBE_LVNX16_IP(vY0,pYrd,2*BBE_SIMD_WIDTH);

    vXsel= BBE_SELNX16(vX10, vX00, vsX0);
    vYsel= BBE_SHFLNX16I(vY0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_1);
    wvZ= BBE_MULRNX16C(vXsel, vYsel, q);

    vXsel= BBE_SELNX16(vX10, vX00, vsX1);
    vYsel= BBE_SHFLNX16I(vY0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_2);
    BBE_MULANX16C(wvZ, vXsel, vYsel);

    vXsel= BBE_SELNX16(vX10, vX00, vsX2);
    /* pack 8-th elements of X0 and X1 */
    vTmp= BBE_SELNX16I(vX11,vX01,BBE_SELI_PACK_8);

    vXsel= BBE_SELNX16I(vTmp, vXsel,BBE_SELI_INTERLEAVE_4_EVEN);
    vYsel= BBE_SHFLNX16I(vY0, BBE_SHFLI_MMC4X4X4X4_M1_STEP_3);
    BBE_MULANX16C(wvZ, vXsel, vYsel);

    /* pack and save results */
    vZ = BBE_PACKVNX40(wvZ, q);
    BBE_SVNX16_IP(vZ, pZwr, 2*BBE_SIMD_WIDTH);
  }
} /* cmatvmul3x3n() */
#endif

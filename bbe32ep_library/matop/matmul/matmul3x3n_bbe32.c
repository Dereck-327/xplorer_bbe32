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
    Real Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, matmul3x3n,(int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q))
#else
/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/      
// pack 3 5-bit into 15 bit: a - low, c - high
#define PCK555(a,b,c) (a)+((b)<<5)+((c)<<10)                                                                                                                                                
void matmul3x3n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q )
{
 // packed selections
  static
  const int16_t ALIGN(32) aSel0[BBE_SIMD_WIDTH]=  {   PCK555(0,0,0),      PCK555(3,4,1),
                                                      PCK555(6,8,2),      PCK555(15,2,15),
                                                      PCK555(1,6,3),      PCK555(4,10,4),
                                                      PCK555(7,16,5),     PCK555(15,20,15),
                                                      PCK555(2,24,6),     PCK555(5,31,7),
                                                      PCK555(8,31,8),     PCK555(1,31,15),
                                                      PCK555(15,31,15),   PCK555(15,31,15),
                                                      PCK555(15,31,15),   PCK555(15,31,15)};
  int l;
  vsaN  q = BBE_MOVVSA32(Q);

  xb_vecNx16 vTmp;

  vselN vsSelN0, vsSelN1, vsConv;

  const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
  const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
        xb_vecNx16 * restrict pZwr = (      xb_vecNx16 *)z;

  xb_vecNx16 vX, vY, vX0, vX1, vX2, vZ, vZ0, vZ1;

  xb_vecNx40 wvAcc0, wvAcc1;

  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));
  NASSERT(Q >= 0 && Q <= 16);

  if (L<=0) return;
  vTmp  = BBE_LVNX16_I((const xb_vecNx16*)aSel0,0);
  vsSelN0 = BBE_MOVVSV(vTmp, 0);
  vsSelN1 = BBE_MOVVSV(vTmp, 5);
  vsConv= BBE_MOVVSV(vTmp, 10);

  const vboolN bLoad= BBE_LTRNI(15);

  __Pragma( "loop_count min=1" )
  for (l=0; l<L; l++)
  {
      /* load entire X matrix */
      BBE_LVNX16T_IP(vX, pXrd, 2*BBE_SIMD_WIDTH,bLoad);
      /* transform 3x3 to 4x4 */
      vX= BBE_SHFLNX16(vX,vsConv);

      vX = BBE_CONJSNX16C(vX);

      vX0 = BBE_SHFLNX16I(vX, BBE_SHFLI_REP_0X4);
      vX1 = BBE_SHFLNX16I(vX, BBE_SHFLI_REP_1X4);
      vX2 = BBE_SHFLNX16I(vX, BBE_SHFLI_REP_2X4);
        
      /* load entire Y matrix */
      BBE_LVNX16_IP(vY, pYrd, 2*BBE_SIMD_WIDTH);
      vY = BBE_SHFLNX16(vY, vsSelN0);

      wvAcc0 = BBE_MULRNX16PC_0(vY, vX0, q);
      BBE_MULANX16PC_1(wvAcc0, vY, vX1);
      wvAcc1 = BBE_MULRNX16PC_0(vY, vX2, q);

      vZ0 = BBE_PACKVNX40(wvAcc0, q); 
      vZ1 = BBE_PACKVNX40(wvAcc1, q);
      vZ = BBE_SELNX16(vZ1, vZ0, vsSelN1);

      BBE_SVNX16_IP(vZ, pZwr, 2*BBE_SIMD_WIDTH);
  }
} /* matmul3x3n() */
#endif

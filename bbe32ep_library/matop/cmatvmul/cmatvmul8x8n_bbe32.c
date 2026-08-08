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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatvmul8x8n,(complex_fract16 * restrict z, 
            const complex_fract16 * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
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

/* Block Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     None
*/
void cmatvmul8x8n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q )
{
  /*
  final permutations for mulpc_0:
  x00x01 x20x21 x40x41 x60x61  y0y1y0y1y0y1y0y1
  x02x03 x22x23 x42x43 x62x63  y2y3y2y3y2y3y2y3
  x04x05 x24x25 x44x45 x64x65  y4y5y4y5y4y5y4y5
  x06x07 x26x27 x46x47 x66x67  y6y7y6y7y6y7y6y7
  */
  int l;

  const xb_vecNx16 *px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *py = (const xb_vecNx16 *)y;
  xb_vecNx16 * restrict pz = (xb_vecNx16 *)z;
  xb_vecNx16 X0,X1,X2,X3,X4,X5,X6,X7, Y0,  XX,YY,Z;
  xb_vecNx40 r;
  vsaN  q = BBE_MOVVSA32(Q);

  /* check restrictions */
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L <= 0) return;

  __Pragma( "loop_count min=1" ); 
  for ( l=0; l<L; l++ )
  {
    /* Load input matrix X and vector Y */
    BBE_LVNX16_IP(X0, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X3, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X4, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X5, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X6, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X7, px, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y0, py, 2*BBE_SIMD_WIDTH);
        
    /*
    X0=x00x01 x20x21 x02x03 x22x23
    X2=x04x05 x24x25 x06x07 x26x27 
    X4=x40x41 x60x61 x42x43 x62x63
    X6=x44x45 x64x65 x46x47 x66x67 
    X1=x10x11 x30x31 x12x13 x32x33
    X3=x14x15 x34x35 x16x17 x36x37 
    X5=x50x51 x70x71 x52x53 x72x73
    X7=x54x55 x74x75 x56x57 x76x77 
    */
    BBE_DSELNX16I(X2,X0,X2,X0,BBE_DSELI_INTERLEAVE_4); 
    BBE_DSELNX16I(X6,X4,X6,X4,BBE_DSELI_INTERLEAVE_4); 
    BBE_DSELNX16I(X3,X1,X3,X1,BBE_DSELI_INTERLEAVE_4); 
    BBE_DSELNX16I(X7,X5,X7,X5,BBE_DSELI_INTERLEAVE_4); 
    /* multiply:
     x00x01 x20x21 x40x41 x60x61  y0y1y0y1y0y1y0y1 */
    XX=BBE_SELNX16I(X4,X0,BBE_SELI_EXTRACT_LO_HALVES);
    YY=BBE_SHFLNX16I(Y0,BBE_SHFLI_REP_0X4);
    r=BBE_MULRNX16PC_0(XX,YY,q);
    XX=BBE_SELNX16I(X5,X1,BBE_SELI_EXTRACT_LO_HALVES);
    BBE_MULANX16PC_1(r,XX,YY);
    /* x02x03 x22x23 x42x43 x62x63  y2y3y2y3y2y3y2y3 */
    XX=BBE_SELNX16I(X4,X0,BBE_SELI_EXTRACT_HI_HALVES);
    YY=BBE_SHFLNX16I(Y0,BBE_SHFLI_REP_1X4);
    BBE_MULANX16PC_0(r,XX,YY);
    XX=BBE_SELNX16I(X5,X1,BBE_SELI_EXTRACT_HI_HALVES);
    BBE_MULANX16PC_1(r,XX,YY);
    /* x04x05 x24x25 x44x45 x64x65  y4y5y4y5y4y5y4y5 */
    XX=BBE_SELNX16I(X6,X2,BBE_SELI_EXTRACT_LO_HALVES);
    YY=BBE_SHFLNX16I(Y0,BBE_SHFLI_REP_2X4);
    BBE_MULANX16PC_0(r,XX,YY);
    XX=BBE_SELNX16I(X7,X3,BBE_SELI_EXTRACT_LO_HALVES);
    BBE_MULANX16PC_1(r,XX,YY);
    /* x06x07 x26x27 x46x47 x66x67  y6y7y6y7y6y7y6y7 */
    XX=BBE_SELNX16I(X6,X2,BBE_SELI_EXTRACT_HI_HALVES);
    YY=BBE_SHFLNX16I(Y0,BBE_SHFLI_REP_3X4);
    BBE_MULANX16PC_0(r,XX,YY);
    XX=BBE_SELNX16I(X7,X3,BBE_SELI_EXTRACT_HI_HALVES);
    BBE_MULANX16PC_1(r,XX,YY);
    /* packing */
    Z=BBE_PACKVNX40(r,q);
    BBE_SVNX16_IP(Z, pz, 2*BBE_SIMD_WIDTH);
  }
} /* cmatvmul8x8n() */
#endif
